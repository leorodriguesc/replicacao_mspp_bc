library(rbcb)
library(mFilter)
library(dplyr)
library(lubridate)
library(seasonal)
library(ggplot2)
library(tidyr)

inicio    <- as.Date("2002-01-01")
fim       <- as.Date("2025-12-31")
media_ini <- as.Date("2003-10-01")   # 2003T4: início da amostra do boxe
media_fim <- as.Date("2019-10-01")   # 2019T4: fim do pré-Covid

# Nuci: "auto" decide pela amplitude sazonal; TRUE/FALSE força a escolha
NUCI_SA    <- "auto"
LIMIAR_SAZ <- 0.5     # p.p.; abaixo disso a série já é tratada como dessazonalizada

# X-13 trimestral, com fallback para X-11; devolve a série ajustada e a amplitude sazonal
dessaz_tri <- function(x, inicio_serie, nome = "", ...) {
  x_ts <- ts(x, start = c(year(inicio_serie), quarter(inicio_serie)), frequency = 4)
  m <- tryCatch(seas(x_ts, ...), error = function(e) {
    message(nome, ": SEATS falhou (", conditionMessage(e), "); usando X-11")
    seas(x_ts, x11 = "", ...)
  })
  sa  <- as.numeric(final(m))
  fat <- tapply(as.numeric(x_ts) - sa, as.numeric(cycle(x_ts)), mean)
  cat(sprintf("%-16s amplitude do padrão sazonal: %.2f\n", nome, diff(range(fat))))
  list(sa = sa, amplitude = diff(range(fat)))
}

media_janela <- function(x, trimestre, ini, fim) mean(x[trimestre >= ini & trimestre <= fim], na.rm = TRUE)

# Série de fim de trimestre (mar, jun, set, dez) a partir de uma série mensal do SGS
fim_de_trimestre <- function(df, coluna) {
  df |>
    filter(month(date) %in% c(3, 6, 9, 12)) |>
    mutate(trimestre = floor_date(date, "quarter")) |>
    select(trimestre, valor = all_of(coluna)) |>
    filter(trimestre >= inicio, trimestre <= fim) |>
    arrange(trimestre)
}

# Desocupação: PNAD Contínua (trimestre móvel datado pelo último mês)
desoc_tri <- get_series(c(pnad = 24369), start_date = inicio, end_date = fim) |>
  fim_de_trimestre("pnad")
desoc_tri$desoc <- dessaz_tri(desoc_tri$valor, min(desoc_tri$trimestre), "PNAD-C",
                              transform.function = "none")$sa
media_desoc <- media_janela(desoc_tri$desoc, desoc_tri$trimestre, min(desoc_tri$trimestre), media_fim)

print(desoc_tri |>
        pivot_longer(c(valor, desoc), names_to = "serie", values_to = "taxa") |>
        ggplot(aes(trimestre, taxa, color = serie)) +
        geom_line(linewidth = 0.8) +
        geom_hline(yintercept = media_desoc, linetype = "dashed", alpha = 0.5) +
        scale_color_manual(values = c(valor = "grey60", desoc = "darkblue"),
                           labels = c(valor = "Bruta", desoc = "Dessazonalizada")) +
        labs(title = "Taxa de desocupação — PNAD Contínua",
             subtitle = "Tracejado: média 2012T1–2019T4 (o zero de femp)",
             y = "% da força de trabalho", x = NULL, color = NULL) +
        theme_minimal())

desoc_tri <- desoc_tri |> transmute(trimestre, femp = -(desoc - media_desoc))

# Nuci
nuci_tri <- get_series(c(nuci = 28561), start_date = inicio, end_date = fim) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(nuci_bruta = mean(nuci), n = n(), .groups = "drop") |>
  filter(n == 3, trimestre >= inicio, trimestre <= fim) |>
  select(-n) |>
  arrange(trimestre)

ds_nuci <- dessaz_tri(nuci_tri$nuci_bruta, min(nuci_tri$trimestre), "Nuci", transform.function = "none")
usar_sa_nuci <- if (identical(NUCI_SA, "auto")) ds_nuci$amplitude > LIMIAR_SAZ else isTRUE(NUCI_SA)
nuci_tri$nuci <- if (usar_sa_nuci) ds_nuci$sa else nuci_tri$nuci_bruta
cat("Nuci:", if (usar_sa_nuci) "série dessazonalizada" else "série original (sem sazonalidade relevante)", "\n")

print(nuci_tri |>
        mutate(ajustada = ds_nuci$sa) |>
        pivot_longer(c(nuci_bruta, ajustada), names_to = "serie", values_to = "valor") |>
        ggplot(aes(trimestre, valor, color = serie)) +
        geom_line(linewidth = 0.8) +
        scale_color_manual(values = c(nuci_bruta = "grey60", ajustada = "darkgreen"),
                           labels = c(nuci_bruta = "Original", ajustada = "Com X-13")) +
        labs(title = "Nuci — média trimestral",
             subtitle = paste("Usada no modelo:", if (usar_sa_nuci) "com X-13" else "original"),
             y = "%", x = NULL, color = NULL) +
        theme_minimal())
nuci_tri <- nuci_tri |>
  transmute(trimestre, fnuci = nuci - media_janela(nuci, trimestre, media_ini, media_fim))

# Caged: estoque no fim do trimestre, X-13 sobre 100·log e filtro HP
caged_tri <- get_series(c(caged = 28763), start_date = inicio, end_date = fim) |>
  fim_de_trimestre("caged")
caged_sa <- dessaz_tri(100 * log(caged_tri$valor), min(caged_tri$trimestre), "Caged (100·log)",
                       transform.function = "none")$sa
caged_hp <- hpfilter(caged_sa, freq = 1600, type = "lambda")

print(tibble(trimestre = caged_tri$trimestre, sa = caged_sa, tendencia = as.numeric(caged_hp$trend)) |>
        pivot_longer(-trimestre, names_to = "serie", values_to = "valor") |>
        ggplot(aes(trimestre, valor, color = serie)) +
        geom_line(linewidth = 0.8) +
        scale_color_manual(values = c(sa = "darkorange", tendencia = "darkred"),
                           labels = c(sa = "100·log do estoque, dessazonalizado", tendencia = "Tendência HP")) +
        labs(title = "Caged — estoque de empregos formais",
             subtitle = "Um salto em 2020T1, antes da pandemia, indicaria quebra do Novo Caged",
             y = "100·log", x = NULL, color = NULL) +
        theme_minimal() + theme(legend.position = "bottom"))

caged_tri <- caged_tri |> transmute(trimestre, fcaged = as.numeric(caged_hp$cycle))

# Resultado primário: NFSP com sinal invertido (superávit positivo)
rp_tri <- get_series(c(rp = 5783), start_date = inicio, end_date = fim) |>
  mutate(rp = -rp) |>
  fim_de_trimestre("rp")
rp_hp  <- hpfilter(rp_tri$valor, freq = 1600, type = "lambda")

print(rp_tri |>
        mutate(tendencia = as.numeric(rp_hp$trend)) |>
        pivot_longer(c(valor, tendencia), names_to = "serie", values_to = "pib") |>
        ggplot(aes(trimestre, pib, color = serie)) +
        geom_line(linewidth = 0.8) +
        scale_color_manual(values = c(valor = "darkblue", tendencia = "darkred"),
                           labels = c(valor = "Resultado primário (12 meses)", tendencia = "Tendência HP")) +
        labs(title = "Resultado primário — governo federal e BC (superávit positivo)",
             y = "% do PIB", x = NULL, color = NULL) +
        theme_minimal() + theme(legend.position = "bottom"))

rp_tri <- rp_tri |> transmute(trimestre, rp_hiato = valor - as.numeric(rp_hp$trend))

# Consolidação e checagens
dados2b <- desoc_tri |>
  full_join(nuci_tri,  by = "trimestre") |>
  full_join(caged_tri, by = "trimestre") |>
  full_join(rp_tri,    by = "trimestre") |>
  filter(trimestre >= inicio, trimestre <= fim) |>
  arrange(trimestre)

print(dados2b |>
        pivot_longer(-trimestre, names_to = "serie", values_to = "valor") |>
        ggplot(aes(trimestre, valor)) +
        geom_line(color = "darkblue", linewidth = 0.7, na.rm = TRUE) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        facet_wrap(~ serie, scales = "free_y", ncol = 2) +
        labs(title = "Componentes cíclicos que entram no modelo", y = NULL, x = NULL) +
        theme_minimal())

janela <- dados2b |> filter(trimestre >= as.Date("2002-07-01"), trimestre <= as.Date("2023-10-01"))
cat("NAs na janela 2002T3–2023T4 (esperado: só femp, antes de 2012):\n")
print(colSums(is.na(janela[, -1])))
cat("Correlações entre os indicadores pró-cíclicos (devem ser positivas):\n")
print(round(cor(janela[, c("femp", "fnuci", "fcaged")], use = "pairwise.complete.obs"), 2))

write.csv(dados2b, "dados2b.csv", row.names = FALSE)
cat("Gravado: dados2b.csv\n")