library(rbcb)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(seasonal)
library(ggplot2)

inicio <- as.Date("2001-01-01")
fim    <- as.Date("2025-12-31")

# Série diária do SGS em blocos (a API recusa janelas longas para séries diárias)
baixar_sgs_diaria <- function(codigo, nome, inicio, fim, anos_bloco = 2) {
  cortes <- seq(inicio, fim, by = paste(anos_bloco, "years"))
  fins   <- c(cortes[-1] - 1, fim)
  Map(function(a, b) get_series(setNames(codigo, nome), start_date = a, end_date = b),
      cortes, fins) |>
    bind_rows() |>
    distinct(date, .keep_all = TRUE) |>
    arrange(date)
}

# Meta relevante para a expectativa de inflação 12 meses à frente.
# A expectativa Focus do mês m cobre os meses m+1 a m+12: (12 − m) meses no ano
# corrente e m meses no ano seguinte. Média dos três meses do trimestre.
meta_no_horizonte <- function(trimestres, metas) {
  metas_seg <- metas |> transmute(ano = ano - 1, meta_seguinte = meta)
  tidyr::expand_grid(trimestre = trimestres, k = 0:2) |>
    mutate(mes_data = trimestre %m+% months(k),
           ano      = year(mes_data),
           f        = (12 - month(mes_data)) / 12) |>
    left_join(metas,     by = "ano") |>
    left_join(metas_seg, by = "ano") |>
    mutate(meta_mes = f * meta + (1 - f) * meta_seguinte) |>
    group_by(trimestre) |>
    summarise(meta_12m = mean(meta_mes), .groups = "drop")
}

# Metas de inflação (CMN). Inclui 2026: a expectativa 12 meses de 2025 já cobre 2026.
metas_anuais <- tribble(
  ~ano,  ~meta,
  2001,  4.0,  2002,  3.5,  2003,  4.0,  2004,  5.5,  2005,  4.5,  2006,  4.5,
  2007,  4.5,  2008,  4.5,  2009,  4.5,  2010,  4.5,  2011,  4.5,  2012,  4.5,
  2013,  4.5,  2014,  4.5,  2015,  4.5,  2016,  4.5,  2017,  4.5,  2018,  4.5,
  2019,  4.25, 2020,  4.0,  2021,  3.75, 2022,  3.5,  2023,  3.25, 2024,  3.0,
  2025,  3.0,  2026,  3.0
)

trimestres <- seq(floor_date(inicio, "quarter"), floor_date(fim, "quarter"), by = "quarter")

metas_tri <- tibble(trimestre = trimestres, ano = year(trimestres)) |>
  left_join(metas_anuais, by = "ano") |>
  left_join(meta_no_horizonte(trimestres, metas_anuais), by = "trimestre") |>
  select(trimestre, meta, meta_12m)

# Selic: meta (SGS 432), média diária do trimestre
selic_tri <- baixar_sgs_diaria(432, "meta", inicio, fim) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(selic = mean(meta, na.rm = TRUE), .groups = "drop")

print(ggplot(selic_tri, aes(trimestre, selic)) +
        geom_line(color = "darkblue", linewidth = 0.8) +
        labs(title = "Meta para a taxa Selic — média trimestral", y = "% a.a.", x = NULL) +
        theme_minimal())

# IPCA cheio e livres: ajuste sazonal mensal (X-13) e composição trimestral
ipca_raw <- get_series(c(ipca_cheio = 433, ipca_livres = 11428),
                       start_date = inicio, end_date = fim) |>
  reduce(full_join, by = "date") |>
  arrange(date)

inicio_ts <- c(year(min(ipca_raw$date)), month(min(ipca_raw$date)))
dessaz_mensal <- function(x) {
  as.numeric(final(seas(ts(x, start = inicio_ts, frequency = 12), transform.function = "none")))
}

ipca_tri <- ipca_raw |>
  mutate(ipca_cheio_sa  = dessaz_mensal(ipca_cheio),
         ipca_livres_sa = dessaz_mensal(ipca_livres),
         trimestre      = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(ipca_cheio_orig  = (prod(1 + ipca_cheio / 100) - 1) * 100,
            ipca_livres_orig = (prod(1 + ipca_livres / 100) - 1) * 100,
            ipca_cheio_sa    = (prod(1 + ipca_cheio_sa / 100) - 1) * 100,
            ipca_livres_sa   = (prod(1 + ipca_livres_sa / 100) - 1) * 100,
            n_meses = n(), .groups = "drop") |>
  filter(n_meses == 3) |>
  select(-n_meses)

# Conferência do ajuste sazonal: a série ajustada deve perder o padrão de cada trimestre
for (s in c("cheio", "livres")) print(
  ipca_tri |>
    select(trimestre, Original = paste0("ipca_", s, "_orig"), Dessazonalizado = paste0("ipca_", s, "_sa")) |>
    pivot_longer(-trimestre, names_to = "tipo", values_to = "valor") |>
    ggplot(aes(trimestre, valor, colour = tipo)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = c(Original = "gray60", Dessazonalizado = "steelblue")) +
    labs(title = paste("IPCA", s, "— trimestralização após dessazonalização mensal"),
         y = "% no trimestre", x = NULL, color = NULL) +
    theme_minimal())
ipca_tri <- ipca_tri |> select(trimestre, ipca_cheio_sa, ipca_livres_sa)

# Câmbio (SGS 3698): variação logarítmica da média trimestral
cambio_tri <- get_series(c(cambio = 3698), start_date = inicio, end_date = fim) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(cambio = mean(cambio, na.rm = TRUE), .groups = "drop") |>
  arrange(trimestre) |>
  mutate(var_cambio = (log(cambio) - lag(log(cambio))) * 100)

print(ggplot(cambio_tri, aes(trimestre, cambio)) +
        geom_line(color = "darkgreen", linewidth = 0.8) +
        labs(title = "Câmbio R$/US$ — média trimestral", y = "R$/US$", x = NULL) +
        theme_minimal())
print(cambio_tri |>
        left_join(metas_tri, by = "trimestre") |>
        filter(!is.na(var_cambio)) |>
        ggplot(aes(trimestre)) +
        geom_line(aes(y = var_cambio, color = "Variação trimestral"), linewidth = 0.6) +
        geom_line(aes(y = (meta - 2) / 4, color = "Variação de longo prazo (PPC)"), linewidth = 0.9) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        scale_color_manual(values = c("Variação trimestral" = "darkorange",
                                      "Variação de longo prazo (PPC)" = "darkblue")) +
        labs(title = "Câmbio: variação trimestral e referência da PPC, (meta − 2)/4",
             y = "% no trimestre", x = NULL, color = NULL) +
        theme_minimal() + theme(legend.position = "bottom"))
cambio_tri <- cambio_tri |> select(trimestre, var_cambio)

# IC-Br agregado (SGS 27574): variação logarítmica menos a meta trimestralizada
icbr_tri <- get_series(c(icbr = 27574), start_date = inicio, end_date = fim) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(icbr = mean(icbr, na.rm = TRUE), .groups = "drop") |>
  arrange(trimestre) |>
  left_join(metas_tri |> select(trimestre, meta), by = "trimestre") |>
  mutate(icbr_desvio = (log(icbr) - log(lag(icbr))) * 100 - meta / 4) |>
  select(trimestre, icbr_desvio)

print(ggplot(icbr_tri |> filter(!is.na(icbr_desvio)), aes(trimestre, icbr_desvio)) +
        geom_line(color = "darkblue", linewidth = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        labs(title = "Inflação importada: variação do IC-Br menos a meta trimestral",
             y = "p.p. no trimestre", x = NULL) +
        theme_minimal())

print(ggplot(metas_tri, aes(trimestre)) +
        geom_step(aes(y = meta, color = "Meta do ano-calendário"), linewidth = 0.8) +
        geom_line(aes(y = meta_12m, color = "Meta nos 12 meses à frente"), linewidth = 0.8) +
        scale_color_manual(values = c("Meta do ano-calendário" = "grey50",
                                      "Meta nos 12 meses à frente" = "darkred")) +
        labs(title = "Meta de inflação", y = "% a.a.", x = NULL, color = NULL) +
        theme_minimal() + theme(legend.position = "bottom"))

# Consolidação e checagem
dados1 <- metas_tri |>
  full_join(selic_tri,  by = "trimestre") |>
  full_join(ipca_tri,   by = "trimestre") |>
  full_join(cambio_tri, by = "trimestre") |>
  full_join(icbr_tri,   by = "trimestre") |>
  arrange(trimestre) |>
  filter(trimestre >= as.Date("2002-01-01"), trimestre <= fim)

janela <- dados1 |> filter(trimestre >= as.Date("2002-07-01"), trimestre <= as.Date("2023-10-01"))
na_janela <- colSums(is.na(janela))
if (any(na_janela > 0)) warning("NAs na janela 2002T3–2023T4: ",
                                paste(names(na_janela)[na_janela > 0], collapse = ", "))

write.csv(dados1, "dados1.csv", row.names = FALSE)
cat("Gravado: dados1.csv\n")