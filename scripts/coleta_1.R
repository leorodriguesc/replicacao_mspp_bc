# ============================================================
# BLOCO 1 — Selic, IPCA, câmbio e IC-Br (versão 2)
# ============================================================
# MUDANÇAS EM RELAÇÃO À VERSÃO ANTERIOR
#  [1] Selic: `selic` passa a ser a META (SGS 432, média diária do trimestre),
#      como no boxe (i_t = meta para a taxa Selic). A taxa efetiva (SGS 4189)
#      fica em `selic_efetiva`, só para comparação.
#  [2] Câmbio: novas colunas da PPC (nota 5 e eq. 4.2 do boxe)
#        dcambio_ppc    = (meta − 2) / 4           variação de longo prazo
#        dcambio_desvio = var_cambio − dcambio_ppc  Δê (entra na Phillips em t−1)
#  [3] Metas: tabela trimestral única, com
#        meta             meta do ano-calendário        (IC-Br, PPC)
#        meta_trimestral  meta / 4                      (mesma convenção do boxe)
#        meta_12m         meta média dos 12 meses cobertos pela expectativa
#                         Focus (regra de Taylor e eq. 5 de expectativas)
#  [4] IPCA: as colunas dessazonalizadas entram uma vez só (a versão anterior
#      gerava ipca_livres_sa.x/.y e ipca_cheio_sa.x/.y no dados1).
#  [5] IC-Br: desvio em relação a meta / 4 (antes: meta composta trimestral).
# ============================================================

library(rbcb)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(purrr)
library(seasonal)

inicio <- as.Date("2001-01-01")
fim    <- as.Date("2025-12-31")

# ------------------------------------------------------------
# Funções
# ------------------------------------------------------------

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
#  metas: ano, meta
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

# ------------------------------------------------------------
# METAS DE INFLAÇÃO (CMN)
# ------------------------------------------------------------
# Inclui 2026: a expectativa 12 meses de 2025 já cobre meses de 2026
metas_anuais <- tribble(
  ~ano,  ~meta,
  2001,  4.0,
  2002,  3.5,
  2003,  4.0,
  2004,  5.5,
  2005,  4.5,
  2006,  4.5,
  2007,  4.5,
  2008,  4.5,
  2009,  4.5,
  2010,  4.5,
  2011,  4.5,
  2012,  4.5,
  2013,  4.5,
  2014,  4.5,
  2015,  4.5,
  2016,  4.5,
  2017,  4.5,
  2018,  4.5,
  2019,  4.25,
  2020,  4.0,
  2021,  3.75,
  2022,  3.5,
  2023,  3.25,
  2024,  3.0,
  2025,  3.0,
  2026,  3.0
)

trimestres <- seq(floor_date(inicio, "quarter"), floor_date(fim, "quarter"), by = "quarter")

metas_tri <- tibble(trimestre = trimestres, ano = year(trimestres)) |>
  left_join(metas_anuais, by = "ano") |>
  mutate(meta_trimestral = meta / 4,
         dcambio_ppc     = (meta - 2) / 4) |>            # eq. 4.2: (π^meta − π*ss)/4, π*ss = 2
  left_join(meta_no_horizonte(trimestres, metas_anuais), by = "trimestre") |>
  select(-ano)

#SELIC
# ---- Meta (SGS 432): média diária no trimestre ----
meta_diaria <- baixar_sgs_diaria(
  432, "meta",
  as.Date("2001-01-01"),
  fim,
  anos_bloco = 2
)

selic_tri <- meta_diaria |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(selic = mean(meta, na.rm = TRUE), .groups = "drop")

# ---- Efetiva (SGS 4189), só para comparação ----
selic_efetiva_tri <- get_series(c(selic_efetiva = 4189), start_date = inicio, end_date = fim) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(selic_efetiva = mean(selic_efetiva, na.rm = TRUE), .groups = "drop")

selic_tri <- selic_tri |> left_join(selic_efetiva_tri, by = "trimestre")

dif_selic <- with(selic_tri, selic_efetiva - selic)
cat(sprintf("\nSelic efetiva − meta (média trimestral): média %+.2f p.p. | mín %+.2f | máx %+.2f\n",
            mean(dif_selic, na.rm = TRUE), min(dif_selic, na.rm = TRUE), max(dif_selic, na.rm = TRUE)))
cat("Esperado: em torno de −0,10 p.p. (a Selic over roda logo abaixo da meta).\n")

ggplot(selic_tri, aes(x = trimestre)) +
  geom_line(aes(y = selic, color = "Meta (usada no modelo)"), linewidth = 0.8) +
  geom_line(aes(y = selic_efetiva, color = "Efetiva"), linewidth = 0.6, linetype = "dashed") +
  scale_color_manual(values = c("Meta (usada no modelo)" = "darkblue", "Efetiva" = "grey50")) +
  labs(title = "Selic — média trimestral",
       y = "% a.a.", x = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

#IPCA
ipca_raw <- get_series(
  c(ipca_cheio = 433, ipca_livres = 11428),
  start_date = inicio,
  end_date = fim
) |>
  reduce(full_join, by = "date") |>
  arrange(date)

ano_inicio <- year(min(ipca_raw$date))
mes_inicio <- month(min(ipca_raw$date))

ts_cheio_mensal <- ts(ipca_raw$ipca_cheio, start = c(ano_inicio, mes_inicio), frequency = 12)
ts_livres_mensal <- ts(ipca_raw$ipca_livres, start = c(ano_inicio, mes_inicio), frequency = 12)

#dessazonalizando
x13_cheio_mensal <- seas(ts_cheio_mensal, transform.function = "none")
x13_livres_mensal <- seas(ts_livres_mensal, transform.function = "none")

ipca_raw <- ipca_raw |>
  mutate(
    ipca_cheio_sa = as.numeric(final(x13_cheio_mensal)),
    ipca_livres_sa = as.numeric(final(x13_livres_mensal))
  )

#trimestraliazar
ipca_tri <- ipca_raw |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(
    ipca_cheio_orig  = (prod(1 + ipca_cheio/100) - 1) * 100,
    ipca_livres_orig = (prod(1 + ipca_livres/100) - 1) * 100,
    ipca_cheio_sa    = (prod(1 + ipca_cheio_sa/100) - 1) * 100,
    ipca_livres_sa   = (prod(1 + ipca_livres_sa/100) - 1) * 100,
    
    n_meses = n(),
    .groups = "drop"
  ) |>
  filter(n_meses == 3) |>
  select(-n_meses) |>
  arrange(trimestre) |>
  mutate(ipca_cheio_ma = (lag(ipca_cheio_sa, 1) + lag(ipca_cheio_sa, 2) +    # termo α1I da eq. 1
                            lag(ipca_cheio_sa, 3) + lag(ipca_cheio_sa, 4)) / 4)

plot_livres <- ipca_tri |>
  select(trimestre, ipca_livres_orig, ipca_livres_sa) |>
  pivot_longer(cols = c(ipca_livres_orig, ipca_livres_sa), names_to = "tipo", values_to = "valor")

ggplot(plot_livres, aes(x = trimestre, y = valor, colour = tipo)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("ipca_livres_orig" = "gray60", "ipca_livres_sa" = "steelblue"),
    labels = c("Original", "Dessazonalizado")
  ) +
  labs(title = "IPCA Livres: Trimestralização após Dessazonalização Mensal",
       y = "% no trimestre", x = NULL, color = NULL) +
  theme_minimal()

plot_cheio <- ipca_tri |>
  select(trimestre, ipca_cheio_orig, ipca_cheio_sa) |>
  pivot_longer(cols = c(ipca_cheio_orig, ipca_cheio_sa), names_to = "tipo", values_to = "valor")

ggplot(plot_cheio, aes(x = trimestre, y = valor, colour = tipo)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    values = c("ipca_cheio_orig" = "gray60", "ipca_cheio_sa" = "darkred"),
    labels = c("Original", "Dessazonalizado")
  ) +
  labs(title = "IPCA Cheio: Trimestralização após Dessazonalização Mensal",
       y = "% no trimestre", x = NULL, color = NULL) +
  theme_minimal()

#CAMBIO
cambio_raw <- get_series(c(cambio = 3698),
                         start_date = inicio,
                         end_date = fim)

cambio_tri <- cambio_raw |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(
    valor = mean(cambio, na.rm = TRUE)) |>
  rename(cambio = valor) |>
  arrange(trimestre) |>
  mutate(
    log_cambio   = log(cambio),
    var_cambio   = (log_cambio - lag(log_cambio)) * 100   # var. log * 100 ≈ %
  ) |>
  left_join(metas_tri |> select(trimestre, dcambio_ppc), by = "trimestre") |>
  mutate(dcambio_desvio = var_cambio - dcambio_ppc)       # Δê_t (a Phillips usa t−1)

ggplot(cambio_tri, aes(x = trimestre, y = cambio)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  labs(title = "Câmbio R$/US$ — média trimestral",
       y = "R$/US$", x = NULL) +
  theme_minimal()

ggplot(cambio_tri |> filter(!is.na(dcambio_desvio)), aes(x = trimestre)) +
  geom_line(aes(y = var_cambio, color = "Variação trimestral"), linewidth = 0.6) +
  geom_line(aes(y = dcambio_ppc, color = "Variação de longo prazo (PPC)"), linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("Variação trimestral" = "darkorange",
                                "Variação de longo prazo (PPC)" = "darkblue")) +
  labs(title = "Câmbio: variação trimestral e referência da PPC",
       subtitle = "Δê = variação trimestral − (meta − 2)/4",
       y = "% no trimestre", x = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

#ICBR
icbr_raw <- get_series(
  c(icbr_geral = 27574, icbr_agro = 27575,
    icbr_metal = 27576, icbr_energia = 27577),
  start_date = inicio,
  end_date = fim
)

icbr_tri <- icbr_raw |>
  reduce(full_join, by = "date") |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(across(starts_with("icbr_"), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(trimestre) |>
  left_join(metas_tri |> select(trimestre, meta_trimestral), by = "trimestre") |>
  mutate(
    var_icbr_geral   = (log(icbr_geral)   - log(lag(icbr_geral)))   * 100,
    var_icbr_agro    = (log(icbr_agro)    - log(lag(icbr_agro)))    * 100,
    var_icbr_metal   = (log(icbr_metal)   - log(lag(icbr_metal)))   * 100,
    var_icbr_energia = (log(icbr_energia) - log(lag(icbr_energia))) * 100,
    icbr_desvio         = var_icbr_geral   - meta_trimestral,
    icbr_agro_desvio    = var_icbr_agro    - meta_trimestral,
    icbr_metal_desvio   = var_icbr_metal   - meta_trimestral,
    icbr_energia_desvio = var_icbr_energia - meta_trimestral
  ) |>
  select(-meta_trimestral)   # a meta entra no dados1 pela tabela de metas

icbr_long <- icbr_tri |>
  select(trimestre, starts_with("var_")) |>
  pivot_longer(-trimestre, names_to = "componente", values_to = "valor")

ggplot(icbr_long, aes(x = trimestre, y = valor, color = componente)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ componente, scales = "free_y", ncol = 2) +
  labs(title = "IC-Br — variação trimestral por componente",
       y = "% no trimestre", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(icbr_tri |> filter(!is.na(icbr_desvio)),
       aes(x = trimestre, y = icbr_desvio)) +
  geom_line(color = "darkblue", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Inflação importada (desvio da variação IC-Br vs meta)",
       subtitle = "Variável π̂* da Curva de Phillips",
       y = "p.p. no trimestre", x = NULL) +
  theme_minimal()

# Meta do ano vs meta no horizonte das expectativas
ggplot(metas_tri, aes(x = trimestre)) +
  geom_step(aes(y = meta, color = "Meta do ano-calendário"), linewidth = 0.8) +
  geom_line(aes(y = meta_12m, color = "Meta nos 12 meses à frente"), linewidth = 0.8) +
  scale_color_manual(values = c("Meta do ano-calendário" = "grey50",
                                "Meta nos 12 meses à frente" = "darkred")) +
  labs(title = "Meta de inflação",
       subtitle = "meta_12m entra na regra de Taylor e na equação de expectativas",
       y = "% a.a.", x = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# Consolidação
# ------------------------------------------------------------
dados1 <- metas_tri |>
  full_join(selic_tri,  by = "trimestre") |>
  full_join(ipca_tri,   by = "trimestre") |>
  full_join(cambio_tri |> select(-dcambio_ppc), by = "trimestre") |>
  full_join(icbr_tri,   by = "trimestre") |>
  arrange(trimestre) |>
  filter(trimestre >= as.Date("2002-01-01"), trimestre <= fim)

# Checagens: nenhuma coluna duplicada e nada faltando na janela 2002T3–2023T4
duplicadas <- grep("\\.(x|y)$", names(dados1), value = TRUE)
if (length(duplicadas) > 0) stop("Colunas duplicadas no dados1: ", paste(duplicadas, collapse = ", "))

janela <- dados1 |> filter(trimestre >= as.Date("2002-07-01"), trimestre <= as.Date("2023-10-01"))
cat(sprintf("\n=== NAs na janela 2002T3–2023T4 (%d trimestres) ===\n", nrow(janela)))
print(colSums(is.na(janela))[colSums(is.na(janela)) > 0])
cat("(vazio = nenhuma coluna com NA)\n")

cat("\nColunas do dados1:\n", paste(names(dados1), collapse = ", "), "\n")

write.csv(dados1, "dados1.csv",row.names = FALSE)