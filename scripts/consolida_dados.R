# CONSOLIDAÇÃO — base pronta para o 1_modelo.R

library(dplyr)
library(readr)
library(lubridate)

AMOSTRA_INI <- as.Date("2003-04-01")   # 2003T2
AMOSTRA_FIM <- as.Date("2023-10-01")   # 2023T4

ler <- function(arq) read_csv(arq, show_col_types = FALSE) |> mutate(trimestre = as.Date(trimestre))

painel <- ler("dados1.csv") |>
  full_join(ler("dados2a.csv"), by = "trimestre") |>
  full_join(ler("dados2b.csv"), by = "trimestre") |>
  full_join(ler("dados2c.csv"), by = "trimestre") |>
  full_join(ler("dados3.csv"),  by = "trimestre") |>
  arrange(trimestre) |>
  filter(trimestre >= as.Date("2002-01-01"), trimestre <= as.Date("2025-12-31"))

soma_lags    <- function(x, ks) Reduce(`+`, lapply(ks, function(k) dplyr::lag(x, k)))
media_lags   <- function(x, ks) soma_lags(x, ks) / length(ks)
acum4        <- function(x) (exp(soma_lags(log1p(x / 100), 0:3)) - 1) * 100   # % em 4 trimestres
lead_ult     <- function(x, k) coalesce(dplyr::lead(x, k), last(x))           # meta futura (conhecida)
clima_desvio <- function(x) media_lags(x, 0:2) - media_lags(x, 3:5)
f_ano_corrente <- c(10, 7, 4, 1) / 12   # fração do horizonte de 12 meses no ano corrente

# ---- Transformações sobre o painel completo (2002T1–2025T4) ----
painel <- painel |>
  mutate(
    dselic   = selic - dplyr::lag(selic),
    clima_el = clima_desvio(d_elnino * clima_sq * 100),     # décimos de grau
    clima_la = clima_desvio(d_lanina * clima_sq * 100),
    meta_exp = meta_12m + (8.5 - 4.0) * case_when(
      year(trimestre) == 2003 ~ f_ano_corrente[quarter(trimestre)],
      year(trimestre) == 2002 ~ 1 - f_ano_corrente[quarter(trimestre)],
      TRUE                    ~ 0),
    meta_cal = ifelse(year(trimestre) == 2003, 8.5, meta),
    # taxas anualizadas
    ipca_livres_sa_an = ipca_livres_sa * 4,
    icbr_desvio_an    = icbr_desvio * 4,
    var_cambio_an     = var_cambio * 4,
    ipca_an_l0 = ipca_cheio_sa * 4,
    ipca_an_l1 = dplyr::lag(ipca_an_l0, 1),
    ipca_an_l2 = dplyr::lag(ipca_an_l0, 2),
    ipca_an_l3 = dplyr::lag(ipca_an_l0, 3),
    ipca_cheio_ma4_lag1 = media_lags(ipca_an_l0, 1:4),
    ipca_acum4       = acum4(ipca_cheio_sa),
    ipca_4t_seguinte = dplyr::lead(ipca_acum4, 4),
    # metas defasadas e futuras
    meta_exp_lag1 = dplyr::lag(meta_exp),
    meta_exp_f1 = lead_ult(meta_exp, 1), meta_exp_f2 = lead_ult(meta_exp, 2),
    meta_exp_f3 = lead_ult(meta_exp, 3), meta_exp_f4 = lead_ult(meta_exp, 4),
    meta_cal_ma4_lag1 = media_lags(meta_cal, 1:4),
    meta_cal_l0 = meta_cal,
    meta_cal_l1 = dplyr::lag(meta_cal, 1),
    meta_cal_l2 = dplyr::lag(meta_cal, 2),
    meta_cal_l3 = dplyr::lag(meta_cal, 3),
    meta_cal_f1 = lead_ult(meta_cal, 1), meta_cal_f2 = lead_ult(meta_cal, 2),
    meta_cal_f3 = lead_ult(meta_cal, 3)
  )

# ---- Janela do modelo e defasagens das exógenas ----
dados_modelo <- painel |>
  filter(trimestre >= AMOSTRA_INI, trimestre <= AMOSTRA_FIM) |>
  mutate(exp_selic_lag1   = dplyr::lag(exp_selic_12m),
         r_trend_lag1     = dplyr::lag(juro_real_hp),
         ipca_livres_lag1 = dplyr::lag(ipca_livres_sa_an),
         dcambio_lag1     = dplyr::lag(var_cambio_an),
         selic_lag1       = dplyr::lag(selic),
         selic_lag2       = dplyr::lag(selic, 2),
         meta_lag1        = dplyr::lag(meta))

# ---- Checagem: nada faltando na amostra efetiva, exceto a desocupação antes de 2012 ----
efetiva <- dados_modelo |> filter(trimestre >= as.Date("2003-10-01"))
faltam  <- setdiff(names(efetiva)[colSums(is.na(efetiva)) > 0], "femp")
if (length(faltam)) warning("NAs na amostra efetiva (2003T4–2023T4): ", paste(faltam, collapse = ", "))
cat(sprintf("dados_modelo.csv: %s a %s, %d trimestres, %d colunas\n",
            format(min(dados_modelo$trimestre)), format(max(dados_modelo$trimestre)),
            nrow(dados_modelo), ncol(dados_modelo)))

write_csv(dados_modelo, "dados_modelo.csv")