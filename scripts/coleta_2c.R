# ============================================================
# BLOCO 2c — Preços administrados (para a projeção E_t π da eq. 5)
# ============================================================
# A eq. (5) do boxe usa E_t π_{t,t+4}: a inflação do IPCA acumulada quatro
# trimestres à frente, consistente com o modelo. O modelo agregado projeta só
# os preços livres; para chegar ao IPCA:
#   IPCA_{t+j} = w_t · π^L_{t+j} + (1 − w_t) · π^adm_{t+j}
# Este script produz, por trimestre:
#   peso_livres          w_t: peso efetivo dos livres no IPCA, por regressão
#                        móvel (24 meses, só passado) da identidade
#                        IPCA − Adm = w · (Livres − Adm)
#   adm_exp_12m          expectativa Focus para os administrados nos 12 meses
#                        seguintes (média do trimestre), a partir das
#                        expectativas anuais e do acumulado no ano já realizado
#   adm_exp_12m_simples  interpolação f·E_ano + (1 − f)·E_ano+1 (robustez)
#   adm_fonte_focus      cobertura do trimestre pela pesquisa de administrados
#                        (0 = sem pesquisa: usa a expectativa Focus do IPCA)
#   adm_realizado_12m    administrados realizados nos 12 meses seguintes
#                        (só diagnóstico; não entra no modelo)
#   ipca_monit           variação trimestral dos administrados (monitorados)
# ============================================================

library(rbcb)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(purrr)

inicio <- as.Date("2000-01-01")
fim    <- as.Date("2025-12-31")
JANELA_PESO <- 24   # meses na regressão móvel do peso

# ------------------------------------------------------------
# Funções
# ------------------------------------------------------------

# Peso dos livres por regressão sem intercepto em janela móvel (só passado):
# (IPCA − Adm) = w · (Livres − Adm)
peso_movel <- function(ipca, livres, adm, janela = 24) {
  y <- ipca - adm; x <- livres - adm
  w <- rep(NA_real_, length(y))
  for (i in seq_along(y)) if (i >= janela) {
    k <- (i - janela + 1):i
    w[i] <- sum(x[k] * y[k]) / sum(x[k]^2)
  }
  w
}

# Expectativa para os 12 meses seguintes (m+1 a m+12) a partir das expectativas
# anuais. E_ano e E_ano1 em % no ano; ytd = acumulado realizado de janeiro a m−1 (%).
#   restante do ano corrente (meses m..12) = (1 + E_ano)/(1 + ytd) − 1
#   meses m+1..12: fração (12 − m)/(13 − m) do restante, em termos geométricos
#   meses 1..m do ano seguinte: fração m/12 de E_ano1
exp_12m_de_anuais <- function(E_ano, E_ano1, ytd, m) {
  r_resto <- (1 + E_ano / 100) / (1 + ytd / 100) - 1
  r_parte <- (1 + r_resto)^((12 - m) / (13 - m)) - 1
  ((1 + r_parte) * (1 + E_ano1 / 100)^(m / 12) - 1) * 100
}

# ------------------------------------------------------------
# 1. IPCA, livres e administrados (monitorados) — mensal
# ------------------------------------------------------------
mensal <- get_series(c(ipca = 433, livres = 11428, adm = 4449),
                     start_date = inicio, end_date = fim) |>
  reduce(full_join, by = "date") |>
  arrange(date) |>
  filter(!is.na(ipca), !is.na(livres), !is.na(adm))

mensal$peso_livres <- peso_movel(mensal$ipca, mensal$livres, mensal$adm, JANELA_PESO)

ajuste <- with(mensal |> filter(!is.na(peso_livres)),
               1 - sum((ipca - (peso_livres * livres + (1 - peso_livres) * adm))^2) /
                 sum((ipca - mean(ipca))^2))
cat(sprintf("Peso dos livres (janela de %d meses): média %.3f | mín %.3f | máx %.3f | R² da identidade %.4f\n",
            JANELA_PESO, mean(mensal$peso_livres, na.rm = TRUE),
            min(mensal$peso_livres, na.rm = TRUE), max(mensal$peso_livres, na.rm = TRUE), ajuste))
cat("Esperado: peso entre ~0,70 e ~0,80 e R² perto de 1.\n")

# Acumulado no ano dos administrados, de janeiro até o mês m (para m = 0..12)
ytd_adm <- mensal |>
  transmute(ano = year(date), m = month(date), adm) |>
  group_by(ano) |>
  arrange(m, .by_group = TRUE) |>
  mutate(ytd = (cumprod(1 + adm / 100) - 1) * 100) |>
  ungroup() |>
  select(ano, m, ytd)
ytd_adm <- bind_rows(ytd_adm, distinct(ytd_adm, ano) |> mutate(m = 0L, ytd = 0))

# Realizado nos 12 meses seguintes (m+1 a m+12) — só diagnóstico
log_adm <- log(1 + mensal$adm / 100)
mensal$adm_real_12m <- sapply(seq_along(log_adm), function(i) {
  k <- (i + 1):(i + 12)
  if (max(k) > length(log_adm)) NA_real_ else (exp(sum(log_adm[k])) - 1) * 100
})

# ------------------------------------------------------------
# 2. Focus — IPCA Administrados (anuais) e IPCA 12 meses (para antes da cobertura)
# ------------------------------------------------------------
adm_raw <- get_annual_market_expectations("IPCA Administrados",
                                          start_date = as.Date("2001-01-01"), end_date = fim)
if (nrow(adm_raw) == 0 || !any(adm_raw$indic == "IPCA Administrados"))
  stop("Focus sem 'IPCA Administrados'. Indicadores retornados: ",
       paste(unique(adm_raw$indic), collapse = ", "))

focus_adm <- adm_raw |>
  filter(indic == "IPCA Administrados", base == 0) |>
  transmute(data = as_date(date), ano_ref = as.integer(reference_date), mediana = median) |>
  mutate(horiz = case_when(ano_ref == year(data)     ~ "E_ano",
                           ano_ref == year(data) + 1 ~ "E_ano1")) |>
  filter(!is.na(horiz)) |>
  select(data, horiz, mediana) |>
  pivot_wider(names_from = horiz, values_from = mediana, values_fn = function(x) x[length(x)]) |>
  filter(!is.na(E_ano), !is.na(E_ano1)) |>
  mutate(ano = year(data), m = month(data)) |>
  left_join(ytd_adm |> mutate(m = m + 1L), by = c("ano", "m")) |>   # acumulado até m − 1
  mutate(adm_12m         = exp_12m_de_anuais(E_ano, E_ano1, ytd, m),
         adm_12m_simples = ((12 - m) / 12) * E_ano + (m / 12) * E_ano1) |>
  arrange(data)

cat(sprintf("\nFocus IPCA Administrados: %s a %s | %d datas | datas sem acumulado no ano: %d\n",
            format(min(focus_adm$data)), format(max(focus_adm$data)), nrow(focus_adm),
            sum(is.na(focus_adm$ytd))))

ipca12 <- get_twelve_months_inflation_expectations("IPCA", start_date = as.Date("2001-01-01"),
                                                   end_date = fim) |>
  filter(smoothed == "S", base == 0) |>
  transmute(data = as_date(date), ipca_12m = median)

# ------------------------------------------------------------
# 3. Trimestralização (média das datas de pesquisa do trimestre)
# ------------------------------------------------------------
# Cada série é trimestralizada com as suas próprias datas; a expectativa do
# IPCA só substitui a de administrados nos trimestres sem nenhuma pesquisa de
# administrados (antes da cobertura). Assim as duas nunca se misturam dentro
# de um trimestre coberto.
adm_q <- focus_adm |>
  mutate(trimestre = floor_date(data, "quarter")) |>
  group_by(trimestre) |>
  summarise(adm_q = mean(adm_12m), adm_simples_q = mean(adm_12m_simples),
            n_adm = n(), .groups = "drop")

ipca12_q <- ipca12 |>
  mutate(trimestre = floor_date(data, "quarter")) |>
  group_by(trimestre) |>
  summarise(ipca_q = mean(ipca_12m), n_ipca = n(), .groups = "drop")

adm_tri <- full_join(ipca12_q, adm_q, by = "trimestre") |>
  arrange(trimestre) |>
  transmute(trimestre,
            adm_exp_12m         = coalesce(adm_q, ipca_q),
            adm_exp_12m_simples = coalesce(adm_simples_q, ipca_q),
            # fração do trimestre coberta por pesquisas de administrados (0 = usa o IPCA)
            adm_fonte_focus     = pmin(coalesce(n_adm / n_ipca, ifelse(is.na(n_adm), 0, 1)), 1)) |>
  filter(!is.na(adm_exp_12m))

mensal_tri <- mensal |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(peso_livres       = mean(peso_livres),
            ipca_monit        = (prod(1 + adm / 100) - 1) * 100,
            adm_realizado_12m = mean(adm_real_12m),
            n = n(), .groups = "drop") |>
  filter(n == 3) |>
  select(-n)

dados2c <- full_join(mensal_tri, adm_tri, by = "trimestre") |>
  filter(trimestre >= as.Date("2002-01-01"), trimestre <= fim) |>
  arrange(trimestre)

# ------------------------------------------------------------
# 4. Conferência
# ------------------------------------------------------------
cat("\nPrimeiro trimestre com pesquisa de administrados:",
    format(min(dados2c$trimestre[dados2c$adm_fonte_focus > 0], na.rm = TRUE)),
    "| primeiro com cobertura de pelo menos 90% das datas:",
    format(min(dados2c$trimestre[dados2c$adm_fonte_focus >= 0.9], na.rm = TRUE)), "\n")

comp <- dados2c |> filter(adm_fonte_focus > 0)
cat(sprintf("Com acumulado no ano x interpolação simples: cor = %.3f | dif. média = %+.2f p.p. | |dif| máx = %.2f p.p.\n",
            cor(comp$adm_exp_12m, comp$adm_exp_12m_simples),
            mean(comp$adm_exp_12m - comp$adm_exp_12m_simples),
            max(abs(comp$adm_exp_12m - comp$adm_exp_12m_simples))))

janela <- dados2c |> filter(trimestre >= as.Date("2002-07-01"), trimestre <= as.Date("2023-10-01"))
cat("\nNAs na janela 2002T3–2023T4 (peso_livres e adm_exp_12m devem ter 0):\n")
print(colSums(is.na(janela[, c("peso_livres", "adm_exp_12m", "adm_exp_12m_simples", "ipca_monit")])))

print(
  dados2c |>
    select(trimestre,
           `Expectativa Focus (12 meses)`  = adm_exp_12m,
           `Interpolação simples`          = adm_exp_12m_simples,
           `Realizado nos 12 meses seguintes` = adm_realizado_12m) |>
    pivot_longer(-trimestre, names_to = "serie", values_to = "valor") |>
    ggplot(aes(trimestre, valor, color = serie, linetype = serie)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_color_manual(values = c("darkblue", "darkorange", "grey50")) +
    scale_linetype_manual(values = c("solid", "dashed", "solid")) +
    labs(title = "Preços administrados — expectativa para os 12 meses seguintes",
         subtitle = "Antes da cobertura do Focus de administrados, usa a expectativa Focus do IPCA",
         y = "% em 12 meses", x = NULL, color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

print(
  ggplot(dados2c, aes(trimestre, peso_livres)) +
    geom_line(color = "darkgreen", linewidth = 0.8, na.rm = TRUE) +
    labs(title = "Peso efetivo dos preços livres no IPCA",
         subtitle = paste0("Regressão móvel de ", JANELA_PESO, " meses, só passado"),
         y = NULL, x = NULL) +
    theme_minimal()
)

write.csv(dados2c, "dados2c_administrados.csv", row.names = FALSE)
cat("Gravado: dados2c_administrados.csv\n")