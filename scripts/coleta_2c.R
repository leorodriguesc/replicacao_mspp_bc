library(rbcb)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(ggplot2)

inicio      <- as.Date("2000-01-01")
fim         <- as.Date("2025-12-31")
JANELA_PESO <- 24   # meses na regressão móvel do peso

# Peso dos livres por regressão sem intercepto em janela móvel (só passado)
peso_movel <- function(ipca, livres, adm, janela = 24) {
  y <- ipca - adm; x <- livres - adm
  w <- rep(NA_real_, length(y))
  for (i in seq_along(y)) if (i >= janela) {
    k <- (i - janela + 1):i
    w[i] <- sum(x[k] * y[k]) / sum(x[k]^2)
  }
  w
}

# Expectativa para os 12 meses seguintes (m+1 a m+12) a partir das anuais.
# E_ano e E_ano1 em % no ano; ytd = acumulado realizado de janeiro a m−1 (%).
exp_12m_de_anuais <- function(E_ano, E_ano1, ytd, m) {
  r_resto <- (1 + E_ano / 100) / (1 + ytd / 100) - 1
  r_parte <- (1 + r_resto)^((12 - m) / (13 - m)) - 1
  ((1 + r_parte) * (1 + E_ano1 / 100)^(m / 12) - 1) * 100
}

#IPCA, livres e administrados (monitorados), mensal
mensal <- get_series(c(ipca = 433, livres = 11428, adm = 4449), start_date = inicio, end_date = fim) |>
  reduce(full_join, by = "date") |>
  arrange(date) |>
  filter(!is.na(ipca), !is.na(livres), !is.na(adm))

mensal$peso_livres <- peso_movel(mensal$ipca, mensal$livres, mensal$adm, JANELA_PESO)
cat(sprintf("Peso dos livres: média %.3f | mín %.3f | máx %.3f (esperado entre ~0,70 e ~0,80)\n",
            mean(mensal$peso_livres, na.rm = TRUE), min(mensal$peso_livres, na.rm = TRUE),
            max(mensal$peso_livres, na.rm = TRUE)))

# Acumulado no ano dos administrados, de janeiro até o mês m (m = 0..12)
ytd_adm <- mensal |>
  transmute(ano = year(date), m = month(date), adm) |>
  group_by(ano) |>
  arrange(m, .by_group = TRUE) |>
  mutate(ytd = (cumprod(1 + adm / 100) - 1) * 100) |>
  ungroup() |>
  select(ano, m, ytd)
ytd_adm <- bind_rows(ytd_adm, distinct(ytd_adm, ano) |> mutate(m = 0L, ytd = 0))

# Focus: administrados (anuais) e IPCA 12 meses (antes da cobertura)
focus_adm <- get_annual_market_expectations("IPCA Administrados", start_date = as.Date("2001-01-01"),
                                            end_date = fim) |>
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
  mutate(adm_12m = exp_12m_de_anuais(E_ano, E_ano1, ytd, m))

if (nrow(focus_adm) == 0) stop("Focus sem 'IPCA Administrados'.")

ipca12 <- get_twelve_months_inflation_expectations("IPCA", start_date = as.Date("2001-01-01"),
                                                   end_date = fim) |>
  filter(smoothed == "S", base == 0) |>
  transmute(data = as_date(date), ipca_12m = median)

# Trimestralização (média das datas de pesquisa do trimestre)
# A expectativa do IPCA só substitui a de administrados nos trimestres sem
# nenhuma pesquisa de administrados; as duas nunca se misturam num trimestre.
media_tri <- function(df, coluna, nome) {
  df |>
    mutate(trimestre = floor_date(data, "quarter")) |>
    group_by(trimestre) |>
    summarise("{nome}" := mean(.data[[coluna]]), .groups = "drop")
}

adm_tri <- full_join(media_tri(ipca12, "ipca_12m", "ipca_q"),
                     media_tri(focus_adm, "adm_12m", "adm_q"), by = "trimestre") |>
  arrange(trimestre) |>
  transmute(trimestre, adm_exp_12m = coalesce(adm_q, ipca_q)) |>
  filter(!is.na(adm_exp_12m))

peso_tri <- mensal |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(peso_livres = mean(peso_livres), n = n(), .groups = "drop") |>
  filter(n == 3) |>
  select(-n)

dados2c <- full_join(peso_tri, adm_tri, by = "trimestre") |>
  filter(trimestre >= as.Date("2002-01-01"), trimestre <= fim) |>
  arrange(trimestre)

print(ggplot(dados2c, aes(trimestre, adm_exp_12m)) +
        geom_line(color = "darkblue", linewidth = 0.8, na.rm = TRUE) +
        labs(title = "Preços administrados — expectativa Focus para os 12 meses seguintes",
             subtitle = "Antes da pesquisa de administrados, usa a expectativa Focus do IPCA",
             y = "% em 12 meses", x = NULL) +
        theme_minimal())
print(ggplot(dados2c, aes(trimestre, peso_livres)) +
        geom_line(color = "darkgreen", linewidth = 0.8, na.rm = TRUE) +
        labs(title = "Peso efetivo dos preços livres no IPCA",
             subtitle = paste0("Regressão móvel de ", JANELA_PESO, " meses, só passado"), y = NULL, x = NULL) +
        theme_minimal())

janela <- dados2c |> filter(trimestre >= as.Date("2002-07-01"), trimestre <= as.Date("2023-10-01"))
if (any(is.na(janela$peso_livres) | is.na(janela$adm_exp_12m)))
  warning("NAs em peso_livres ou adm_exp_12m na janela 2002T3–2023T4.")

write.csv(dados2c, "dados2c.csv", row.names = FALSE)
cat("Gravado: dados2c.csv\n")