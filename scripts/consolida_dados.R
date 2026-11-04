library(dplyr)
library(readr)

ler <- function(arq) read_csv(arq, show_col_types = FALSE) |> mutate(trimestre = as.Date(trimestre))

onda1   <- ler("dados1.csv")
onda2a  <- ler("dados2a.csv")
onda2b  <- ler("dados2b_atividade_fiscal.csv")
onda2c  <- ler("dados2c_administrados.csv")
onda3   <- ler("dados3_externas.csv")

painel <- onda1 |>
  full_join(onda2a, by = "trimestre") |>
  full_join(onda2b, by = "trimestre") |>
  full_join(onda2c, by = "trimestre") |>
  full_join(onda3,  by = "trimestre") |>
  arrange(trimestre) |>
  filter(trimestre >= as.Date("2002-01-01"),
         trimestre <= as.Date("2025-12-31"))

write_csv(painel, "dados_consolidado.csv")