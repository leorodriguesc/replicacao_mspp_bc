library(fredr)
library(dplyr)
library(lubridate)
library(readr)
library(readxl)
library(ggplot2)
library(tidyr)

inicio <- as.Date("2002-01-01")
fim    <- as.Date("2025-12-31")

#Juro externo e prêmio de risco
fredr_set_key(Sys.getenv("FRED_API_KEY"))
fedfunds_tri <- fredr(series_id = "FEDFUNDS", observation_start = inicio, observation_end = fim) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(fedfunds = mean(value, na.rm = TRUE), .groups = "drop")

cds_tri <- read_excel("dados_novos.xlsx", sheet = "cds") |>          # CDS diário -> média do trimestre
  rename(data = 1, cds = 2) |>
  mutate(trimestre = floor_date(as.Date(data), "quarter")) |>
  group_by(trimestre) |>
  summarise(cds = mean(cds, na.rm = TRUE), .groups = "drop")

externas <- fedfunds_tri |>
  left_join(cds_tri, by = "trimestre") |>
  arrange(trimestre) |>
  mutate(juros_externos  = fedfunds + cds / 100,
         delta_juros_ext = juros_externos - lag(juros_externos))

print(externas |>
        select(trimestre, `Fed funds (% a.a.)` = fedfunds, `CDS 5 anos (p.b.)` = cds,
               `Juro externo + CDS (% a.a.)` = juros_externos) |>
        pivot_longer(-trimestre, names_to = "serie", values_to = "valor") |>
        ggplot(aes(trimestre, valor)) +
        geom_line(color = "darkblue", linewidth = 0.8, na.rm = TRUE) +
        facet_wrap(~ serie, scales = "free_y", ncol = 1) +
        labs(title = "Juro externo e prêmio de risco — médias trimestrais", y = NULL, x = NULL) +
        theme_minimal())

#Hiato do produto mundial (Moraes, 2025), já trimestral
hiato_mundial_tri <- read_excel("dados_novos.xlsx", sheet = "hiato") |>
  rename(data = 1, hiato_mundial = 2) |>
  transmute(trimestre = floor_date(as.Date(data), "quarter"), hiato_mundial) |>
  filter(!is.na(trimestre))
if (anyDuplicated(hiato_mundial_tri$trimestre)) stop("Hiato mundial com mais de um valor por trimestre.")

print(ggplot(hiato_mundial_tri, aes(trimestre, hiato_mundial)) +
        geom_line(color = "darkblue", linewidth = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        labs(title = "Hiato do produto mundial — Moraes (2025)", y = "p.p.", x = NULL) +
        theme_minimal())

# ONI (NOAA/CPC): estação que coincide com o trimestre-calendário
#    T1 = JFM, T2 = AMJ, T3 = JAS, T4 = OND
oni_raw <- read_table("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
                      show_col_types = FALSE, progress = FALSE)
if (max(oni_raw$YR) < year(fim))
  warning("O arquivo do ONI termina antes de ", year(fim), ": confira se oni.ascii.txt ainda é atualizado.")

estacao_trimestre <- c(JFM = 1, AMJ = 2, JAS = 3, OND = 4)
oni_tri <- oni_raw |>
  filter(SEAS %in% names(estacao_trimestre)) |>
  transmute(trimestre = make_date(YR, (estacao_trimestre[SEAS] - 1) * 3 + 1, 1),
            oni       = ANOM,
            clima_sq  = ANOM^2,
            d_elnino  = as.integer(ANOM >= 0.5),
            d_lanina  = as.integer(ANOM <= -0.5)) |>
  filter(trimestre >= inicio, trimestre <= fim) |>
  arrange(trimestre)

print(ggplot(oni_tri, aes(trimestre, oni)) +
        annotate("rect", xmin = min(oni_tri$trimestre), xmax = max(oni_tri$trimestre),
                 ymin = 0.5, ymax = Inf, alpha = 0.1, fill = "red") +
        annotate("rect", xmin = min(oni_tri$trimestre), xmax = max(oni_tri$trimestre),
                 ymin = -Inf, ymax = -0.5, alpha = 0.1, fill = "blue") +
        geom_hline(yintercept = c(-0.5, 0, 0.5), linetype = "dashed", alpha = 0.5) +
        geom_line(color = "darkblue", linewidth = 0.7) +
        geom_point(color = "darkblue", size = 1) +
        labs(title = "ONI trimestral (estação do próprio trimestre: JFM, AMJ, JAS, OND)",
             subtitle = "Vermelho: El Niño (>= 0,5) | Azul: La Niña (<= −0,5)", y = "Anomalia (°C)", x = NULL) +
        theme_minimal())
oni_tri <- oni_tri |> select(-oni)

# Consolidação
dados3 <- externas |>
  full_join(hiato_mundial_tri, by = "trimestre") |>
  full_join(oni_tri, by = "trimestre") |>
  filter(trimestre >= inicio, trimestre <= fim) |>
  arrange(trimestre)

write.csv(dados3, "dados3.csv", row.names = FALSE)
cat("Gravado: dados3.csv\n")