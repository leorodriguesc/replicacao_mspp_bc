# ============================================================
# BLOCO 3 — Variáveis externas e clima (versão 2)
# ============================================================
# MUDANÇA EM RELAÇÃO À VERSÃO ANTERIOR (apenas a seção ONI)
#  O arquivo da NOAA traz médias móveis de 3 meses (DJF, JFM, ...). A versão
#  anterior datava cada uma pelo mês central e tirava a média de três por
#  trimestre, o que no T1 cobria dezembro a abril (5 meses). Agora o ONI
#  trimestral é a estação que coincide com o trimestre-calendário:
#    T1 = JFM, T2 = AMJ, T3 = JAS, T4 = OND
#  e o quadrado e as dummies saem desse mesmo valor:
#    clima_sq = oni^2
#    d_elnino = 1(oni >= 0,5), d_lanina = 1(oni <= −0,5)    limiar NOAA (principal)
#    d_elnino_sinal, d_lanina_sinal = 1(oni > 0), 1(oni < 0) sinal (robustez)
#  Também avisa se o arquivo da NOAA parou de ser atualizado.
#
# Escala da variável climática: a decidir na fase do modelo (ver diagnóstico).
# EMBI e IGREA abaixo são provisórios: CDS e hiato mundial entram por fonte externa.
# ============================================================

library(rbcb)
library(fredr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(readr)
library(ipeadatar)
library(httr)

inicio <- as.Date("2002-01-01")
fim    <- as.Date("2025-12-31")

#FED FUNDS
fredr_set_key(Sys.getenv("FRED_API_KEY"))

fedfunds_raw <- fredr(
  series_id = "FEDFUNDS",
  observation_start = inicio,
  observation_end   = fim
)

fedfunds_tri <- fedfunds_raw |>
  rename(fedfunds = value) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(fedfunds = mean(fedfunds, na.rm = TRUE), .groups = "drop")

ggplot(fedfunds_tri, aes(x = trimestre, y = fedfunds)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Fed Funds Rate — Federal Reserve",
       subtitle = "Média trimestral, anualizada",
       y = "% a.a.", x = NULL) +
  theme_minimal()

#EMBI
embi_raw <- ipeadata(code = "JPM366_EMBI366") |>
  select(date, embi = value) |>
  arrange(date)

embi_tri <- embi_raw |>
  filter(date >= inicio, date <= fim) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(embi = mean(embi, na.rm = TRUE), .groups = "drop")

ggplot(embi_tri, aes(x = trimestre, y = embi)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 250, linetype = "dashed", alpha = 0.5,
             color = "darkgreen") +
  labs(title = "EMBI+ Brasil — risco-país (proxy do CDS 5y)",
       subtitle = "Média trimestral em pontos-base | Linha verde = nível ~estável (250 pb)",
       y = "pontos-base", x = NULL) +
  theme_minimal()

uip <- full_join(embi_tri, fedfunds_tri, by = "trimestre") |>
  mutate(juros_externos = fedfunds + embi/100,
         delta_juros_ext = juros_externos - lag(juros_externos))

#hiato do produto mundial - proxy
gecon_raw <- fredr(
  series_id = "IGREA",
  observation_start = inicio,
  observation_end   = fim
)

hiato_mundial_tri <- gecon_raw |>
  rename(hiato_mundial = value) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(hiato_mundial = mean(hiato_mundial, na.rm = TRUE), .groups = "drop")

ggplot(hiato_mundial_tri, aes(x = trimestre, y = hiato_mundial)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Hiato do produto mundial (proxy: Kilian's Global Real Economic Activity)",
       subtitle = "Média trimestral, índice em desvio da tendência",
       y = "Índice", x = NULL) +
  theme_minimal()

# ============================================================
# ONI — anomalia de temperatura do Pacífico (NOAA/CPC)
# ============================================================

# ONI trimestral a partir do arquivo da NOAA (colunas SEAS, YR, TOTAL, ANOM).
# Usa só as estações que coincidem com os trimestres-calendário.
oni_trimestral <- function(oni_raw, limiar = 0.5) {
  estacao_trimestre <- c(JFM = 1, AMJ = 2, JAS = 3, OND = 4)
  oni_raw |>
    filter(SEAS %in% names(estacao_trimestre)) |>
    mutate(trimestre = make_date(YR, (estacao_trimestre[SEAS] - 1) * 3 + 1, 1)) |>
    transmute(trimestre,
              oni            = ANOM,
              clima_sq       = oni^2,
              d_elnino       = as.integer(oni >= limiar),
              d_lanina       = as.integer(oni <= -limiar),
              d_elnino_sinal = as.integer(oni > 0),
              d_lanina_sinal = as.integer(oni < 0)) |>
    arrange(trimestre)
}

oni_raw <- read_table(
  "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
  show_col_types = FALSE,
  progress = FALSE
)

# ---- Atualização do arquivo ----
ultima_estacao <- oni_raw |> slice_tail(n = 1)
cat(sprintf("\nONI: última estação no arquivo = %s/%d\n", ultima_estacao$SEAS, ultima_estacao$YR))
if (ultima_estacao$YR < year(fim))
  warning("O arquivo do ONI termina antes de ", year(fim),
          ". A NOAA passou a divulgar o RONI; confira se oni.ascii.txt ainda é atualizado.")

oni_tri <- oni_trimestral(oni_raw) |>
  filter(trimestre >= inicio, trimestre <= fim)

cat("Trimestres por classificação (limiar ±0,5):\n")
print(table(classe = with(oni_tri, case_when(d_elnino == 1 ~ "El Niño",
                                             d_lanina == 1 ~ "La Niña",
                                             TRUE          ~ "Neutro"))))

ggplot(oni_tri, aes(x = trimestre, y = oni)) +
  annotate("rect", xmin = min(oni_tri$trimestre), xmax = max(oni_tri$trimestre),
           ymin = 0.5, ymax = Inf, alpha = 0.1, fill = "red") +
  annotate("rect", xmin = min(oni_tri$trimestre), xmax = max(oni_tri$trimestre),
           ymin = -Inf, ymax = -0.5, alpha = 0.1, fill = "blue") +
  geom_hline(yintercept = c(-0.5, 0, 0.5), linetype = "dashed", alpha = 0.5) +
  geom_line(color = "darkblue", linewidth = 0.7) +
  geom_point(color = "darkblue", size = 1) +
  labs(title = "ONI trimestral (estação do próprio trimestre: JFM, AMJ, JAS, OND)",
       subtitle = "Vermelho: El Niño (≥ 0,5) | Azul: La Niña (≤ −0,5)",
       y = "Anomalia (°C)", x = NULL) +
  theme_minimal()

#consolidando
dados3 <- uip |>
  full_join(hiato_mundial_tri, by = "trimestre") |>
  full_join(oni_tri, by = "trimestre") |>
  arrange(trimestre)

write.csv(dados3, "dados3_externas.csv", row.names = FALSE)