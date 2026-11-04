library(rbcb)
library(sidrar)
library(seasonal)
library(mFilter)
library(dplyr)       # >= 1.1.0 (join_by com desigualdades)
library(lubridate)
library(ggplot2)
library(tidyr)
library(jsonlite)
library(readr)

inicio <- as.Date("2002-01-01")
fim    <- as.Date("2025-12-31")

#PIB
pib_raw <- get_sidra(
  api = "/t/1620/n1/all/v/583/p/all/c11255/90707/d/v583%202"
)

pib_tri <- pib_raw |>
  as_tibble() |>
  select(periodo = `Trimestre (Código)`, valor = Valor) |>
  mutate(
    ano       = as.integer(substr(periodo, 1, 4)),
    trim_num  = as.integer(substr(periodo, 6, 6)),
    trimestre = make_date(ano, (trim_num - 1) * 3 + 1, 1),
    pib_index = as.numeric(valor)
  ) |>
  select(trimestre, pib_index) |>
  arrange(trimestre) |>
  filter(trimestre <= fim)   # sem corte inicial: X-13 e HP usam a série desde 1996

ggplot(pib_tri, aes(x = trimestre, y = pib_index)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  labs(title = "PIB — Índice de volume trimestral (sem ajuste sazonal)",
       y = "Índice (1995=100)", x = NULL) +
  theme_minimal()

#dessazonalizacao pib
pib_ts <- ts(pib_tri$pib_index,
             start = c(year(min(pib_tri$trimestre)),
                       quarter(min(pib_tri$trimestre))),
             frequency = 4)

pib_x13 <- seas(pib_ts)
summary(pib_x13)
pib_sa <- final(pib_x13)

# Comparação visual: original vs dessazonalizada
df_pib <- data.frame(
  trimestre = pib_tri$trimestre,
  original  = as.numeric(pib_ts),
  dessaz    = as.numeric(pib_sa)
) |>
  pivot_longer(c(original, dessaz), names_to = "tipo", values_to = "valor")

ggplot(df_pib, aes(x = trimestre, y = valor, color = tipo)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("original" = "gray60", "dessaz" = "darkblue")) +
  labs(title = "PIB: original (cinza) vs com ajuste sazonal (azul)",
       y = "Índice", x = NULL, color = NULL) +
  theme_minimal()

#filtro hp
pib_hp <- hpfilter(log(pib_sa), freq = 1600, type = "lambda")

df_hp_pib <- data.frame(
  trimestre = pib_tri$trimestre,
  log_pib   = as.numeric(log(pib_sa)),
  tendencia = as.numeric(pib_hp$trend),
  ciclo     = as.numeric(pib_hp$cycle) * 100   # em pontos percentuais do PIB tendencial
)

ggplot(df_hp_pib, aes(x = trimestre)) +
  geom_line(aes(y = log_pib, color = "log(PIB) dessaz"), linewidth = 0.8) +
  geom_line(aes(y = tendencia, color = "Tendência HP"), linewidth = 0.8) +
  scale_color_manual(values = c("log(PIB) dessaz" = "darkblue",
                                "Tendência HP" = "darkred")) +
  labs(title = "PIB dessazonalizado (azul) e tendência HP (vermelho) — em log",
       y = "log(PIB)", x = NULL, color = NULL) +
  theme_minimal()

ggplot(df_hp_pib, aes(x = trimestre, y = ciclo)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Componente cíclico do PIB (filtro HP, λ=1600)",
       subtitle = "Esta é a variável f_pib do modelo do BCB",
       y = "% do PIB tendencial", x = NULL) +
  theme_minimal()

cat(sprintf("\nPIB: X-13 e HP de %s a %s (%d trimestres)\n",
            format(min(df_hp_pib$trimestre)), format(max(df_hp_pib$trimestre)), nrow(df_hp_pib)))

pib_final <- df_hp_pib |>
  filter(trimestre >= inicio) |>   # corte só depois do filtro
  select(trimestre, pib_index = log_pib, pib_tendencia = tendencia, fpib = ciclo)

#EXPECTATIVAS IPCA
focus_ipca_raw <- get_twelve_months_inflation_expectations(
  indic = "IPCA",
  start_date = inicio,
  end_date = fim
)

focus_ipca_tri <- focus_ipca_raw |>
  filter(smoothed == "S", base == 0) |>      # suavizada, agregado total
  select(date, valor = median) |>
  arrange(date) |>                            # garantir ordem cronológica
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(exp_ipca_12m = mean(valor), .groups = "drop")   # boxe: média ao longo do trimestre

ggplot(focus_ipca_tri, aes(x = trimestre, y = exp_ipca_12m)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 3, linetype = "dashed", alpha = 0.5) +
  labs(title = "Focus — Expectativa IPCA 12 meses à frente",
       subtitle = "Mediana suavizada, média do trimestre",
       y = "% a.a.", x = NULL) +
  theme_minimal()

# ============================================================
# EXPECTATIVA DA SELIC — MÉDIA AO LONGO DOS PRÓXIMOS 12 MESES
# ============================================================

HORIZONTE      <- 365    # dias corridos
MAX_PREENCHIDO <- 0.25   # fração máxima do horizonte que pode vir de "herança"
# (reunião sem mediana no Focus herda o degrau anterior);
# acima disso a data usa a aproximação pelas expectativas anuais

# ------------------------------------------------------------
# Funções
# ------------------------------------------------------------

# Série diária do SGS em blocos (a API recusa janelas acima de 10 anos)
baixar_sgs_diaria <- function(codigo, nome, inicio, fim, anos_bloco = 2) {
  cortes <- seq(inicio, fim, by = paste(anos_bloco, "years"))
  fins   <- c(cortes[-1] - 1, fim)
  Map(function(a, b) get_series(setNames(codigo, nome), start_date = a, end_date = b),
      cortes, fins) |>
    bind_rows() |>
    distinct(date, .keep_all = TRUE) |>
    arrange(date)
}

# Meta vigente em cada data (último valor conhecido até a data, inclusive)
meta_vigente <- function(datas, meta_diaria) {
  i <- findInterval(as.numeric(datas), as.numeric(meta_diaria$date))
  ifelse(i > 0, meta_diaria$meta[pmax(i, 1)], NA_real_)
}

# Defasagem (dias) entre cada mudança da meta e a reunião imediatamente anterior
defasagem_vigencia <- function(datas_mudanca, datas_reuniao) {
  r  <- sort(datas_reuniao)
  i  <- findInterval(as.numeric(datas_mudanca), as.numeric(r))
  ok <- i > 0
  as.numeric(datas_mudanca[ok] - r[i[ok]])
}

# (a) Média da trajetória em degraus a partir do Focus por reunião.
#  focus_reuniao: data, reuniao_codigo, mediana
#  calendario:    reuniao_codigo, vigencia (primeiro dia com a nova meta)
#  Retorna, por data de pesquisa: s0, media_reuniao, ponto_reuniao, frac_herdada
media_trajetoria_reunioes <- function(focus_reuniao, calendario, meta_diaria,
                                      horizonte = 365, max_preenchido = 0.25) {
  datas <- focus_reuniao |>
    distinct(data) |>
    mutate(lim = data + horizonte,
           s0  = meta_vigente(data, meta_diaria))
  
  # Reuniões que passam a valer dentro de (d, d + horizonte)
  degraus <- datas |>
    inner_join(calendario |> select(reuniao_codigo, vigencia),
               by = join_by(data < vigencia, lim > vigencia)) |>
    left_join(focus_reuniao, by = c("data", "reuniao_codigo")) |>
    arrange(data, vigencia) |>
    group_by(data) |>
    mutate(herdado = is.na(mediana), valor = mediana) |>
    fill(valor, .direction = "down") |>
    mutate(valor   = coalesce(valor, s0),
           fim_deg = lead(vigencia, default = first(lim)),
           dur     = as.numeric(fim_deg - vigencia)) |>
    summarise(soma          = sum(valor * dur),
              dias_herdados = sum(dur[herdado]),
              primeira      = first(vigencia),   # já ordenado por vigência
              ponto         = last(valor),
              .groups = "drop")
  
  datas |>
    left_join(degraus, by = "data") |>
    mutate(dias_s0       = as.numeric(coalesce(primeira, lim) - data),
           media_reuniao = (coalesce(soma, 0) + s0 * dias_s0) / horizonte,
           ponto_reuniao = coalesce(ponto, s0),
           frac_herdada  = coalesce(dias_herdados, 0) / horizonte,
           valida        = frac_herdada <= max_preenchido,
           media_reuniao = if_else(valida, media_reuniao, NA_real_),
           ponto_reuniao = if_else(valida, ponto_reuniao, NA_real_)) |>
    select(data, s0, media_reuniao, ponto_reuniao, frac_herdada)
}

# (b) Média da trajetória linear S0 -> S1 (31/dez do ano corrente) -> S2 (31/dez do seguinte)
#  f = fração do horizonte até 31/dez. Integrando a trajetória em [0, 1]:
#   média = f (S0 + S1)/2 + (1 - f) S1 + (1 - f)^2 (S2 - S1)/2
media_trapezio <- function(s0, s1, s2, f) {
  f * (s0 + s1) / 2 + (1 - f) * s1 + (1 - f)^2 * (s2 - s1) / 2
}

#  anuais: data, ano_ref, mediana (Selic de fim de ano esperada para ano_ref)
media_trajetoria_anual <- function(anuais, meta_diaria, horizonte = 365) {
  anuais |>
    mutate(horiz = case_when(ano_ref == year(data)     ~ "s1",
                             ano_ref == year(data) + 1 ~ "s2")) |>
    filter(!is.na(horiz)) |>
    select(data, horiz, mediana) |>
    pivot_wider(names_from = horiz, values_from = mediana,
                values_fn = function(x) x[length(x)]) |>
    filter(!is.na(s1), !is.na(s2)) |>
    mutate(s0          = meta_vigente(data, meta_diaria),
           f           = pmin(pmax(as.numeric(make_date(year(data), 12, 31) - data) / horizonte, 0), 1),
           media_anual = media_trapezio(s0, s1, s2, f),
           ponto_anual = f * s1 + (1 - f) * s2) |>
    select(data, media_anual, ponto_anual) |>
    arrange(data)
}

# Emenda diária: por reunião onde é válida, anual no restante
emendar_fontes <- function(reun, anual) {
  full_join(reun |> select(data, media_reuniao, ponto_reuniao),
            anual, by = "data") |>
    arrange(data) |>
    mutate(media = coalesce(media_reuniao, media_anual),
           ponto = if_else(!is.na(media_reuniao), ponto_reuniao, ponto_anual),
           fonte = case_when(!is.na(media_reuniao) ~ "reuniao",
                             !is.na(media_anual)   ~ "anual"))
}

# ------------------------------------------------------------
# 1. Meta Selic diária (SGS 432)
# ------------------------------------------------------------
meta_diaria <- baixar_sgs_diaria(432, "meta", as.Date("2001-01-01"), fim)

# ------------------------------------------------------------
# 2. Calendário do Copom — passadas + futuras
# ------------------------------------------------------------
copom_passadas <- "https://www.bcb.gov.br/api/servico/sitebcb/copom/atas?quantidade=9999" |>
  fromJSON() |> purrr::chuck("conteudo") |> as_tibble()

copom_futuras <- paste0(
  "https://www.bcb.gov.br/api/servico/sitebcb/agendas/itemmaisvelho?lista=",
  "Reuni%C3%B5es%20do%20Copom"
) |> fromJSON() |> purrr::chuck("conteudo") |> as_tibble()

copom_passadas_clean <- copom_passadas |>
  mutate(data_reuniao = as_date(dataReferencia), ano = year(data_reuniao)) |>
  arrange(data_reuniao) |>
  select(data_reuniao, ano)

copom_futuras_clean <- copom_futuras |>
  filter(evento == "Reunião do Copom") |>
  mutate(data_reuniao = as_date(dataEvento), ano = year(data_reuniao)) |>
  group_by(ano_mes = floor_date(data_reuniao, "month")) |>
  filter(data_reuniao == max(data_reuniao)) |>
  ungroup() |>
  select(data_reuniao, ano) |>
  filter(data_reuniao > max(copom_passadas_clean$data_reuniao))

copom_calendario <- bind_rows(copom_passadas_clean, copom_futuras_clean) |>
  distinct(data_reuniao, .keep_all = TRUE) |>
  arrange(data_reuniao) |>
  group_by(ano) |>
  mutate(numero_reuniao = row_number()) |>
  ungroup() |>
  mutate(reuniao_codigo = paste0("R", numero_reuniao, "/", ano))

# ---- 2.1 Quando a nova meta passa a valer? ----
# Mede, desde 2006, a distância entre cada mudança da meta (SGS 432) e a
# reunião anterior. O esperado é 1 dia (decisão à noite, vigência no dia
# seguinte). Se `dataReferencia` fosse o 1º dia da reunião, apareceria 2.
datas_mudanca <- meta_diaria |>
  filter(meta != lag(meta), date >= as.Date("2006-01-01")) |>
  pull(date)

dif_vig <- defasagem_vigencia(datas_mudanca, copom_calendario$data_reuniao)
cat("\n=== Dias entre a reunião e a mudança da meta (desde 2006) ===\n")
print(table(cut(dif_vig, c(-Inf, 0, 1, 2, 3, 7, Inf),
                labels = c("<=0", "1", "2", "3", "4-7", ">7"))))

dif_ok <- dif_vig[dif_vig >= 1 & dif_vig <= 7]
DIAS_ATE_VIGENCIA <- if (length(dif_ok) > 0) {
  as.integer(names(which.max(table(dif_ok))))
} else 1L
cat("Vigência adotada: reunião +", DIAS_ATE_VIGENCIA, "dia(s)\n")
if (DIAS_ATE_VIGENCIA != 1)
  warning("Defasagem modal diferente de 1 dia: confira o campo dataReferencia das atas.")

copom_calendario <- copom_calendario |>
  mutate(vigencia = data_reuniao + DIAS_ATE_VIGENCIA)

# ------------------------------------------------------------
# 3. Focus — Selic por reunião
# ------------------------------------------------------------
url_selic_expec <- paste0(
  "https://olinda.bcb.gov.br/olinda/servico/Expectativas/versao/v1/odata/",
  "ExpectativasMercadoSelic?$filter=baseCalculo%20eq%200&$format=text/csv"
)

selic_expec_raw <- read_csv(
  url_selic_expec,
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE,
  progress = FALSE
)

focus_reuniao <- selic_expec_raw |>
  transmute(data           = as_date(Data),
            reuniao_codigo = Reuniao,
            mediana        = as.numeric(Mediana)) |>
  filter(!is.na(mediana)) |>
  distinct(data, reuniao_codigo, .keep_all = TRUE)

# ---- 3.1 Checagens de cobertura e de numeração das reuniões ----
cat("\n=== Focus por reunião ===\n")
cat("Cobertura:", format(min(focus_reuniao$data)), "a", format(max(focus_reuniao$data)),
    "|", n_distinct(focus_reuniao$data), "datas de pesquisa\n")

sem_calendario <- focus_reuniao |>
  anti_join(copom_calendario, by = "reuniao_codigo") |>
  distinct(reuniao_codigo)
cat("Códigos de reunião do Focus sem par no calendário:", nrow(sem_calendario), "\n")
if (nrow(sem_calendario) > 0) print(sem_calendario, n = 20)

# Pesquisa feita depois da vigência da própria reunião indica numeração trocada
apos_vigencia <- focus_reuniao |>
  inner_join(copom_calendario |> select(reuniao_codigo, vigencia), by = "reuniao_codigo") |>
  filter(data >= vigencia)
cat("Linhas com data de pesquisa >= vigência da reunião:", nrow(apos_vigencia),
    "(esperado: 0 ou muito poucas)\n")
if (nrow(apos_vigencia) > 0)
  print(count(apos_vigencia, ano = year(vigencia), name = "linhas"))

# ------------------------------------------------------------
# 4. Média pela trajetória em degraus (Focus por reunião)
# ------------------------------------------------------------
selic_reuniao <- media_trajetoria_reunioes(focus_reuniao, copom_calendario, meta_diaria,
                                           horizonte = HORIZONTE,
                                           max_preenchido = MAX_PREENCHIDO)

cat("\nFração do horizonte preenchida por herança (por ano):\n")
selic_reuniao |>
  group_by(ano = year(data)) |>
  summarise(media_herdada = round(mean(frac_herdada), 3),
            datas_descartadas = sum(frac_herdada > MAX_PREENCHIDO),
            .groups = "drop") |>
  print(n = Inf)

# ------------------------------------------------------------
# 5. Média pela trajetória linear (expectativas de fim de ano)
# ------------------------------------------------------------
anuais_raw <- get_annual_market_expectations("Selic", start_date = inicio, end_date = fim)

if (!any(anuais_raw$indic == "Selic"))
  stop("Indicador 'Selic' não encontrado. Nomes disponíveis: ",
       paste(unique(anuais_raw$indic), collapse = ", "))

anuais_selic <- anuais_raw |>
  filter(indic == "Selic", base == 0) |>
  transmute(data    = as_date(date),
            ano_ref = as.integer(reference_date),
            mediana = median) |>
  arrange(data)

selic_anual <- media_trajetoria_anual(anuais_selic, meta_diaria, horizonte = HORIZONTE)

# ------------------------------------------------------------
# 6. Emenda e trimestralização
# ------------------------------------------------------------
selic_diaria <- emendar_fontes(selic_reuniao, selic_anual) |>
  filter(data >= inicio, data <= fim)

# ---- 6.1 Validação: as duas aproximações onde ambas existem ----
sobrep_d <- selic_diaria |> filter(!is.na(media_reuniao), !is.na(media_anual))
sobrep_t <- sobrep_d |>
  mutate(trimestre = floor_date(data, "quarter")) |>
  group_by(trimestre) |>
  summarise(reuniao = mean(media_reuniao), anual = mean(media_anual), .groups = "drop")

cat(sprintf(paste0("\nSobreposição degraus x linear: %d dias / %d trimestres\n",
                   "  diária:     cor = %.3f | dif. média = %+.2f p.p. | |dif| p95 = %.2f p.p.\n",
                   "  trimestral: cor = %.3f | dif. média = %+.2f p.p. | |dif| máx = %.2f p.p.\n"),
            nrow(sobrep_d), nrow(sobrep_t),
            cor(sobrep_d$media_reuniao, sobrep_d$media_anual),
            mean(sobrep_d$media_reuniao - sobrep_d$media_anual),
            quantile(abs(sobrep_d$media_reuniao - sobrep_d$media_anual), 0.95),
            cor(sobrep_t$reuniao, sobrep_t$anual),
            mean(sobrep_t$reuniao - sobrep_t$anual),
            max(abs(sobrep_t$reuniao - sobrep_t$anual))))

# ---- 6.2 Série trimestral (substitui a anterior) ----
selic_tri_diag <- selic_diaria |>
  filter(!is.na(media)) |>
  mutate(trimestre = floor_date(data, "quarter")) |>
  group_by(trimestre) |>
  summarise(exp_selic_12m       = mean(media),             # média ao longo de 12 meses (principal)
            exp_selic_12m_ponto = mean(ponto, na.rm = TRUE), # ponto em t+12m (alternativa)
            frac_reuniao        = mean(fonte == "reuniao"),
            .groups = "drop")

selic_tri_diag <- selic_tri_diag |>
  mutate(exp_selic_12m_ponto = ifelse(is.nan(exp_selic_12m_ponto), NA_real_, exp_selic_12m_ponto))

# Principal (média) e alternativa (ponto) seguem para o painel
selic_focus_tri <- selic_tri_diag |> select(trimestre, exp_selic_12m, exp_selic_12m_ponto)

# ---- 6.3 Quanto a correção muda a série ----
cat("\nMédia (nova) − ponto em t+12m (antiga), trimestral:\n")
cat(sprintf("  |dif| média = %.2f p.p. | trimestres com |dif| > 0,5 p.p.: %d de %d\n",
            mean(abs(selic_tri_diag$exp_selic_12m - selic_tri_diag$exp_selic_12m_ponto)),
            sum(abs(selic_tri_diag$exp_selic_12m - selic_tri_diag$exp_selic_12m_ponto) > 0.5),
            nrow(selic_tri_diag)))
cat("  Maiores diferenças:\n")
selic_tri_diag |>
  mutate(dif = exp_selic_12m - exp_selic_12m_ponto) |>
  slice_max(abs(dif), n = 8) |>
  mutate(across(where(is.numeric), \(x) round(x, 2))) |>
  print()

print(
  selic_tri_diag |>
    select(trimestre,
           `Média em 12 meses (nova)` = exp_selic_12m,
           `Ponto em t+12m (antiga)`  = exp_selic_12m_ponto) |>
    pivot_longer(-trimestre, names_to = "medida", values_to = "valor") |>
    ggplot(aes(trimestre, valor, color = medida)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("darkblue", "grey55")) +
    labs(title = "Focus — Selic esperada para os próximos 12 meses",
         subtitle = "A diferença entre as duas aparece nos ciclos de aperto e afrouxamento",
         y = "% a.a.", x = NULL, color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

# 5. Juro real ex-ante (diferença simples, eq. 2.1 do boxe)
# ---- HP que ignora os NAs das pontas (o hpfilter não aceita NA) ----
hp_seguro <- function(x, lambda = 1600) {
  ok <- which(!is.na(x))
  if (length(ok) == 0) stop("Série inteiramente NA")
  if (any(diff(ok) != 1)) stop("NA no meio da série — o HP exige um trecho contínuo")
  tendencia <- rep(NA_real_, length(x))
  tendencia[ok] <- as.numeric(hpfilter(x[ok], freq = lambda, type = "lambda")$trend)
  tendencia
}

focus_consolidado <- focus_ipca_tri |>
  full_join(selic_focus_tri, by = "trimestre") |>
  arrange(trimestre) |>                                    # HP exige ordem cronológica
  mutate(
    juro_real_exante =  exp_selic_12m - exp_ipca_12m,
    juro_real_hp     = hp_seguro(juro_real_exante),
    # Alternativa: Selic esperada no ponto de 12 meses. A tendência HP é refeita
    # com ela, porque a Selic esperada entra na IS e também na tendência do juro
    # real que compõe a taxa neutra.
    juro_real_exante_ponto = exp_selic_12m_ponto - exp_ipca_12m,
    juro_real_hp_ponto     = hp_seguro(juro_real_exante_ponto)
  )

cat(sprintf(paste0("\nTendência HP do juro real — ponto menos média: dif. média = %+.2f p.p. | ",
                   "|dif| máx = %.2f p.p.\n"),
            mean(focus_consolidado$juro_real_hp_ponto - focus_consolidado$juro_real_hp, na.rm = TRUE),
            max(abs(focus_consolidado$juro_real_hp_ponto - focus_consolidado$juro_real_hp), na.rm = TRUE)))

ggplot(focus_consolidado |> filter(!is.na(juro_real_hp)), aes(x = trimestre)) +
  geom_line(aes(y = juro_real_exante, color = "Juro real ex-ante"),
            linewidth = 0.7) +
  geom_line(aes(y = juro_real_hp, color = "Tendência HP"),
            linewidth = 0.9) +
  scale_color_manual(values = c("Juro real ex-ante" = "darkblue",
                                "Tendência HP" = "darkred")) +
  labs(title = "Juro real ex-ante e tendência HP",
       subtitle = "r̄^trend será usado como input do modelo",
       y = "% a.a.", x = NULL, color = NULL) +
  theme_minimal()

dados2a <- pib_final |>
  full_join(focus_consolidado, by = "trimestre") |>
  arrange(trimestre)

write.csv(dados2a,"dados2a.csv", row.names = FALSE)