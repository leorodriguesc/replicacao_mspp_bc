library(rbcb)
library(sidrar)
library(seasonal)
library(mFilter)
library(dplyr)       
library(lubridate)
library(tidyr)
library(jsonlite)
library(readr)
library(readxl)
library(ggplot2)

inicio <- as.Date("2002-01-01")
fim    <- as.Date("2025-12-31")

# PIB: ajuste sazonal e filtro HP sobre a série desde 1996
pib_tri <- get_sidra(api = "/t/1620/n1/all/v/583/p/all/c11255/90707/d/v583%202") |>
  as_tibble() |>
  transmute(trimestre = make_date(as.integer(substr(`Trimestre (Código)`, 1, 4)),
                                  (as.integer(substr(`Trimestre (Código)`, 6, 6)) - 1) * 3 + 1, 1),
            pib_index = as.numeric(Valor)) |>
  arrange(trimestre) |>
  filter(trimestre <= fim)   # sem corte inicial: X-13 e HP usam a série desde 1996

pib_ts <- ts(pib_tri$pib_index, start = c(year(min(pib_tri$trimestre)), quarter(min(pib_tri$trimestre))),
             frequency = 4)
pib_sa <- final(seas(pib_ts))
pib_hp <- hpfilter(log(pib_sa), freq = 1600, type = "lambda")

df_pib <- tibble(trimestre = pib_tri$trimestre,
                 original  = as.numeric(pib_ts),
                 dessaz    = as.numeric(pib_sa),
                 log_pib   = as.numeric(log(pib_sa)),
                 tendencia = as.numeric(pib_hp$trend),
                 fpib      = as.numeric(pib_hp$cycle) * 100)   # % do PIB tendencial

print(df_pib |>
        pivot_longer(c(original, dessaz), names_to = "tipo", values_to = "valor") |>
        ggplot(aes(trimestre, valor, color = tipo)) +
        geom_line(linewidth = 0.8) +
        scale_color_manual(values = c(original = "gray60", dessaz = "darkblue"),
                           labels = c(original = "Original", dessaz = "Com ajuste sazonal")) +
        labs(title = "PIB — índice de volume trimestral", y = "Índice (1995 = 100)", x = NULL, color = NULL) +
        theme_minimal())
print(ggplot(df_pib, aes(trimestre)) +
        geom_line(aes(y = log_pib, color = "log(PIB) dessazonalizado"), linewidth = 0.8) +
        geom_line(aes(y = tendencia, color = "Tendência HP"), linewidth = 0.8) +
        scale_color_manual(values = c("log(PIB) dessazonalizado" = "darkblue", "Tendência HP" = "darkred")) +
        labs(title = "PIB dessazonalizado e tendência HP — em log", y = "log(PIB)", x = NULL, color = NULL) +
        theme_minimal())
print(ggplot(df_pib, aes(trimestre, fpib)) +
        geom_line(color = "darkblue", linewidth = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        labs(title = "Componente cíclico do PIB (filtro HP, lambda = 1.600): fpib",
             y = "% do PIB tendencial", x = NULL) +
        theme_minimal())

pib_final <- df_pib |>
  select(trimestre, fpib) |>
  filter(trimestre >= inicio)                                       # corte só depois do filtro

# Focus: IPCA 12 meses à frente
focus_ipca_tri <- get_twelve_months_inflation_expectations(indic = "IPCA", start_date = inicio,
                                                           end_date = fim) |>
  filter(smoothed == "S", base == 0) |>      # suavizada, agregado total
  select(date, valor = median) |>
  arrange(date) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(exp_ipca_12m = mean(valor), .groups = "drop")

print(ggplot(focus_ipca_tri, aes(trimestre, exp_ipca_12m)) +
        geom_line(color = "darkred", linewidth = 0.8) +
        labs(title = "Focus — expectativa do IPCA 12 meses à frente",
             subtitle = "Mediana suavizada, média do trimestre", y = "% a.a.", x = NULL) +
        theme_minimal())

# Focus: Selic média esperada nos 12 meses seguintes
# Para cada data de pesquisa, a trajetória esperada fica na meta vigente até a
# próxima reunião do Copom e, depois de cada reunião, na mediana esperada para ela.
# Antes de nov/2004 (sem Focus por reunião), a trajetória é linear entre a meta
# vigente e as expectativas para o fim do ano corrente e do seguinte.
HORIZONTE      <- 365    # dias corridos
MAX_PREENCHIDO <- 0.25   # fração máxima do horizonte que pode vir de "herança"
# (reunião sem mediana herda o degrau anterior); acima disso vale a trajetória linear

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

# (a) Média da trajetória em degraus a partir do Focus por reunião
media_trajetoria_reunioes <- function(focus_reuniao, calendario, meta_diaria,
                                      horizonte = 365, max_preenchido = 0.25) {
  datas <- focus_reuniao |>
    distinct(data) |>
    mutate(lim = data + horizonte,
           s0  = meta_vigente(data, meta_diaria))
  
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
              primeira      = first(vigencia),
              .groups = "drop")
  
  datas |>
    left_join(degraus, by = "data") |>
    mutate(dias_s0       = as.numeric(coalesce(primeira, lim) - data),
           media_reuniao = (coalesce(soma, 0) + s0 * dias_s0) / horizonte,
           frac_herdada  = coalesce(dias_herdados, 0) / horizonte,
           media_reuniao = if_else(frac_herdada <= max_preenchido, media_reuniao, NA_real_)) |>
    select(data, media_reuniao)
}

# (b) Média da trajetória linear S0 -> S1 (31/dez do ano corrente) -> S2 (31/dez do seguinte)
#  f = fração do horizonte até 31/dez. Integrando a trajetória em [0, 1]:
#   média = f (S0 + S1)/2 + (1 - f) S1 + (1 - f)^2 (S2 - S1)/2
media_trapezio <- function(s0, s1, s2, f) {
  f * (s0 + s1) / 2 + (1 - f) * s1 + (1 - f)^2 * (s2 - s1) / 2
}

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
           media_anual = media_trapezio(s0, s1, s2, f)) |>
    select(data, media_anual) |>
    arrange(data)
}

# ---- 3.1 Meta Selic diária e calendário do Copom ----
meta_diaria <- baixar_sgs_diaria(432, "meta", as.Date("2001-01-01"), fim)

copom_passadas <- "https://www.bcb.gov.br/api/servico/sitebcb/copom/atas?quantidade=9999" |>
  fromJSON() |> purrr::chuck("conteudo") |> as_tibble() |>
  mutate(data_reuniao = as_date(dataReferencia), ano = year(data_reuniao)) |>
  arrange(data_reuniao) |>
  select(data_reuniao, ano)

copom_futuras <- paste0("https://www.bcb.gov.br/api/servico/sitebcb/agendas/itemmaisvelho?lista=",
                        "Reuni%C3%B5es%20do%20Copom") |>
  fromJSON() |> purrr::chuck("conteudo") |> as_tibble() |>
  filter(evento == "Reunião do Copom") |>
  mutate(data_reuniao = as_date(dataEvento), ano = year(data_reuniao)) |>
  group_by(ano_mes = floor_date(data_reuniao, "month")) |>
  filter(data_reuniao == max(data_reuniao)) |>
  ungroup() |>
  select(data_reuniao, ano) |>
  filter(data_reuniao > max(copom_passadas$data_reuniao))

copom_calendario <- bind_rows(copom_passadas, copom_futuras) |>
  distinct(data_reuniao, .keep_all = TRUE) |>
  arrange(data_reuniao) |>
  group_by(ano) |>
  mutate(reuniao_codigo = paste0("R", row_number(), "/", ano)) |>
  ungroup()

# Quando a nova meta passa a valer: defasagem modal entre a reunião e a mudança
# da meta no SGS (esperado: 1 dia, decisão à noite e vigência no dia seguinte)
datas_mudanca <- meta_diaria |>
  filter(meta != lag(meta), date >= as.Date("2006-01-01")) |>
  pull(date)
dif_ok <- defasagem_vigencia(datas_mudanca, copom_calendario$data_reuniao)
dif_ok <- dif_ok[dif_ok >= 1 & dif_ok <= 7]
DIAS_ATE_VIGENCIA <- if (length(dif_ok) > 0) as.integer(names(which.max(table(dif_ok)))) else 1L
if (DIAS_ATE_VIGENCIA != 1)
  warning("Defasagem modal diferente de 1 dia: confira o campo dataReferencia das atas.")
copom_calendario <- copom_calendario |> mutate(vigencia = data_reuniao + DIAS_ATE_VIGENCIA)

# ---- 3.2 Focus por reunião e expectativas de fim de ano ----
focus_reuniao <- read_csv(
  paste0("https://olinda.bcb.gov.br/olinda/servico/Expectativas/versao/v1/odata/",
         "ExpectativasMercadoSelic?$filter=baseCalculo%20eq%200&$format=text/csv"),
  locale = locale(decimal_mark = ","), show_col_types = FALSE, progress = FALSE) |>
  transmute(data = as_date(Data), reuniao_codigo = Reuniao, mediana = as.numeric(Mediana)) |>
  filter(!is.na(mediana)) |>
  distinct(data, reuniao_codigo, .keep_all = TRUE)

n_sem_par <- n_distinct(anti_join(focus_reuniao, copom_calendario, by = "reuniao_codigo")$reuniao_codigo)
if (n_sem_par > 0) warning(n_sem_par, " códigos de reunião do Focus sem par no calendário do Copom.")

anuais_selic <- get_annual_market_expectations("Selic", start_date = inicio, end_date = fim) |>
  filter(indic == "Selic", base == 0) |>
  transmute(data = as_date(date), ano_ref = as.integer(reference_date), mediana = median) |>
  arrange(data)

# ---- 3.3 Emenda diária e trimestralização ----
selic_diaria <- full_join(
  media_trajetoria_reunioes(focus_reuniao, copom_calendario, meta_diaria, HORIZONTE, MAX_PREENCHIDO),
  media_trajetoria_anual(anuais_selic, meta_diaria, HORIZONTE),
  by = "data") |>
  arrange(data) |>
  mutate(media = coalesce(media_reuniao, media_anual)) |>
  filter(data >= inicio, data <= fim)

# Validação: as duas aproximações onde ambas existem (nota de rodapé da seção 3.2)
sobrep <- selic_diaria |> filter(!is.na(media_reuniao), !is.na(media_anual))
cat(sprintf("Selic esperada, degraus x linear: %d dias de sobreposição | cor = %.3f | dif. média = %+.2f p.p.\n",
            nrow(sobrep), cor(sobrep$media_reuniao, sobrep$media_anual),
            mean(sobrep$media_reuniao - sobrep$media_anual)))

selic_focus_tri <- selic_diaria |>
  filter(!is.na(media)) |>
  mutate(trimestre = floor_date(data, "quarter")) |>
  group_by(trimestre) |>
  summarise(exp_selic_12m = mean(media), .groups = "drop")

meta_tri <- meta_diaria |>
  group_by(trimestre = floor_date(date, "quarter")) |>
  summarise(meta = mean(meta), .groups = "drop")
print(ggplot(left_join(selic_focus_tri, meta_tri, by = "trimestre"), aes(trimestre)) +
        geom_line(aes(y = meta, color = "Meta Selic vigente"), linewidth = 0.6) +
        geom_line(aes(y = exp_selic_12m, color = "Selic esperada (média de 12 meses)"), linewidth = 0.8) +
        scale_color_manual(values = c("Meta Selic vigente" = "grey55",
                                      "Selic esperada (média de 12 meses)" = "darkblue")) +
        labs(title = "Focus — Selic esperada para os 12 meses seguintes", y = "% a.a.", x = NULL, color = NULL) +
        theme_minimal() + theme(legend.position = "bottom"))

# Juro real ex-ante e tendência
# Tendência recalculada com a informação disponível na data do boxe: a série
# observada vai até FIM_INFO e é estendida pelo cenário do RI de jun/2024
# (Gráfico 2.2.2, juro real Focus quatro trimestres à frente, até 2026T4), como o
# BCB faz com a sua medida (Focus 1 ano à frente, filtro HP: 5,2 em 2023T2).
# O juro real é a taxa composta, a mesma do Gráfico 2.2.2. Na IS, o juro real
# continua a diferença simples (equação do boxe), porque pi^e é um estado e o
# modelo precisa ser linear nele.
ARQ_ANEXO <- "ri202406anp.xlsx"
FIM_INFO  <- as.Date("2024-01-01")

como_num  <- function(x) suppressWarnings(as.numeric(x))
serie_anexo <- function(nome, nomes) {
  a <- as.data.frame(suppressMessages(read_excel(ARQ_ANEXO, sheet = nome, col_names = FALSE)))
  a <- a[!is.na(como_num(a[[1]])) & como_num(a[[1]]) > 30000, ]
  d <- data.frame(trimestre = floor_date(as.Date(como_num(a[[1]]), origin = "1899-12-30"), "quarter"))
  for (j in seq_along(nomes)) d[[nomes[j]]] <- como_num(a[[j + 1]])
  d
}
filtro_hp <- function(x, lambda = 1600) {
  n <- length(x); D <- diff(diag(n), differences = 2)
  as.numeric(solve(diag(n) + lambda * crossprod(D), x))
}

focus_consolidado <- focus_ipca_tri |>
  full_join(selic_focus_tri, by = "trimestre") |>
  arrange(trimestre) |>
  mutate(juro_real_exante = exp_selic_12m - exp_ipca_12m)

cen      <- serie_anexo("Graf 2.2.2", c("ri", "ri_anterior"))
composta <- with(focus_consolidado, ((1 + exp_selic_12m / 100) / (1 + exp_ipca_12m / 100) - 1) * 100)
obs      <- which(!is.na(composta) & focus_consolidado$trimestre <= FIM_INFO)
tend     <- filtro_hp(c(composta[obs], cen$ri[cen$trimestre > FIM_INFO]))
focus_consolidado$juro_real_hp <- NA_real_
focus_consolidado$juro_real_hp[obs] <- tend[seq_along(obs)]

com <- match(cen$trimestre[cen$trimestre <= FIM_INFO], focus_consolidado$trimestre)
cat(sprintf("Juro real composto x Gráfico 2.2.2 até 2024T1: dif. média %+.2f p.p. (deve ficar perto de 0)\n",
            mean(composta[com] - cen$ri[cen$trimestre <= FIM_INFO])))
cat(sprintf("Tendência do juro real em 2023T2: %.2f (BCB: 5,2)\n",
            focus_consolidado$juro_real_hp[focus_consolidado$trimestre == as.Date("2023-04-01")]))

print(ggplot(filter(mutate(focus_consolidado, composta = composta), !is.na(juro_real_exante)), aes(trimestre)) +
        geom_line(aes(y = juro_real_exante, color = "Juro real ex-ante (diferença simples)"), linewidth = 0.5) +
        geom_line(aes(y = composta, color = "Juro real ex-ante (taxa composta)"), linewidth = 0.5) +
        geom_line(aes(y = juro_real_hp, color = "Tendência HP (composta)"), linewidth = 0.9) +
        geom_point(data = cen, aes(trimestre, ri, color = "Cenário do RI (Gráfico 2.2.2)"), size = 1) +
        scale_color_manual(values = c("Juro real ex-ante (diferença simples)" = "grey65",
                                      "Juro real ex-ante (taxa composta)" = "darkblue",
                                      "Tendência HP (composta)" = "darkred",
                                      "Cenário do RI (Gráfico 2.2.2)" = "black")) +
        labs(title = "Juro real ex-ante e tendência HP", y = "% a.a.", x = NULL, color = NULL) +
        theme_minimal() + theme(legend.position = "bottom") +
        guides(color = guide_legend(nrow = 2)))

# Consolidação
dados2a <- pib_final |>
  full_join(focus_consolidado, by = "trimestre") |>
  arrange(trimestre)

write.csv(dados2a, "dados2a.csv", row.names = FALSE)
cat("Gravado: dados2a.csv\n")