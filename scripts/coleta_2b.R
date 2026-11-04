# ============================================================
# BLOCO 2b — Atividade e fiscal (versão 4)
# Desocupação, Nuci, Caged e resultado primário
# ============================================================
# MUDANÇAS DAS VERSÕES 3 e 4 (desocupação e Nuci)
#  [6] Série principal: PNAD-C pura (NA antes de 2012T1; o filtro de Kalman
#      trata a ausência), desvio da média 2012T1–2019T4. A série encadeada
#      com a PME fica como robustez (sufixo _encadeada).
#      Motivo (v4): a desocupação cai de ~14% para ~7% entre 2003 e 2014, em
#      boa parte por mudança estrutural. Com desvio de uma média constante, a
#      encadeada marca 2004–07 como folga forte (femp ≈ −1 a −3), enquanto Nuci
#      e Caged indicam economia neutra ou aquecida. As duas médias quase não
#      diferem (≈ 0,25 p.p.), então a PNAD-C pura não perde nada no "zero".
#  [7] Encadeamento: inclinação pela regressão na sobreposição (2012T1–2015T4)
#      e nível ancorado na média de 2012, para não haver degrau na emenda.
#  [8] Média pré-Covid da Nuci (e da encadeada) em 2003T4–2019T4, início da
#      amostra de estimação do boxe.
#  [9] Coluna femp = −fdesoc (pró-cíclica). A eq. (8) do boxe pede um
#      observável que suba com o hiato (γ_emp ∈ [0;3]). fdesoc mantém o sinal
#      antigo porque o script do modelo hoje inverte o sinal; quando a inversão
#      sair de lá, use femp diretamente.
# [10] Checagens na janela de estimação com defasagens (2002T3–2023T4),
#      checagem de sinais entre os observáveis de atividade e comparação da
#      encadeada antes e depois de 2012.
#
# MUDANÇAS DA VERSÃO 2 (mantidas)
#  [1] Dessazonalização (X-13) de desocupação, Caged e — se necessário — Nuci.
#  [2] PNAD Contínua: a série 24369 é de trimestre MÓVEL datado pelo último
#      mês; o trimestre-calendário é o valor de mar/jun/set/dez.
#  [4] Caged e resultado primário: valor de fim de trimestre.
#  [5] Checagem de quebra do Novo Caged (jan/2020).
#
# Coleta até 2025T4 de propósito, para afastar a ponta do filtro HP.
# ============================================================

library(rbcb)
library(sidrar)
library(mFilter)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(forecast)
library(tsibble)
library(seasonal)

inicio        <- as.Date("2002-01-01")
fim           <- as.Date("2025-12-31")

# Janela das médias pré-Covid (define o "zero" de Nuci e desocupação)
media_ini     <- as.Date("2003-10-01")   # 2003T4: início da amostra do boxe
media_fim     <- as.Date("2019-10-01")   # 2019T4: fim do pré-Covid

# Janela de checagem: amostra do boxe (2003T4–2023T4) mais as defasagens
# do modelo (o termo climático vai até t−5 -> 2002T3)
checa_ini     <- as.Date("2002-07-01")
checa_fim     <- as.Date("2023-10-01")

# Nuci: "auto" decide pela amplitude sazonal; TRUE/FALSE força a escolha
NUCI_SA       <- "auto"
LIMIAR_SAZ    <- 0.5     # p.p.; abaixo disso a série já é tratada como dessazonalizada

# ------------------------------------------------------------
# Auxiliares
# ------------------------------------------------------------

# X-13 trimestral, com fallback para X-11
# Retorna a série ajustada e o padrão sazonal retirado (bruta − ajustada)
dessaz_tri <- function(x, inicio_serie, nome = "", ...) {
  x_ts <- ts(x, start = c(year(inicio_serie), quarter(inicio_serie)), frequency = 4)
  m <- tryCatch(
    seas(x_ts, ...),
    error = function(e) {
      message(nome, ": SEATS falhou (", conditionMessage(e), "); usando X-11")
      seas(x_ts, x11 = "", ...)
    }
  )
  sa  <- as.numeric(final(m))
  fat <- tapply(as.numeric(x_ts) - sa, as.numeric(cycle(x_ts)), mean)
  amp <- diff(range(fat))
  cat(sprintf("%-16s padrão sazonal (T1 T2 T3 T4): %s | amplitude %.2f\n",
              nome, paste(sprintf("%+.2f", fat), collapse = "  "), amp))
  list(sa = sa, fatores = fat, amplitude = amp)
}

# Média de x nos trimestres da janela [ini, fim]
media_janela <- function(x, trimestre, ini, fim) {
  mean(x[trimestre >= ini & trimestre <= fim], na.rm = TRUE)
}

# Encadeia a PME (regiões metropolitanas, até fev/2016) à PNAD-C (Brasil, desde 2012).
#  pnad: trimestre, desoc_pnadc | pme: trimestre, pme_sa (ambas dessazonalizadas)
#  - inclinação b: regressão desoc_pnadc ~ pme_sa na sobreposição
#    (corrige a diferença de amplitude entre metrópoles e Brasil);
#  - nível a: escolhido para que o resíduo médio no ano-âncora seja zero
#    (evita degrau na emenda; o intercepto da regressão pesa mais 2015, quando
#    as duas séries se descolam).
encadear_desocupacao <- function(pnad, pme, ano_ancora = 2012) {
  sobrep <- inner_join(pnad, pme, by = "trimestre")
  reg    <- lm(desoc_pnadc ~ pme_sa, data = sobrep)
  b      <- unname(coef(reg)[2])
  anc    <- sobrep |> filter(year(trimestre) == ano_ancora)
  a      <- mean(anc$desoc_pnadc - b * anc$pme_sa)
  emenda <- min(pnad$trimestre[!is.na(pnad$desoc_pnadc)])
  prim   <- sobrep |> filter(trimestre == emenda)
  
  serie <- full_join(pme, pnad, by = "trimestre") |>
    arrange(trimestre) |>
    mutate(pme_calibrada = a + b * pme_sa,
           desoc         = if_else(trimestre < emenda, pme_calibrada, desoc_pnadc))
  
  list(serie            = serie,
       a = a, b = b,
       r2               = summary(reg)$r.squared,
       n_sobreposicao   = nrow(sobrep),
       emenda           = emenda,
       degrau_regressao = unname(coef(reg)[1] + b * prim$pme_sa - prim$desoc_pnadc),
       degrau_ancorado  = a + b * prim$pme_sa - prim$desoc_pnadc)
}

cat("=== Padrões sazonais retirados ===\n")

# ============================================================
# 1. DESOCUPAÇÃO
# ============================================================
pnad_raw <- get_series(c(pnad = 24369), start_date = inicio, end_date = fim)
pme_raw  <- get_series(c(pme  = 10777), start_date = inicio, end_date = fim)

# ---- 1.1 PNAD Contínua: trimestre-calendário exato ----
pnad_tri <- pnad_raw |>
  filter(month(date) %in% c(3, 6, 9, 12)) |>          # mar = jan-fev-mar, etc.
  mutate(trimestre = floor_date(date, "quarter")) |>
  select(trimestre, desoc_bruta = pnad) |>
  filter(trimestre >= inicio, trimestre <= fim) |>
  arrange(trimestre)

ds_pnad <- dessaz_tri(pnad_tri$desoc_bruta, min(pnad_tri$trimestre), "PNAD-C",
                      transform.function = "none")
pnad_tri$desoc_pnadc <- ds_pnad$sa

# ---- 1.2 PME: taxa mensal (não móvel), média dos três meses; só trimestres completos ----
pme_tri <- pme_raw |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(pme = mean(pme), n = n(), .groups = "drop") |>
  filter(n == 3, trimestre >= inicio) |>
  select(-n) |>
  arrange(trimestre)

ds_pme <- dessaz_tri(pme_tri$pme, min(pme_tri$trimestre), "PME",
                     transform.function = "none")
pme_tri$pme_sa <- ds_pme$sa

# ---- 1.3 Encadeamento ----
enc <- encadear_desocupacao(pnad_tri |> select(trimestre, desoc_pnadc),
                            pme_tri  |> select(trimestre, pme_sa))

cat(sprintf(paste0("\nEncadeamento PME -> PNAD-C: %d trimestres de sobreposição | R² = %.3f",
                   " | b = %.3f | a (âncora 2012) = %.3f\n"),
            enc$n_sobreposicao, enc$r2, enc$b, enc$a))
cat(sprintf("Diferença na emenda (%s): intercepto da regressão %+.2f p.p. | âncora 2012 %+.2f p.p.\n",
            format(enc$emenda), enc$degrau_regressao, enc$degrau_ancorado))

# ---- 1.4 Componentes cíclicos ----
# Principal: PNAD-C pura (média 2012T1–2019T4). Robustez: encadeada (média 2003T4–2019T4).
media_desoc     <- media_janela(enc$serie$desoc_pnadc, enc$serie$trimestre,
                                enc$emenda, media_fim)
media_desoc_enc <- media_janela(enc$serie$desoc, enc$serie$trimestre, media_ini, media_fim)

desoc_tri <- enc$serie |>
  rename(desoc_encadeada = desoc) |>
  mutate(desoc            = desoc_pnadc,
         fdesoc           = desoc - media_desoc,                  # principal (contracíclico)
         femp             = -fdesoc,                              # pró-cíclico: entra direto na eq. (8)
         fdesoc_encadeada = desoc_encadeada - media_desoc_enc,    # robustez
         femp_encadeada   = -fdesoc_encadeada)

cat(sprintf("Média pré-Covid — PNAD-C pura (%s a %s): %.2f%% | encadeada (%s a %s): %.2f%%\n\n",
            format(enc$emenda), format(media_fim), media_desoc,
            format(media_ini), format(media_fim), media_desoc_enc))

# ============================================================
# 2. NUCI
# ============================================================
nuci_raw <- get_series(c(nuci = 28561), start_date = inicio, end_date = fim)

nuci_tri <- nuci_raw |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  group_by(trimestre) |>
  summarise(nuci_bruta = mean(nuci), n = n(), .groups = "drop") |>
  filter(n == 3, trimestre >= inicio, trimestre <= fim) |>
  select(-n) |>
  arrange(trimestre)

ds_nuci <- dessaz_tri(nuci_tri$nuci_bruta, min(nuci_tri$trimestre), "Nuci",
                      transform.function = "none")

usar_sa_nuci <- if (identical(NUCI_SA, "auto")) ds_nuci$amplitude > LIMIAR_SAZ else isTRUE(NUCI_SA)
nuci_tri$nuci <- if (usar_sa_nuci) ds_nuci$sa else nuci_tri$nuci_bruta
cat("Nuci:", if (usar_sa_nuci) "sazonalidade relevante — usando a série DESSAZONALIZADA"
    else "sazonalidade desprezível — fonte já dessazonalizada, mantida a original", "\n")

if (min(nuci_tri$trimestre) > media_ini)
  warning("A Nuci começa depois de ", format(media_ini),
          ": a média pré-Covid cobre uma janela menor que a da desocupação.")

media_nuci_pre <- media_janela(nuci_tri$nuci, nuci_tri$trimestre, media_ini, media_fim)
nuci_tri <- nuci_tri |> mutate(fnuci = nuci - media_nuci_pre)
cat(sprintf("Média pré-Covid da Nuci (%s a %s): %.2f%%\n",
            format(media_ini), format(media_fim), media_nuci_pre))

# ============================================================
# 3. CAGED (estoque de empregos formais)
# ============================================================
caged_raw <- get_series(c(caged_estoque = 28763), start_date = inicio, end_date = fim)

caged_tri <- caged_raw |>
  filter(month(date) %in% c(3, 6, 9, 12)) |>          # estoque no fim do trimestre
  mutate(trimestre = floor_date(date, "quarter")) |>
  select(trimestre, caged_estoque) |>
  filter(trimestre >= inicio, trimestre <= fim) |>
  arrange(trimestre)

# ---- 3.1 Checagem de quebra metodológica (Novo Caged a partir de jan/2020) ----
cat("\nChecagem de quebra — variação trimestral do estoque bruto (%):\n")
caged_tri |>
  mutate(var_pct = round(100 * (caged_estoque / lag(caged_estoque) - 1), 2)) |>
  filter(trimestre >= as.Date("2019-01-01"), trimestre <= as.Date("2020-12-31")) |>
  print()
cat("Um salto atípico em 2020T1 (antes da pandemia, que veio em abr–mai) indica quebra.\n\n")

# ---- 3.2 Dessazonalizar em escala 100·log, depois HP ----
ds_caged <- dessaz_tri(100 * log(caged_tri$caged_estoque), min(caged_tri$trimestre),
                       "Caged (100·log)", transform.function = "none")
caged_tri$log100_sa <- ds_caged$sa

caged_hp <- hpfilter(caged_tri$log100_sa, freq = 1600, type = "lambda")
caged_tri <- caged_tri |>
  mutate(log100_tendencia = as.numeric(caged_hp$trend),
         fcaged           = as.numeric(caged_hp$cycle))   # já em %

# ============================================================
# 4. RESULTADO PRIMÁRIO
# ============================================================
# Série acumulada em 12 meses (% do PIB): usar o valor de fim de trimestre.
# Sinal invertido: NFSP positiva = déficit -> rp positivo = superávit.
# PENDENTE (item 4 do diagnóstico): governo central, correção de outliers e
# ajuste cíclico.
rp_raw <- get_series(c(rp = 5783), start_date = inicio, end_date = fim)

rp_tri <- rp_raw |>
  mutate(rp = -rp) |>
  filter(month(date) %in% c(3, 6, 9, 12)) |>
  mutate(trimestre = floor_date(date, "quarter")) |>
  select(trimestre, rp) |>
  filter(trimestre >= inicio, trimestre <= fim) |>
  arrange(trimestre)

rp_hp <- hpfilter(rp_tri$rp, freq = 1600, type = "lambda")
rp_tri <- rp_tri |>
  mutate(rp_tendencia = as.numeric(rp_hp$trend),
         rp_hiato     = rp - rp_tendencia)

# ============================================================
# 5. GRÁFICOS DE CONFERÊNCIA
# ============================================================
print(
  pnad_tri |>
    pivot_longer(c(desoc_bruta, desoc_pnadc), names_to = "serie", values_to = "taxa") |>
    ggplot(aes(trimestre, taxa, color = serie)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c(desoc_bruta = "grey60", desoc_pnadc = "darkblue"),
                       labels = c(desoc_bruta = "Bruta", desoc_pnadc = "Dessazonalizada")) +
    labs(title = "Taxa de desocupação — PNAD Contínua",
         subtitle = "O dente de serra do primeiro trimestre deve sumir na série ajustada",
         y = "% da força de trabalho", x = NULL, color = NULL) +
    theme_minimal()
)

# A PME calibrada deve colar na PNAD-C na sobreposição, sem degrau em 2012T1
print(
  desoc_tri |>
    select(trimestre,
           `PNAD-C dessazonalizada`      = desoc_pnadc,
           `PME calibrada`               = pme_calibrada,
           `Encadeada (robustez)`        = desoc_encadeada) |>
    pivot_longer(-trimestre, names_to = "serie", values_to = "taxa") |>
    ggplot(aes(trimestre, taxa, color = serie, linetype = serie)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = enc$emenda, linetype = "dotted", alpha = 0.6) +
    scale_color_manual(values = c("darkblue", "darkorange", "grey45")) +
    scale_linetype_manual(values = c("solid", "dashed", "solid")) +
    labs(title = "Desocupação — encadeamento PME -> PNAD-C",
         subtitle = paste0("Linha pontilhada: emenda em ", format(enc$emenda)),
         y = "% da força de trabalho", x = NULL, color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

print(
  desoc_tri |>
    select(trimestre, fdesoc, fdesoc_encadeada) |>
    pivot_longer(-trimestre, names_to = "serie", values_to = "valor") |>
    ggplot(aes(trimestre, valor, color = serie)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c(fdesoc = "darkblue", fdesoc_encadeada = "darkred"),
                       labels = c(fdesoc = "PNAD-C pura, média 2012T1–2019T4 (principal)",
                                  fdesoc_encadeada = "Encadeada, média 2003T4–2019T4 (robustez)")) +
    labs(title = "Componente cíclico da desocupação (contracíclico)",
         subtitle = "No modelo entra femp = −fdesoc",
         y = "p.p.", x = NULL, color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

print(
  ggplot(nuci_tri, aes(trimestre, fnuci)) +
    geom_line(color = "darkgreen", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = "Componente cíclico da Nuci",
         subtitle = paste0("Média 2003T4–2019T4 = ", round(media_nuci_pre, 2), "%"),
         y = "p.p.", x = NULL) +
    theme_minimal()
)

print(
  ggplot(caged_tri, aes(trimestre, fcaged)) +
    geom_line(color = "darkorange", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = "Componente cíclico do Caged (dessazonalizado, HP λ=1600)",
         y = "% do estoque tendencial", x = NULL) +
    theme_minimal()
)

print(
  ggplot(rp_tri, aes(trimestre, rp_hiato)) +
    geom_line(color = "darkblue", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = "Desvio do resultado primário em relação à tendência",
         y = "p.p. do PIB", x = NULL) +
    theme_minimal()
)

# ============================================================
# 6. CONSOLIDAÇÃO
# ============================================================
# Colunas da desocupação:
#   desoc             PNAD-C dessazonalizada (NA antes de 2012)
#   fdesoc            desoc − média 2012T1–2019T4 (contracíclico; o modelo hoje inverte)
#   femp              −fdesoc (pró-cíclico; usar direto quando a inversão sair do modelo)
#   desoc_encadeada   PME calibrada até 2011T4, PNAD-C depois (robustez)
#   fdesoc_encadeada  desoc_encadeada − média 2003T4–2019T4
#   femp_encadeada    −fdesoc_encadeada
dados2b <- desoc_tri |>
  select(trimestre, desoc, fdesoc, femp,
         desoc_encadeada, fdesoc_encadeada, femp_encadeada) |>
  full_join(nuci_tri  |> select(trimestre, nuci_bruta, nuci, fnuci), by = "trimestre") |>
  full_join(caged_tri |> select(trimestre, caged_estoque, fcaged),   by = "trimestre") |>
  full_join(rp_tri    |> select(trimestre, rp, rp_hiato),            by = "trimestre") |>
  filter(trimestre >= inicio, trimestre <= fim) |>
  arrange(trimestre)

# ---- Checagens finais ----
janela <- dados2b |> filter(trimestre >= checa_ini, trimestre <= checa_fim)

cat(sprintf("\n=== NAs na janela de estimação com defasagens (%s a %s, %d trimestres) ===\n",
            format(checa_ini), format(checa_fim), nrow(janela)))
print(colSums(is.na(janela)))
cat("Esperado: desoc, fdesoc e femp com 38 NAs (2002T3–2011T4); todas as demais com 0.\n")

cat("\n=== Resumo dos componentes cíclicos na janela ===\n")
print(round(sapply(janela[, c("fdesoc", "fdesoc_encadeada", "fnuci", "fcaged", "rp_hiato")],
                   function(x) c(media = mean(x, na.rm = TRUE), dp = sd(x, na.rm = TRUE),
                                 min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))), 2))

# Os observáveis pró-cíclicos devem andar juntos (correlações positivas).
# Uma correlação negativa com femp indicaria sinal trocado.
cat("\n=== Correlações entre os observáveis pró-cíclicos (janela, pares disponíveis) ===\n")
print(round(cor(janela[, c("femp", "femp_encadeada", "fnuci", "fcaged")],
                use = "pairwise.complete.obs"), 2))

# Por que a encadeada não é a principal: antes de 2012 ela se descola de Nuci e Caged
cat("\n=== Encadeada por subperíodo: médias e correlações com Nuci e Caged ===\n")
janela |>
  mutate(periodo = if_else(trimestre < as.Date("2012-01-01"), "2002T3-2011T4", "2012T1-2023T4")) |>
  group_by(periodo) |>
  summarise(n                = n(),
            media_femp_enc   = mean(femp_encadeada),
            media_fnuci      = mean(fnuci),
            media_fcaged     = mean(fcaged),
            cor_femp_fnuci   = cor(femp_encadeada, fnuci),
            cor_femp_fcaged  = cor(femp_encadeada, fcaged),
            .groups = "drop") |>
  mutate(across(where(is.double), \(x) round(x, 2))) |>
  as.data.frame() |>
  print()

cat("\nÚltimo trimestre disponível:", format(max(dados2b$trimestre)), "\n")

write.csv(dados2b, "dados2b_atividade_fiscal.csv", row.names = FALSE)
cat("Gravado: dados2b_atividade_fiscal.csv\n")