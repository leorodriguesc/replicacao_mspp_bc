library(KFAS)
library(dplyr)
library(lubridate)
library(readxl)
library(ggplot2)

# dados -------------------------------------------------------------------
AMOSTRA_INI    <- as.Date("2003-04-01")   
AMOSTRA_FIM    <- as.Date("2023-10-01")   
INICIO_EFETIVO <- 3                       

soma_lags  <- function(x, ks) Reduce(`+`, lapply(ks, function(k) dplyr::lag(x, k)))
media_lags <- function(x, ks) soma_lags(x, ks) / length(ks)
acum4      <- function(x) (exp(soma_lags(log1p(x / 100), 0:3)) - 1) * 100   # % em 4 trimestres
lead_ult   <- function(x, k) coalesce(dplyr::lead(x, k), last(x))           # meta futura (conhecida)
clima_desvio <- function(x) media_lags(x, 0:2) - media_lags(x, 3:5)       
f_ano_corrente <- c(10, 7, 4, 1) / 12   # fração do horizonte de 12 meses no ano corrente

ler_brutos <- function() {
  painel <- read.csv("dados_consolidado.csv") |>
    as_tibble() |> mutate(trimestre = as.Date(trimestre)) |> arrange(trimestre)
  cds_trim <- read_excel("dados_novos.xlsx", sheet = "cds") |>          # CDS diário -> média do trimestre
    rename(data = 1, cds = 2) |>
    mutate(trimestre = floor_date(as.Date(data), "quarter")) |>
    group_by(trimestre) |> summarise(cds = mean(cds, na.rm = TRUE), .groups = "drop")
  hiato_mundo <- read_excel("dados_novos.xlsx", sheet = "hiato") |>
    rename(data = 1, hiato_novo = 2) |>
    transmute(trimestre = floor_date(as.Date(data), "quarter"), hiato_novo)
  list(painel = painel, cds_trim = cds_trim, hiato_mundo = hiato_mundo)
}

preparar_dados <- function(e, brutos = BRUTOS) {
  painel <- brutos$painel
  if (e$selic_esperada == "ponto") {    
    painel$exp_selic_12m    <- painel$exp_selic_12m_ponto
    painel$juro_real_hp     <- painel$juro_real_hp_ponto
    painel$juro_real_exante <- painel$juro_real_exante_ponto
  }
  painel_anual <- painel |>
    select(-any_of(c("cds", "hiato_novo"))) |>
    left_join(brutos$cds_trim, by = "trimestre") |>
    left_join(brutos$hiato_mundo, by = "trimestre") |>
    arrange(trimestre) |>
    mutate(
      juros_externos  = fedfunds + cds / 100,
      delta_juros_ext = juros_externos - dplyr::lag(juros_externos),
      hiato_mundial   = hiato_novo,
      dselic          = selic - dplyr::lag(selic),
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
      # administrados na projeção E_t π
      ipca_12m_passado = if (e$adm_projecao == "focus") adm_exp_12m else acum4(ipca_cheio_sa),
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
  
  janela <- painel_anual |> filter(trimestre >= AMOSTRA_INI, trimestre <= AMOSTRA_FIM)
  
  #a desocupação entra com sinal invertido pois a equação é de emprego
  dados_obs <- janela |>
    transmute(trimestre, fpib, fnuci, fdesoc = -fdesoc, fcaged,
              pi_L = ipca_livres_sa_an, selic, dcambio = var_cambio_an, pi_focus = exp_ipca_12m)
  
  dados_exog <- janela |>
    mutate(exp_selic_lag1   = dplyr::lag(exp_selic_12m),
           r_trend_lag1     = dplyr::lag(juro_real_hp),
           ipca_livres_lag1 = dplyr::lag(ipca_livres_sa_an),
           dcambio_lag1     = dplyr::lag(var_cambio_an),
           selic_lag1       = dplyr::lag(selic),
           selic_lag2       = dplyr::lag(selic, 2),
           meta_lag1        = dplyr::lag(meta))
  
  efetivo <- INICIO_EFETIVO:nrow(dados_exog)
  
  # AR(1) das exógenas usadas na projeção de E_t π
  exog_ar  <- c("rp_hiato", "hiato_mundial", "icbr_desvio_an", "clima_el", "clima_la")
  ar_coefs <- setNames(lapply(exog_ar, function(v) {
    x <- dados_exog[[v]][efetivo]
    list(mu = mean(x, na.rm = TRUE), rho = cor(x[-1], x[-length(x)], use = "complete.obs"))
  }), exog_ar)
  
  list(painel_anual = painel_anual, dados_obs = dados_obs, dados_exog = dados_exog,
       efetivo = efetivo, ar_coefs = ar_coefs)
}

prev_ar <- function(v, x_t, j) ar_coefs[[v]]$mu + ar_coefs[[v]]$rho^j * (x_t - ar_coefs[[v]]$mu)

# secao 2 -----------------------------------------------------------------
#essa seção reúne num só lugar todas as escolhas que o boxe deixa em aberto 
#e cria as funções que trocam uma escolha por vez, para o 2_lacunas.R medir 
#quanto cada uma pesa nos resultados.

ESPEC_PRINCIPAL <- list(
  selic_esperada = "media",        
  adm_projecao   = "ipca12m",      
  sigma_eps_h    = 0.15,           # choque ordinário da IS
  sigma_pi_e     = 0.25,           # choque da equação de expectativas
  escala_ph      = 0.10,           # escala da variância da Phillips na pandemia
  crise_escala   = 1,     # [própria] escala dos termos de crise da IS: variâncias elevadas, como
  #   pede o boxe; os dados prefeririam 0,5, mas os parâmetros quase não mudam
  sigma_rr       = 0.15,  # [própria] choques das taxas neutras: o critério do boxe não é atingível,
  #   porque a variância da taxa neutra vem quase toda da tendência HP
  p1_neutro      = 1,     # [própria] variância inicial das taxas neutras: ancora a taxa neutra
  #   na tendência HP do juro real
  sigma_L        = 2.0264,         # erro da Phillips
  sigma_i        = 0.5124,         # erro da regra de Taylor
  sigma_e        = 26.1441,        # erro da paridade descoberta (câmbio anualizado)
  sigma_pi_e_obs = 0.01,  # [própria] erro de medida da Focus quase nulo: a eq. (5) fica com um
  #   único choque, sigma_pi_e
  pi_ss_externo  = 2.0    # [própria] inflação externa de equilíbrio: meta de 2% do Fed
)

# Especificação principal com algumas escolhas trocadas: espec(sigma_rr = 0.25)
espec <- function(...) {
  mud <- list(...)
  stopifnot(all(names(mud) %in% names(ESPEC_PRINCIPAL)))
  modifyList(ESPEC_PRINCIPAL, mud)
}

calibrados <- function(e) c(sigma_eps_h = e$sigma_eps_h, sigma_L = e$sigma_L, sigma_i = e$sigma_i,
                            sigma_e = e$sigma_e, sigma_pi_e = e$sigma_pi_e,
                            sigma_pi_e_obs = e$sigma_pi_e_obs,
                            sigma_rr_is = e$sigma_rr, sigma_rr_taylor = e$sigma_rr,
                            pi_ss_externo = e$pi_ss_externo)

ativar <- function(e = ESPEC_PRINCIPAL) {
  d <- preparar_dados(e)
  g <- globalenv()
  assign("ESPEC", e, envir = g)
  for (k in names(d)) assign(k, d[[k]], envir = g)
  assign("s_perfil", construir_perfis(d$dados_obs$trimestre, e), envir = g)
  assign("par_calibrados", calibrados(e), envir = g)
  invisible(e)
}

com_espec <- function(e, expr) {
  anterior <- ESPEC
  on.exit(ativar(anterior), add = TRUE)
  ativar(e)
  force(expr)
}

# ---- Rodar: lê e prepara os dados da especificação principal ----
BRUTOS <- ler_brutos()
ESPEC  <- ESPEC_PRINCIPAL
list2env(preparar_dados(ESPEC), envir = globalenv())   # painel_anual, dados_obs, dados_exog, efetivo, ar_coefs
cat("Amostra efetiva:", format(dados_obs$trimestre[efetivo[1]]), "a",
    format(max(dados_obs$trimestre)), "|", length(efetivo), "trimestres\n")
print(summary(dados_obs[efetivo, -1]))

# equacoes e espaco de estados --------------------------------------------
# Estados: h, r̂r^IS, r̂r^T, s, π^e, h_{t−1}, π^e_{t−1}, r̂r^IS_{t−1} e uma
# constante, que carrega os termos exógenos.
m <- 8; p <- 8; m_full <- 9
idx <- list(h = 1, r_is = 2, r_taylor = 3, s_h = 4, pi_e = 5,
            h_lag1 = 6, pi_e_lag1 = 7, r_is_lag1 = 8, const = 9)
obs <- list(fpib = 1, fnuci = 2, fdesoc = 3, fcaged = 4,
            pi_L = 5, selic = 6, dcambio = 7, pi_focus = 8)
nomes_obs <- names(obs)

# Matrizes base. Transição na convenção do KFAS: α_{t+1} = T α_t + R η_t.
construir_modelo <- function(par, ex) {
  pp <- as.list(par)
  coef_exp      <- 1 - pp$alpha_1L - pp$alpha_1I   # peso das expectativas na Phillips
  lambda_taylor <- 1 - pp$theta_1 - pp$theta_2
  
  # Observação: eqs. (6)–(9), Phillips, Taylor e Focus
  Z <- matrix(0, p, m)
  Z[obs$fpib,     idx$h]        <- 1
  Z[obs$fnuci,    idx$h]        <- pp$gamma_nuci
  Z[obs$fdesoc,   idx$h_lag1]   <- pp$gamma_emp
  Z[obs$fcaged,   idx$h_lag1]   <- pp$gamma_caged
  Z[obs$pi_L,     idx$h]        <- pp$alpha_4
  Z[obs$pi_L,     idx$pi_e]     <- coef_exp
  Z[obs$selic,    idx$r_taylor] <- lambda_taylor
  Z[obs$selic,    idx$pi_e]     <- lambda_taylor * pp$theta_3
  Z[obs$pi_focus, idx$pi_e]     <- 1
  
  # Erros de medida dos indicadores de atividade: desvio σ_h na série dividida por γ
  # (eqs. 6–9 do boxe), ou seja, γ σ_h na série original
  g_err <- c(1, pp$gamma_nuci, pp$gamma_emp, pp$gamma_caged)
  H <- diag(c((g_err * pp$sigma_h)^2, pp$sigma_L^2, pp$sigma_i^2, pp$sigma_e^2,
              pp$sigma_pi_e_obs^2))
  
  # Transição: IS, taxas neutras (passeio aleatório), eq. (5) e defasagens
  T_mat <- matrix(0, m, m)
  T_mat[idx$h, idx$h]               <- pp$beta1
  T_mat[idx$h, idx$pi_e]            <- pp$beta2 / 4
  T_mat[idx$h, idx$r_is]            <- pp$beta2 / 4
  T_mat[idx$h, idx$s_h]             <- pp$beta5
  T_mat[idx$s_h, idx$s_h]           <- pp$beta5
  T_mat[idx$r_is, idx$r_is]         <- 1
  T_mat[idx$r_taylor, idx$r_taylor] <- 1
  T_mat[idx$pi_e, idx$pi_e]         <- pp$phi_1
  T_mat[idx$h_lag1, idx$h]          <- 1
  T_mat[idx$pi_e_lag1, idx$pi_e]    <- 1
  T_mat[idx$r_is_lag1, idx$r_is]    <- 1
  
  # Choques: 1 IS (em s e no hiato), 2 neutra da IS, 3 neutra da Taylor,
  # 4 eq. (5), 5 termo de crise da IS (variância só nas janelas de crise)
  R_mat <- matrix(0, m, 5)
  R_mat[idx$s_h, 1] <- 1; R_mat[idx$h, 1] <- 1
  R_mat[idx$r_is, 2] <- 1
  R_mat[idx$r_taylor, 3] <- 1
  R_mat[idx$pi_e, 4] <- 1
  R_mat[idx$h, 5] <- 1
  Q_mat <- diag(c(pp$sigma_eps_h^2, pp$sigma_rr_is^2, pp$sigma_rr_taylor^2, pp$sigma_pi_e^2, 0))
  
  # Termos exógenos de cada equação (base anual)
  c_h <- -(pp$beta2 / 4) * ex$exp_selic_lag1 + (pp$beta2 / 4) * ex$r_trend_lag1 -
    pp$beta3 * ex$rp_hiato + pp$beta4 * ex$hiato_mundial
  c_pi <- pp$alpha_1L * ex$ipca_livres_lag1 + pp$alpha_1I * ex$ipca_cheio_ma4_lag1 +
    pp$alpha_2 * ex$icbr_desvio_an +
    pp$alpha_3 * (ex$dcambio_lag1 - (ex$meta_lag1 - pp$pi_ss_externo)) +   # câmbio em desvio da PPC
    pp$alpha_5 * ex$clima_el + pp$alpha_6 * ex$clima_la
  c_i <- pp$theta_1 * ex$selic_lag1 + pp$theta_2 * ex$selic_lag2 +
    lambda_taylor * ex$juro_real_hp + lambda_taylor * (1 - pp$theta_3) * ex$meta_exp
  c_e <- (ex$meta - pp$pi_ss_externo) - pp$delta * (ex$dselic - ex$delta_juros_ext)
  c_pi_e <- pp$phi_3 * ex$ipca_cheio_ma4_lag1 +         
    ((1 - pp$phi_2) * ex$meta_exp - pp$phi_1 * ex$meta_exp_lag1 - pp$phi_3 * ex$meta_cal_ma4_lag1)
  
  list(Z = Z, H = H, T = T_mat, R = R_mat, Q = Q_mat,
       c_h = c_h, c_pi = c_pi, c_i = c_i, c_e = c_e, c_pi_e = c_pi_e)
}

# Expectativa consistente com o modelo: E_t π_{t,t+4} = d_t + g_t'α_t.
# Itera o modelo quatro trimestres à frente com a informação de t. Cada grandeza é
# afim no estado (colunas: intercepto | estados | a própria E_t); como E aparece dos
# dois lados, E = A + b·E => E = A / (1 − b).
proj_afim <- function(par, ex, y, mod, H = 4) {
  pp <- as.list(par)
  n  <- nrow(ex); K <- m + 2
  cst  <- function(v) { M <- matrix(0, n, K); M[, 1] <- v; M }
  est  <- function(k) { M <- matrix(0, n, K); M[, 1 + k] <- 1; M }
  prev <- function(v, j) prev_ar(v, ex[[v]], j)
  
  coef_exp <- 1 - pp$alpha_1L - pp$alpha_1I
  w   <- ex$peso_livres
  adm <- ex$ipca_12m_passado
  
  h  <- est(idx$h);    s  <- est(idx$s_h)
  rr <- est(idx$r_is); pe <- est(idx$pi_e)
  E  <- matrix(0, n, K); E[, K] <- 1
  
  piL      <- coef_exp * pe + pp$alpha_4 * h + cst(mod$c_pi)
  ipca_ult <- list(cst(ex$ipca_an_l0 - w * y$pi_L) + w * piL,           # IPCA de t … t−3
                   cst(ex$ipca_an_l1), cst(ex$ipca_an_l2), cst(ex$ipca_an_l3))
  meta_f   <- cbind(ex$meta_exp_f1, ex$meta_exp_f2, ex$meta_exp_f3, ex$meta_exp_f4)
  meta_ant <- cbind(ex$meta_exp, meta_f[, 1:3])
  meta_cal <- cbind(ex$meta_cal_l3, ex$meta_cal_l2, ex$meta_cal_l1, ex$meta_cal_l0,
                    ex$meta_cal_f1, ex$meta_cal_f2, ex$meta_cal_f3)
  dcamb_t  <- ex$var_cambio_an - (ex$meta - pp$pi_ss_externo)
  
  soma <- matrix(0, n, K)
  for (j in seq_len(H)) {
    ma4 <- (ipca_ult[[1]] + ipca_ult[[2]] + ipca_ult[[3]] + ipca_ult[[4]]) / 4
    c_h <- -(pp$beta2 / 4) * ex$exp_selic_12m + (pp$beta2 / 4) * ex$juro_real_hp -
      pp$beta3 * prev("rp_hiato", j) + pp$beta4 * prev("hiato_mundial", j)
    h_novo <- pp$beta1 * h + (pp$beta2 / 4) * (pe + rr) + pp$beta5 * s + cst(c_h)
    s      <- pp$beta5 * s
    termo_meta <- (1 - pp$phi_2) * meta_f[, j] - pp$phi_1 * meta_ant[, j] -
      pp$phi_3 * rowMeans(meta_cal[, j:(j + 3)])
    pe_novo <- pp$phi_1 * pe + pp$phi_2 * E + pp$phi_3 * ma4 + cst(termo_meta)
    piL <- pp$alpha_1L * piL + pp$alpha_1I * ma4 + coef_exp * pe_novo + pp$alpha_4 * h_novo +
      cst(pp$alpha_2 * prev("icbr_desvio_an", j) + pp$alpha_3 * (if (j == 1) dcamb_t else 0) +
            pp$alpha_5 * prev("clima_el", j) + pp$alpha_6 * prev("clima_la", j))
    ipca_j   <- w * piL + cst((1 - w) * adm)
    soma     <- soma + ipca_j
    ipca_ult <- c(list(ipca_j), ipca_ult[1:3])
    h <- h_novo; pe <- pe_novo
  }
  
  A  <- soma / H
  b  <- A[, K]
  ok <- is.finite(rowSums(A))
  if (any(1 - b[ok] <= 0.05)) stop("Expectativa consistente sem solução estável.")
  A <- A[, -K, drop = FALSE] / (1 - b)
  G <- A[, -1, drop = FALSE]; colnames(G) <- names(idx)[1:m]
  list(d = A[, 1], G = G, b = b, ok = ok)
}

# A transição passa a ser A α_{t+1} = T α_t + R η, com A = I − φ2 e_πe g'.
# Por Sherman–Morrison, só a linha de π^e muda: linha(T*) = linha(T) + v'T, v = φ2 g / (1 − φ2 g_πe).
aplicar_expectativa <- function(T_full, R_full, proj, phi_2) {
  phi_2 <- unname(phi_2); ie <- idx$pi_e
  for (t in 1:dim(T_full)[3]) {
    tn <- min(t + 1, dim(T_full)[3])                  
    if (!proj$ok[tn]) next
    g   <- proj$G[tn, ]
    den <- 1 - phi_2 * g[ie]
    if (den <= 0.05) stop("Expectativa explosiva: 1 − φ2·g_πe = ", signif(den, 3))
    v <- phi_2 * g / den
    T_full[ie, idx$const, t] <- T_full[ie, idx$const, t] + phi_2 * proj$d[tn]
    T_full[ie, , t] <- T_full[ie, , t] + as.numeric(v %*% T_full[1:m, , t])
    R_full[ie, , t] <- R_full[ie, , t] + as.numeric(v %*% R_full[1:m, , t])
  }
  list(T = T_full, R = R_full)
}

# Variâncias que mudam no tempo (ideia de Lenza e Primiceri, 2022, citada no boxe)
construir_perfis <- function(datas, e = ESPEC) {
  crise <- c("2008-10-01" = 5, "2009-01-01" = 3, "2009-04-01" = 3, "2009-07-01" = 3,
             "2009-10-01" = 3, "2020-01-01" = 5, "2020-04-01" = 50, "2020-07-01" = 30,
             "2020-10-01" = 15)
  CR <- unname(crise[as.character(datas)]); CR[is.na(CR)] <- 0
  PH <- rep(1, length(datas))
  PH[datas == as.Date("2020-04-01")] <- 50
  PH[datas == as.Date("2020-07-01")] <- 30
  PH[datas == as.Date("2020-10-01")] <- 15
  decai <- datas >= as.Date("2021-01-01") & datas <= as.Date("2022-10-01")
  j     <- (year(datas) - 2020) * 4 + quarter(datas) - 4        # trimestres desde 2020T4
  PH[decai] <- 1 + (15 - 1) * 0.7^j[decai]
  list(CR = e$crise_escala * CR, PH = 1 + e$escala_ph * (PH - 1))
}

# Modelo do KFAS na amostra efetiva
construir_kfas_lp <- function(par, dados_exog, dados_obs, s_perfil) {
  mod <- construir_modelo(par, dados_exog)
  n <- nrow(dados_obs); k <- ncol(mod$R)
  Zt <- array(0, c(p, m_full, n)); Tt <- array(0, c(m_full, m_full, n))
  Rt <- array(0, c(m_full, k, n)); Qt <- array(mod$Q, c(k, k, n)); Ht <- array(mod$H, c(p, p, n))
  for (t in 1:n) {
    tn <- min(t + 1, n)                               # período de destino da transição
    Zt[, 1:m, t] <- mod$Z
    Zt[obs$pi_L,    idx$const, t] <- mod$c_pi[t]
    Zt[obs$selic,   idx$const, t] <- mod$c_i[t]
    Zt[obs$dcambio, idx$const, t] <- mod$c_e[t]
    Tt[1:m, 1:m, t] <- mod$T
    Tt[idx$const, idx$const, t] <- 1
    Tt[idx$h,    idx$const, t]  <- mod$c_h[tn]
    Tt[idx$pi_e, idx$const, t]  <- mod$c_pi_e[tn]
    Rt[1:m, , t] <- mod$R
    Ht[obs$pi_L, obs$pi_L, t] <- mod$H[obs$pi_L, obs$pi_L] * s_perfil$PH[t]^2
    Qt[k, k, t] <- s_perfil$CR[tn]^2
  }
  proj <- proj_afim(par, dados_exog, dados_obs, mod)
  aj   <- aplicar_expectativa(Tt, Rt, proj, par["phi_2"])
  
  # Condições iniciais [própria]: π^e começa na Focus do 1º trimestre (variância 4,
  # desvio de 2 p.p.); hiato e choque da IS começam em zero com variâncias grandes
  # (100 e 10), para que os dados dominem; taxas neutras com variância p1_neutro
  ia <- INICIO_EFETIVO:n
  a1 <- rep(0, m_full); a1[c(idx$pi_e, idx$pi_e_lag1)] <- dados_obs$pi_focus[ia[1]]; a1[idx$const] <- 1
  p1r <- ESPEC$p1_neutro
  P1 <- diag(c(100, p1r, p1r, 10, 4, 100, 4, p1r, 0))
  
  y_mat <- as.matrix(dados_obs[ia, nomes_obs])
  modelo <- SSModel(
    y_mat ~ -1 + SSMcustom(Z = Zt[, , ia], T = aj$T[, , ia], R = aj$R[, , ia, drop = FALSE],
                           Q = Qt[, , ia], a1 = a1, P1 = P1, P1inf = diag(0, m_full)),
    H = Ht[, , ia])
  list(modelo = modelo, idx_amostra = ia, proj = proj)
}

# ---- Rodar: perfis de variância e calibrados da especificação principal ----
s_perfil       <- construir_perfis(dados_obs$trimestre)
par_calibrados <- calibrados(ESPEC)
perfis <- data.frame(trimestre = dados_obs$trimestre, crise_IS = s_perfil$CR, phillips = round(s_perfil$PH, 2))
print(perfis[perfis$crise_IS > 0 | perfis$phillips > 1, ], row.names = FALSE)   # variância elevada

# priori e posteriori -----------------------------------------------------
# α4, α5, α6 e δ têm o suporte multiplicado por 4 (inflação e câmbio anualizados)
unif <- function(a, b) list(tipo = "unif", min = a, max = b)
prioris <- list(
  alpha_1L = unif(0, 1), alpha_1I = unif(0, 1), alpha_2 = unif(0, 1), alpha_3 = unif(0, 1),
  alpha_4  = unif(0, 4), alpha_5 = unif(0, 0.04), alpha_6 = unif(0, 0.04),
  beta1 = unif(0, 1), beta2 = unif(0, 1), beta3 = list(tipo = "beta", mean = 0.03, sd = 0.002),
  beta4 = unif(0, 1), beta5 = unif(0, 1),
  theta_1 = unif(0, 2), theta_2 = unif(-1, 1), theta_3 = unif(0, 4),
  phi_1 = unif(0, 1), phi_2 = unif(0, 1), phi_3 = unif(0, 1),
  delta = unif(0, 40), gamma_nuci = unif(0, 3), gamma_emp = unif(0, 3), gamma_caged = unif(0, 3),
  sigma_h = unif(0, 2))

log_prior <- function(par) {
  lp <- 0
  for (nome in names(prioris)) {
    pr <- prioris[[nome]]; x <- par[[nome]]
    if (pr$tipo == "unif") {
      if (x < pr$min || x > pr$max) return(-Inf)
      lp <- lp + dunif(x, pr$min, pr$max, log = TRUE)
    } else {
      k <- pr$mean * (1 - pr$mean) / pr$sd^2 - 1
      if (x <= 0 || x >= 1) return(-Inf)
      lp <- lp + dbeta(x, pr$mean * k, (1 - pr$mean) * k, log = TRUE)
    }
  }
  # Restrições implícitas: pesos não negativos; regra de Taylor estacionária
  if (par[["alpha_1L"]] + par[["alpha_1I"]] > 1)            return(-Inf)
  if (par[["phi_1"]] + par[["phi_2"]] + par[["phi_3"]] > 1) return(-Inf)
  if (par[["theta_1"]] + par[["theta_2"]] >= 1)             return(-Inf)
  if (par[["theta_2"]] - par[["theta_1"]] >= 1)             return(-Inf)
  lp
}

# Modas e intervalos de 90% do BCB (Tabela 1 do boxe, unidades do boxe)
par_canonico <- c(
  alpha_1L = 0.24, alpha_1I = 0.38, alpha_2 = 0.023, alpha_3 = 0.011,
  alpha_4 = 0.120, alpha_5 = 0.0012, alpha_6 = 0.0007,
  beta1 = 0.85, beta2 = 0.44, beta3 = 0.030, beta4 = 0.054, beta5 = 0.84,
  theta_1 = 1.48, theta_2 = -0.58, theta_3 = 2.03,
  phi_1 = 0.75, phi_2 = 0.11, phi_3 = 0.021,
  delta = 1.90, gamma_nuci = 1.87, gamma_emp = 1.10, gamma_caged = 0.69, sigma_h = 1.09)
nomes_theta <- names(par_canonico)
IC_BCB <- rbind(
  inf = c(0.02, 0.13, 0.006, 0, 0.072, 0.0004, 0, 0.70, 0.21, 0.027, 0, 0.59,
          1.41, -0.63, 1.47, 0.68, 0.06, 0, 0.77, 1.65, 0.96, 0.61, 0.98),
  sup = c(0.38, 0.69, 0.039, 0.025, 0.198, 0.0019, 0.0021, 0.95, 0.66, 0.032, 0.23, 0.99,
          1.54, -0.52, 2.64, 0.82, 0.13, 0.049, 3.22, 2.12, 1.25, 0.79, 1.21))
colnames(IC_BCB) <- nomes_theta

# Unidades do código (anualizadas) -> unidades do boxe: dividir por `escala`
escala <- setNames(ifelse(nomes_theta %in% c("alpha_4", "alpha_5", "alpha_6", "delta"), 4, 1),
                   nomes_theta)
par_bcb_base <- par_canonico * escala

combinar_par <- function(theta, pc = par_calibrados) c(theta[nomes_theta], pc)

log_posterior <- function(theta, pc = par_calibrados, sp = s_perfil) {
  names(theta) <- nomes_theta
  lp <- log_prior(theta)
  if (!is.finite(lp)) return(-Inf)
  ll <- tryCatch(as.numeric(logLik(construir_kfas_lp(combinar_par(theta, pc),
                                                     dados_exog, dados_obs, sp)$modelo)),
                 error = function(e) -Inf)
  if (!is.finite(ll)) return(-Inf)
  lp + ll
}

# ---- Rodar: log-posteriori nas modas do BCB ----
cat(sprintf("Log-posteriori nas modas do BCB: %.2f\n", log_posterior(par_bcb_base)))

# estimacao ---------------------------------------------------------------
# O modo é buscado por BFGS num espaço sem restrições, θ = g(u):
#   caixa [a, b]                 x = a + (b − a)·logística(u)
#   (α1L, α1I) e (φ1, φ2, φ3)    quebra de bastão: somas ≤ 1
#   (θ1, θ2)                     autocorrelações parciais: Taylor estacionária
eps <- 1e-6
limites_inf <- sapply(prioris, function(pr) if (pr$tipo == "unif") pr$min else 0.001)[nomes_theta] + eps
limites_sup <- sapply(prioris, function(pr) if (pr$tipo == "unif") pr$max else 0.999)[nomes_theta] - eps
limites_inf["sigma_h"] <- 0.001 + eps

logist  <- function(u) pmin(pmax(1 / (1 + exp(-u)), 1e-12), 1 - 1e-12)
logit   <- function(q) log(q / (1 - q))
apertar <- function(q, e = 1e-9) pmin(pmax(q, e), 1 - e)
BASTAO  <- list(c("alpha_1L", "alpha_1I"), c("phi_1", "phi_2", "phi_3"))
CAIXA   <- setdiff(nomes_theta, c(unlist(BASTAO), "theta_1", "theta_2"))

u_para_theta <- function(u) {
  u <- setNames(as.numeric(u), nomes_theta); th <- u
  a <- limites_inf[CAIXA]; b <- limites_sup[CAIXA]
  th[CAIXA] <- a + (b - a) * logist(u[CAIXA])
  for (g in BASTAO) { resto <- 1; for (k in g) { th[k] <- resto * logist(u[k]); resto <- resto - th[k] } }
  r1 <- logist(u["theta_1"]); r2 <- 2 * logist(u["theta_2"]) - 1
  th["theta_2"] <- r2; th["theta_1"] <- r1 * (1 - r2)
  th
}

theta_para_u <- function(th) {
  th <- setNames(as.numeric(th[nomes_theta]), nomes_theta); u <- th
  a <- limites_inf[CAIXA]; b <- limites_sup[CAIXA]
  u[CAIXA] <- logit(apertar((th[CAIXA] - a) / (b - a)))
  for (g in BASTAO) {
    resto <- 1
    for (k in g) { s_k <- apertar(th[k] / resto); u[k] <- logit(s_k); resto <- resto - s_k * resto }
  }
  r2 <- max(min(th["theta_2"], 1 - 1e-9), -1 + 1e-9)
  u["theta_2"] <- logit(apertar((r2 + 1) / 2)); u["theta_1"] <- logit(apertar(th["theta_1"] / (1 - r2)))
  u
}

# Menos a log-posteriori em u (o que o BFGS minimiza)
menos_lp_u <- function(u, pc = par_calibrados, sp = s_perfil) {
  v <- log_posterior(u_para_theta(u), pc, sp)
  if (is.finite(v)) -v else 1e10
}

otimizar <- function(theta_ini, pc = par_calibrados, sp = s_perfil, maxit = 2000) {
  o <- optim(theta_para_u(theta_ini), menos_lp_u, pc = pc, sp = sp, method = "BFGS",
             control = list(maxit = maxit, reltol = 1e-10))
  list(par = u_para_theta(o$par), value = o$value, convergence = o$convergence)
}

# ---- Modo ----
# A posteriori da IS é bimodal (persistência no hiato ou no choque): parte das
# modas do BCB e reinicia nas duas regiões de (β1, β2, β5).
estimar_modo <- function() {
  o1 <- otimizar(pmax(pmin(par_bcb_base, limites_sup), limites_inf))
  partidas <- list(modo_inicial = o1$par, regiao_bcb = o1$par, choque_persistente = o1$par)
  partidas$regiao_bcb[c("beta1", "beta2", "beta5")]         <- c(0.85, 0.44, 0.84)
  partidas$choque_persistente[c("beta1", "beta2", "beta5")] <- c(0.10, 0.30, 0.90)
  buscas <- lapply(partidas, otimizar)
  lp <- sapply(buscas, function(o) -o$value)
  print(t(sapply(buscas, function(o) round(c(o$par[c("beta1", "beta2", "beta5")],
                                             logpost = -o$value, codigo = o$convergence), 3))))
  melhor <- buscas[[which.max(lp)]]
  # Modo alternativo (outra região da IS): ponto de partida da segunda cadeia
  regiao <- function(th) th[["beta1"]] > th[["beta5"]]
  outras <- Filter(function(o) regiao(o$par) != regiao(melhor$par), buscas)
  alt <- if (length(outras)) outras[[which.max(sapply(outras, function(o) -o$value))]]$par else NULL
  saveRDS(list(theta_modo = melhor$par, theta_alt = alt, par_calibrados = par_calibrados,
               s_perfil = s_perfil, log_post_modo = -melhor$value), "modo_fase33.rds")
  cat("Gravado: modo_fase33.rds\n")
}

# Rodar: modo
if (!file.exists("modo_fase33.rds")) estimar_modo()
theta_modo <- readRDS("modo_fase33.rds")$theta_modo
cat(sprintf("Log-posteriori no modo: %.2f\n", log_posterior(theta_modo)))
print(round(data.frame(modo = theta_modo / escala, moda_bcb = par_canonico), 4))   # unidades do boxe

# ---- Posteriori: Metropolis adaptativo ----
dp_priori <- sapply(prioris[nomes_theta], function(pr)
  if (pr$tipo == "unif") (pr$max - pr$min) / sqrt(12) else pr$sd)

# Parâmetro na fronteira (|u| > 5): a partida é trazida para dentro até a
# log-posteriori cair 0,5; o deslocamento d mede a escala da posteriori nessa direção.
afastar_fronteira <- function(theta) {
  lp_u <- function(u) log_posterior(u_para_theta(u))
  u <- theta_para_u(theta); desloc <- list()
  for (k in which(abs(u) > 5)) {
    lp_ref <- lp_u(u)
    passo  <- function(s) replace(u, k, u[k] + s * (sign(u[k]) * 2 - u[k]))
    lo <- 1
    if (lp_u(passo(1)) < lp_ref - 0.5) {                   # bissecção no segmento
      lo <- 0; hi <- 1
      for (i in 1:20) { meio <- (lo + hi) / 2; if (lp_u(passo(meio)) >= lp_ref - 0.5) lo <- meio else hi <- meio }
    }
    desloc[[nomes_theta[k]]] <- u_para_theta(passo(lo)) - u_para_theta(u)
    u <- passo(lo)
  }
  list(theta = setNames(u_para_theta(u), nomes_theta), desloc = desloc)
}

# Proposta inicial: inversa da hessiana no modo (em u), levada a θ pelo jacobiano,
# mais (2d)(2d)' nas direções de fronteira. Direções ainda singulares recebem 10%
# do desvio a priori.
proposta_inicial <- function(theta_modo) {
  u0 <- theta_para_u(theta_modo)
  Hs <- optimHess(u0, menos_lp_u); Hs <- (Hs + t(Hs)) / 2
  ev <- eigen(Hs, symmetric = TRUE)
  if (max(ev$values) <= 0) stop("Hessiana sem curvatura positiva no modo.")
  if (any(ev$values <= 1e-8))
    Hs <- ev$vectors %*% diag(pmax(ev$values, max(ev$values) * 1e-6)) %*% t(ev$vectors)
  J <- sapply(seq_along(u0), function(k) {
    e <- replace(numeric(length(u0)), k, 1e-5)
    (u_para_theta(u0 + e) - u_para_theta(u0 - e)) / 2e-5
  })
  W  <- J %*% solve(Hs) %*% t(J)
  af <- afastar_fronteira(theta_modo)
  for (d in af$desloc) W <- W + tcrossprod(2 * d)
  
  C   <- ((W + t(W)) / 2) / outer(dp_priori, dp_priori)
  ev  <- eigen(C, symmetric = TRUE)
  lam <- ifelse(ev$values <= max(ev$values) * 1e-10, 0.1^2, ev$values)
  W   <- (ev$vectors %*% diag(lam) %*% t(ev$vectors)) * outer(dp_priori, dp_priori)
  list(Sigma = 2.38^2 / length(u0) * (W + t(W)) / 2,          # escala ótima do passeio aleatório
       inicio = af$theta)
}

# Cadeia A: passada curta com a proposta inicial e passada longa com a covariância
# empírica da primeira. Cadeia B: parte da outra região da IS. As duas são
# combinadas se o R-hat de todos os parâmetros for menor que 1,1.
estimar_posteriori <- function() {
  library(adaptMCMC); library(coda)
  md <- readRDS("modo_fase33.rds")
  lp <- function(th) log_posterior(th)
  prop <- proposta_inicial(md$theta_modo)
  
  set.seed(2024)
  c1 <- MCMC(lp, n = 100000, init = prop$inicio, scale = prop$Sigma, adapt = 25000, acc.rate = 0.234)
  Sigma2 <- 2.38^2 / length(nomes_theta) * cov(c1$samples[25001:100000, ])
  
  set.seed(2027)
  cA <- MCMC(lp, n = 250000, init = prop$inicio, scale = Sigma2, adapt = 50000, acc.rate = 0.234)
  ini_b <- if (!is.null(md$theta_alt)) md$theta_alt else
    replace(md$theta_modo, c("beta1", "beta2", "beta5"), c(0.10, 0.30, 0.90))
  set.seed(2031)
  cB <- MCMC(lp, n = 250000, init = afastar_fronteira(ini_b)$theta, scale = Sigma2, adapt = 50000, acc.rate = 0.234)
  
  A <- cA$samples[50001:250000, ]; B <- cB$samples[50001:250000, ]
  rhat <- gelman.diag(mcmc.list(mcmc(A), mcmc(B)), autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
  cat(sprintf("Aceitação: A %.3f | B %.3f | R-hat máximo %.3f\n",
              cA$acceptance.rate, cB$acceptance.rate, max(rhat)))
  cat(sprintf("Fração com β1 > β5: A %.2f | B %.2f\n",
              mean(A[, "beta1"] > A[, "beta5"]), mean(B[, "beta1"] > B[, "beta5"])))
  if (max(rhat) < 1.1) amostras <- rbind(A, B) else {
    warning("R-hat >= 1,1: cadeias não combinadas (fica só a A). Examinar antes de seguir.")
    amostras <- A
  }
  cat("ESS mínimo:", round(min(effectiveSize(mcmc(amostras)))), "\n")
  saveRDS(list(amostras = amostras), "posteriori_final.rds")
  cat("Gravado: posteriori_final.rds\n")
}

# Rodar: posteriori
if (!file.exists("posteriori_final.rds")) estimar_posteriori()
amostras <- readRDS("posteriori_final.rds")$amostras
colnames(amostras) <- nomes_theta
cat("Extrações da posteriori:", nrow(amostras), "\n")

# resultados --------------------------------------------------------------
# Identificação: razão entre os desvios-padrão da posteriori e da priori
# (abaixo de 0,25 = bem identificado)
amostras_boxe <- sweep(amostras, 2, escala, "/")          # unidades do boxe
tab <- data.frame(
  parametro = nomes_theta,
  modo      = theta_modo / escala,
  media     = colMeans(amostras_boxe),
  ic_inf    = apply(amostras_boxe, 2, quantile, .05),
  ic_sup    = apply(amostras_boxe, 2, quantile, .95),
  moda_bcb  = par_canonico,
  ic_bcb    = sprintf("[%s; %s]", IC_BCB["inf", ], IC_BCB["sup", ]),
  contracao = round(apply(amostras, 2, sd) / dp_priori, 2))
tab$compativel <- tab$moda_bcb >= tab$ic_inf & tab$moda_bcb <= tab$ic_sup
tab$bem_ident  <- tab$contracao < 0.25 & tab$parametro != "beta3"   # β3 tem priori informativa

tab_imp <- tab; num <- sapply(tab, is.numeric); tab_imp[num] <- signif(tab[num], 3)
print(tab_imp, row.names = FALSE)
cat(sprintf("\nModa do BCB no nosso intervalo: %d de %d | bem identificados: %d de 22 com priori uniforme",
            sum(tab$compativel), nrow(tab), sum(tab$bem_ident)))
cat(sprintf(" | compatíveis e bem identificados: %d\n", sum(tab$compativel & tab$bem_ident)))
cat(sprintf("Correlação a posteriori entre θ1 e θ2: %.2f\n", cor(amostras[, "theta_1"], amostras[, "theta_2"])))

# ---- Tabela 4.2: combinações de parâmetros ----
combinacoes <- function(th) cbind(
  lambda            = 1 - th[, "theta_1"] - th[, "theta_2"],
  peso_exp_phillips = 1 - th[, "alpha_1L"] - th[, "alpha_1I"],
  persistencia_IS   = 1 - (1 - th[, "beta1"]) * (1 - th[, "beta5"]),
  efeito_lp_juro    = (th[, "beta2"] / 4) / (1 - th[, "beta1"]),
  peso_meta         = 1 - th[, "phi_1"] - th[, "phi_2"] - th[, "phi_3"])
ca <- combinacoes(amostras)
tab_comb <- data.frame(combinacao = colnames(ca), mediana = apply(ca, 2, median),
                       q05 = apply(ca, 2, quantile, .05), q95 = apply(ca, 2, quantile, .95),
                       bcb = as.numeric(combinacoes(t(par_canonico))))
print(tab_comb, digits = 3, row.names = FALSE)

# ---- Hiato e taxa neutra com bandas de 90% (Figuras 4.4 e 4.5) ----
# Faixa escura: incerteza dos parâmetros (500 sorteios da posteriori).
# Faixa clara: parâmetros e filtragem (soma a variância do estado suavizado).
set.seed(42)
n_draws  <- 500
sorteios <- sample(nrow(amostras), n_draws)
datas    <- dados_obs$trimestre[efetivo]
tend     <- dados_exog$juro_real_hp[efetivo]
h_par <- h_tot <- r_par <- r_tot <- matrix(NA_real_, n_draws, length(datas))
for (i in seq_len(n_draws)) {
  kf <- KFS(construir_kfas_lp(combinar_par(amostras[sorteios[i], ]), dados_exog, dados_obs, s_perfil)$modelo,
            smoothing = "state")
  h_par[i, ] <- kf$alphahat[, idx$h]
  h_tot[i, ] <- h_par[i, ] + sqrt(pmax(kf$V[idx$h, idx$h, ], 0)) * rnorm(length(datas))
  r_par[i, ] <- tend + kf$alphahat[, idx$r_is]
  r_tot[i, ] <- r_par[i, ] + sqrt(pmax(kf$V[idx$r_is, idx$r_is, ], 0)) * rnorm(length(datas))
}
faixas <- function(x_par, x_tot) data.frame(
  trimestre = datas, mediana = apply(x_par, 2, median),
  q05_par = apply(x_par, 2, quantile, .05), q95_par = apply(x_par, 2, quantile, .95),
  q05_tot = apply(x_tot, 2, quantile, .05), q95_tot = apply(x_tot, 2, quantile, .95))
serie_hiato  <- faixas(h_par, h_tot)
serie_neutra <- faixas(r_par, r_tot)
serie_neutra$juro_ex_ante <- dados_exog$juro_real_exante[efetivo]
serie_neutra$tendencia    <- tend
# Largura média das bandas e parcela da variância que vem da filtragem
for (nm in c("hiato", "neutra")) {
  sr <- if (nm == "hiato") serie_hiato else serie_neutra
  lp <- mean(sr$q95_par - sr$q05_par); lt <- mean(sr$q95_tot - sr$q05_tot)
  cat(sprintf("%-6s banda de 90%%: parâmetros %.2f | total %.2f p.p. | filtragem %.0f%% da variância\n",
              nm, lp, lt, 100 * (1 - (lp / lt)^2)))
}

# Projeção do modelo contra a Focus e o IPCA realizado
# Estado filtrado, no modo
r_modo <- construir_kfas_lp(combinar_par(theta_modo), dados_exog, dados_obs, s_perfil)
at     <- KFS(r_modo$modelo, filtering = "state", smoothing = "none")$att
ipca_seguinte <- dplyr::lead(acum4(painel_anual$ipca_cheio_sa), 4)
serie_projecao <- data.frame(
  trimestre = datas,
  projecao  = r_modo$proj$d[efetivo] + rowSums(r_modo$proj$G[efetivo, ] * at[, 1:m]),
  focus     = dados_obs$pi_focus[efetivo],
  ipca_seguinte = ipca_seguinte[match(datas, painel_anual$trimestre)])
cat(sprintf("Projeção x Focus: correlação %.2f\n", cor(serie_projecao$projecao, serie_projecao$focus)))
reqm <- function(x) sqrt(mean((x - serie_projecao$ipca_seguinte)^2, na.rm = TRUE))   # raiz do erro quadrático médio
cat(sprintf("Erro contra o IPCA dos 4 trimestres seguintes (REQM): projeção %.2f | Focus %.2f p.p.\n",
            reqm(serie_projecao$projecao), reqm(serie_projecao$focus)))

#Figuras 4.1 a 4.5 (PDF vetorial; título e fonte ficam no LaTeX)
FONTE <- "Times New Roman"; TAM <- 10
AZUL <- "#2C6FD1"; CINZA <- "grey35"; CLARO <- "grey70"; FAIXA <- "grey92"
virgula <- function(x) format(x, decimal.mark = ",", trim = TRUE)
salvar <- function(g, arquivo, altura) {
  ggsave(arquivo, g, width = 16, height = altura, units = "cm", device = cairo_pdf,
         family = FONTE, symbolfamily = cairoSymbolFont(FONTE, usePUA = FALSE))
  cat("Gravado:", arquivo, "\n")
}
tema <- theme_minimal(base_size = TAM, base_family = FONTE) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.title = element_blank(), plot.margin = margin(4, 8, 4, 4))
anos   <- scale_x_date(breaks = as.Date(paste0(seq(2004, 2020, 4), "-01-01")), date_labels = "%Y")
bandas <- scale_fill_manual(values = c("Parâmetros e filtragem" = adjustcolor(AZUL, 0.15),
                                       "Incerteza dos parâmetros" = adjustcolor(AZUL, 0.38)))

#nosso intervalo e a moda do BCB na escala do intervalo do BCB (0 = limite inferior, 1 = superior)
rotulo <- c(alpha_1L = "alpha[1*L]", alpha_1I = "alpha[1*I]", alpha_2 = "alpha[2]", alpha_3 = "alpha[3]",
            alpha_4 = "alpha[4]", alpha_5 = "alpha[5]", alpha_6 = "alpha[6]", beta1 = "beta[1]",
            beta2 = "beta[2]", beta3 = "beta[3]", beta4 = "beta[4]", beta5 = "beta[5]",
            theta_1 = "theta[1]", theta_2 = "theta[2]", theta_3 = "theta[3]", phi_1 = "varphi[1]",
            phi_2 = "varphi[2]", phi_3 = "varphi[3]", delta = "delta", gamma_nuci = "gamma[nuci]",
            gamma_emp = "gamma[emp]", gamma_caged = "gamma[caged]", sigma_h = "sigma^h")
bloco <- rep(c("Phillips", "IS", "Taylor", "Expectativas", "Paridade", "Observação"), c(7, 5, 3, 3, 1, 4))
largura <- IC_BCB["sup", ] - IC_BCB["inf", ]
f1 <- data.frame(rotulo = ifelse(tab$compativel, rotulo, paste0("bold(", rotulo, ")")),
                 bloco  = factor(bloco, unique(bloco)),
                 lo     = (tab$ic_inf   - IC_BCB["inf", ]) / largura,
                 hi     = (tab$ic_sup   - IC_BCB["inf", ]) / largura,
                 media  = (tab$media    - IC_BCB["inf", ]) / largura,
                 moda   = (tab$moda_bcb - IC_BCB["inf", ]) / largura)
f1$rotulo <- factor(f1$rotulo, rev(f1$rotulo))
g1 <- ggplot(f1, aes(y = rotulo)) +
  geom_rect(data = data.frame(xmin = 0, xmax = 1), inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = "Intervalo de 90% do BCB")) +
  geom_point(aes(x = moda, shape = "Moda do BCB"), size = 3.2) +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = "Nosso intervalo de 90% e média"), width = 0,
                orientation = "y", linewidth = 0.7) +
  geom_point(aes(x = media, colour = "Nosso intervalo de 90% e média"), size = 1.6) +
  facet_grid(bloco ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = FAIXA) + scale_shape_manual(values = 124) + scale_colour_manual(values = AZUL) +
  guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2),
         colour = guide_legend(order = 3, override.aes = list(linetype = 1, shape = 16))) +
  scale_y_discrete(labels = function(l) parse(text = l)) +
  scale_x_continuous(breaks = -1:2, labels = virgula) +
  labs(x = "Posição na escala do intervalo de 90% do BCB (0 = limite inferior; 1 = superior)",
       y = NULL, caption = "Em negrito: parâmetros em que a moda do BCB fica fora do nosso intervalo de 90%.") +
  tema + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
               panel.grid.major.y = element_blank(), plot.caption = element_text(hjust = 0, size = TAM - 1))
salvar(g1, "fig4_1_intervalos.pdf", 16)

#posteriori conjunta de β1 e β5, com a curva de persistência constante
set.seed(7)
pts   <- as.data.frame(amostras[sample(nrow(amostras), min(3000, nrow(amostras))), c("beta1", "beta5")])
pers  <- median(1 - (1 - amostras[, "beta1"]) * (1 - amostras[, "beta5"]))
curva <- data.frame(beta1 = seq(0.5, 1 - (1 - pers) / 0.55, length.out = 200))
curva$beta5 <- 1 - (1 - pers) / (1 - curva$beta1)
g2 <- ggplot(pts, aes(beta1, beta5)) +
  geom_point(colour = CLARO, size = 0.4, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = CINZA, linewidth = 0.3) +
  geom_line(data = curva, colour = AZUL, linewidth = 0.9) +
  annotate("point", x = 0.85, y = 0.84, size = 2.2) +
  annotate("text", x = 0.86, y = 0.62, label = "Moda do BCB (0,85; 0,84)", hjust = 0,
           size = TAM / 3.2, family = FONTE) +
  annotate("segment", x = 0.85, y = 0.835, xend = 0.865, yend = 0.635, linewidth = 0.3) +
  annotate("text", x = 0.51, y = 0.985, label = paste0("persistência = ", virgula(round(pers, 3))),
           hjust = 0, colour = AZUL, size = TAM / 3.2, family = FONTE) +
  annotate("text", x = 0.995, y = 0.93, label = "beta[1] == beta[5]", parse = TRUE, hjust = 1,
           colour = CINZA, size = TAM / 3.4) +
  coord_cartesian(xlim = c(0.5, 1), ylim = c(0.55, 1), expand = FALSE) +
  scale_x_continuous(labels = virgula) + scale_y_continuous(labels = virgula) +
  labs(x = expression(beta[1] * ", inércia própria do hiato"), y = expression(beta[5] * ", persistência do choque")) +
  tema
salvar(g2, "fig4_2_beta1_beta5.pdf", 11)

#projeção do modelo, Focus e IPCA dos quatro trimestres seguintes
f3 <- rbind(
  data.frame(trimestre = datas, valor = serie_projecao$projecao, serie = "Projeção do modelo"),
  data.frame(trimestre = datas, valor = serie_projecao$focus, serie = "Expectativa Focus"),
  data.frame(trimestre = datas, valor = serie_projecao$ipca_seguinte,
             serie = "IPCA efetivo nos quatro trimestres seguintes"))
f3$serie <- factor(f3$serie, unique(f3$serie))
g3 <- ggplot(f3, aes(trimestre, valor, colour = serie, linetype = serie, linewidth = serie)) +
  geom_line() +
  scale_colour_manual(values = c(AZUL, CINZA, CLARO)) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  scale_linewidth_manual(values = c(0.9, 0.5, 0.5)) +
  anos + scale_y_continuous(labels = virgula, limits = c(0, 12.5)) +
  labs(x = NULL, y = "% em quatro trimestres") + tema
salvar(g3, "fig4_3_projecao_focus.pdf", 8.5)

#hiato, com as crises sombreadas e alguns valores marcados
crises <- data.frame(ini = as.Date(c("2008-08-15", "2019-11-15")), fim = as.Date(c("2009-11-15", "2020-11-15")))
marcas <- serie_hiato[serie_hiato$trimestre %in%
                        as.Date(c("2009-01-01", "2013-10-01", "2016-10-01", "2020-04-01", "2023-10-01")), ]
marcas$rot <- paste0(year(marcas$trimestre), "T", quarter(marcas$trimestre), ": ",
                     ifelse(marcas$mediana > 0, "+", ""), virgula(round(marcas$mediana, 1)), "%")
g4 <- ggplot(serie_hiato, aes(trimestre)) +
  geom_rect(data = crises, aes(xmin = ini, xmax = fim, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = FAIXA) +
  geom_hline(yintercept = 0, colour = CINZA, linewidth = 0.3) +
  geom_ribbon(aes(ymin = q05_tot, ymax = q95_tot, fill = "Parâmetros e filtragem")) +
  geom_ribbon(aes(ymin = q05_par, ymax = q95_par, fill = "Incerteza dos parâmetros")) +
  geom_line(aes(y = mediana, colour = "Mediana"), linewidth = 0.8) +
  geom_point(data = marcas, aes(y = mediana), colour = AZUL, size = 1.4) +
  geom_text(data = marcas, aes(y = mediana, label = rot), size = TAM / 3.4, family = FONTE,
            vjust = ifelse(marcas$mediana > 0, -1.2, 1.9),
            hjust = ifelse(marcas$trimestre > as.Date("2022-01-01"), 1, 0.5)) +
  bandas + scale_colour_manual(values = c(Mediana = AZUL)) +
  anos + scale_y_continuous(labels = virgula, breaks = seq(-8, 2, 2)) +
  labs(x = NULL, y = "% do produto potencial") + tema
salvar(g4, "fig4_4_hiato.pdf", 9)

#taxa neutra, juro real ex-ante e sua tendência HP
g5 <- ggplot(serie_neutra, aes(trimestre)) +
  geom_hline(yintercept = 0, colour = CINZA, linewidth = 0.3) +
  geom_ribbon(aes(ymin = q05_tot, ymax = q95_tot, fill = "Parâmetros e filtragem")) +
  geom_ribbon(aes(ymin = q05_par, ymax = q95_par, fill = "Incerteza dos parâmetros")) +
  geom_line(aes(y = juro_ex_ante, colour = "Juro real ex-ante"), linewidth = 0.4) +
  geom_line(aes(y = tendencia, colour = "Tendência HP do juro real"), linewidth = 0.5, linetype = "dashed") +
  geom_line(aes(y = mediana, colour = "Taxa neutra (mediana)"), linewidth = 0.8) +
  bandas +
  scale_colour_manual(values = c("Taxa neutra (mediana)" = AZUL, "Tendência HP do juro real" = CINZA,
                                 "Juro real ex-ante" = CLARO)) +
  anos + scale_y_continuous(labels = virgula, breaks = seq(-2, 14, 2)) +
  labs(x = NULL, y = "% ao ano") +
  tema + theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt")) +
  guides(fill = guide_legend(order = 2), colour = guide_legend(order = 1))
salvar(g5, "fig4_5_neutro.pdf", 10.5)