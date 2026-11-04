library(tidyr)

N_DRAWS  <- 300                                    # sorteios da posteriori para as bandas
H_IRF    <- 17                                     # trimestres 0 a 16, como no boxe
W_LIVRES <- mean(tail(dados_exog$peso_livres, 4))  # peso dos preços livres em 2023
RHO_ICBR <- ar_coefs$icbr_desvio_an$rho            # AR(1) do IC-Br na projeção

# Efeitos máximos do boxe (modelo agregado)
ALVOS <- c(juros = -0.27, cambio_ipca = 0.96, cambio_livres = 0.72, cambio_adm = 1.65, hiato = 0.49)

pico     <- function(x) x[which.max(abs(x))]                     # efeito máximo (com sinal)
tri_pico <- function(df, v) df$trimestre[which.max(abs(df[[v]]))]
imprimir <- function(df) { num <- sapply(df, is.numeric); df[num] <- round(df[num], 3); print(df, row.names = FALSE) }

# Em desvios do equilíbrio, a expectativa do modelo é linear nas variáveis de t:
# E_t π = coef'x. Os coeficientes saem da mesma recursão de proj_afim (1_modelo.R),
# sem os termos de meta e das exógenas, que são nulos em desvios.
ENTRADAS <- c("h", "s", "rr", "pe", "piL", "ipca0", "ipca1", "ipca2", "ipca3",
              "adm_proj", "esel", "de", "ps")

coef_operador <- function(par, w = W_LIVRES, rho = RHO_ICBR, H = 4) {
  pp <- as.list(par); cexp <- 1 - pp$alpha_1L - pp$alpha_1I
  recursao <- function(x, E) {
    h <- x[["h"]]; s <- x[["s"]]; pe <- x[["pe"]]; pl <- x[["piL"]]
    ult <- c(x[["ipca0"]], x[["ipca1"]], x[["ipca2"]], x[["ipca3"]]); soma <- 0
    for (j in seq_len(H)) {
      ma4  <- mean(ult)
      h_n  <- pp$beta1 * h + (pp$beta2 / 4) * (pe + x[["rr"]]) + pp$beta5 * s - (pp$beta2 / 4) * x[["esel"]]
      s    <- pp$beta5 * s
      pe_n <- pp$phi_1 * pe + pp$phi_2 * E + pp$phi_3 * ma4
      pl   <- pp$alpha_1L * pl + pp$alpha_1I * ma4 + cexp * pe_n + pp$alpha_4 * h_n +
        pp$alpha_2 * rho^j * x[["ps"]] + pp$alpha_3 * (if (j == 1) x[["de"]] else 0)
      ip   <- w * pl + (1 - w) * x[["adm_proj"]]
      soma <- soma + ip; ult <- c(ip, ult[1:3]); h <- h_n; pe <- pe_n
    }
    soma / H
  }
  zero <- setNames(numeric(length(ENTRADAS)), ENTRADAS)
  b <- recursao(zero, 1)                         # peso da própria E_t na recursão
  if (1 - b <= 0.05) return(NULL)                # sem solução estável (mesma guarda de proj_afim)
  sapply(ENTRADAS, function(k) recursao(replace(zero, k, 1), 0) / (1 - b))
}

# Choques (t0 = trimestre do choque):
#   juros_regra  Selic +1 p.p. por quatro trimestres, depois a regra de Taylor
#   juros        trajetória fixa da Selic: +1 p.p. por quatro trimestres, depois decai a 0,8
#   cambio       depreciação permanente de 10%
#   hiato        hiato +1 p.p. por um trimestre
# Administrados: indexados ao IPCA dos quatro trimestres anteriores; uma fração
# comb_adm (combustíveis, seção 5.2) passa a seguir o câmbio.
# canal_cambio = FALSE fixa o câmbio; canal_exp = FALSE fixa as expectativas.
simular_irf <- function(par, tipo = c("juros_regra", "juros", "cambio", "hiato"),
                        canal_cambio = TRUE, canal_exp = TRUE, comb_adm = 0,
                        H = H_IRF, w = W_LIVRES, max_it = 1000, tol = 1e-10) {
  tipo <- match.arg(tipo)
  p <- as.list(par)
  cexp <- 1 - p$alpha_1L - p$alpha_1I
  lam  <- 1 - p$theta_1 - p$theta_2
  coef_op <- coef_operador(par, w = w)
  if (canal_exp && is.null(coef_op)) return(NULL)
  
  t0 <- 6; n <- t0 + H + 8; zz <- function() numeric(n)
  imp_h <- zz(); de_exo <- rep(NA_real_, n); i_fix <- rep(NA_real_, n)
  if (tipo == "hiato")  imp_h[t0] <- 1
  if (tipo == "cambio") { de_exo[t0:n] <- 0; de_exo[t0] <- 40 }   # 10% no trimestre = 40 anualizado
  if (tipo %in% c("juros_regra", "juros")) i_fix[t0:(t0 + 3)] <- 1
  if (tipo == "juros") for (t in (t0 + 4):n) i_fix[t] <- 0.8 * i_fix[t - 1]
  
  # A Selic esperada depende da trajetória futura da Selic: itera até ela não mudar
  i_g <- ipca_g <- adm_g <- zz()
  for (it in seq_len(max_it)) {
    s <- h <- pe <- piL <- i <- de <- esel <- adm <- ipca <- E <- zz()
    for (t in t0:n) {
      s[t] <- p$beta5 * s[t - 1]
      h[t] <- p$beta1 * h[t - 1] + (p$beta2 / 4) * (pe[t - 1] - esel[t - 1]) + s[t] + imp_h[t]
      ma4  <- mean(ipca[(t - 4):(t - 1)])
      fut  <- function(k) i_g[min(t + k, n)]
      
      # Todas as variáveis de t, dada a expectativa Focus π^e_t
      bloco <- function(pe_t) {
        i_t  <- if (!is.na(i_fix[t])) i_fix[t] else
          p$theta_1 * i[t - 1] + p$theta_2 * i[t - 2] + lam * p$theta_3 * pe_t
        de_t <- if (!canal_cambio) 0 else if (!is.na(de_exo[t])) de_exo[t] else -p$delta * (i_t - i[t - 1])
        es_t <- (0.5 * i_t + fut(1) + fut(2) + fut(3) + 0.5 * fut(4)) / 4
        ad_t <- (1 - comb_adm) * ma4 + comb_adm * de_t
        pl_t <- p$alpha_1L * piL[t - 1] + p$alpha_1I * ma4 + cexp * pe_t + p$alpha_2 * de_t +
          p$alpha_3 * de[t - 1] + p$alpha_4 * h[t]
        ip_t <- w * pl_t + (1 - w) * ad_t
        E_t  <- 0
        if (canal_exp) {
          x <- c(h = h[t], s = s[t], rr = 0, pe = pe_t, piL = pl_t, ipca0 = ip_t,
                 ipca1 = ipca[t - 1], ipca2 = ipca[t - 2], ipca3 = ipca[t - 3],
                 adm_proj = mean(c(ip_t, ipca[t - 1], ipca[t - 2], ipca[t - 3])),
                 esel = es_t, de = de_t, ps = de_t)
          E_t <- sum(coef_op[ENTRADAS] * x[ENTRADAS])
        }
        list(i = i_t, de = de_t, esel = es_t, adm = ad_t, piL = pl_t, ipca = ip_t, E = E_t,
             pe_impl = p$phi_1 * pe[t - 1] + p$phi_2 * E_t + p$phi_3 * ma4)
      }
      
      # π^e_t resolve π^e = F(π^e), com F afim: solução exata com duas avaliações
      if (canal_exp) {
        b0 <- bloco(0); b1 <- bloco(1)
        pe[t] <- b0$pe_impl / (1 - (b1$pe_impl - b0$pe_impl))
      }
      r <- bloco(pe[t])
      i[t] <- r$i; de[t] <- r$de; esel[t] <- r$esel; adm[t] <- r$adm
      piL[t] <- r$piL; ipca[t] <- r$ipca; E[t] <- r$E
    }
    dif <- max(abs(i - i_g), abs(ipca - ipca_g), abs(adm - adm_g))
    i_g <- i; ipca_g <- ipca; adm_g <- adm
    if (dif < tol) break
  }
  if (dif >= tol) warning("simular_irf(", tipo, "): não convergiu")
  
  acum4 <- function(v) sapply(seq_len(n), function(t) if (t < 4) NA_real_ else sum(v[(t - 3):t]) / 4)
  k <- t0:(t0 + H - 1)
  data.frame(trimestre    = 0:(H - 1),
             selic        = i[k],
             hiato        = h[k],
             cambio_acum  = cumsum(de / 4)[k],        # câmbio, % acumulada desde o choque
             ipca_acum4   = acum4(ipca)[k],           # inflação acumulada em quatro trimestres
             livres_acum4 = acum4(piL)[k],
             adm_acum4    = acum4(adm)[k])
}

#Respostas no nosso modo, com bandas, e nas modas do BCB
p_modo <- combinar_par(theta_modo)
p_bcb  <- combinar_par(par_bcb_base)

CHOQUES <- c(juros = "juros_regra", cambio = "cambio", hiato = "hiato")
VARS    <- c("selic", "hiato", "cambio_acum", "ipca_acum4", "livres_acum4", "adm_acum4")

set.seed(13)
sorteios <- sample(nrow(amostras), N_DRAWS)
irf <- list(); picos_post <- list()
for (nm in names(CHOQUES)) {
  r <- simular_irf(p_modo, CHOQUES[[nm]])
  b <- simular_irf(p_bcb,  CHOQUES[[nm]])
  arr <- array(NA_real_, c(N_DRAWS, H_IRF, length(VARS)), dimnames = list(NULL, NULL, VARS))
  for (d in seq_len(N_DRAWS)) {
    rd <- simular_irf(combinar_par(amostras[sorteios[d], ]), CHOQUES[[nm]])
    if (!is.null(rd)) arr[d, , ] <- as.matrix(rd[, VARS])
  }
  for (v in VARS) {
    r[[paste0(v, "_q05")]] <- apply(arr[, , v], 2, quantile, .05, na.rm = TRUE)
    r[[paste0(v, "_q95")]] <- apply(arr[, , v], 2, quantile, .95, na.rm = TRUE)
    r[[paste0(v, "_bcb")]] <- b[[v]]
  }
  irf[[nm]] <- r
  picos_post[[nm]] <- sapply(VARS, function(v) apply(arr[, , v], 1, pico))
  cat(sprintf("%-7s pronto (%d sorteios sem solução estável)\n", nm, sum(is.na(arr[, 1, 1]))))
}

# Efeitos máximos contra o boxe
linha <- function(nm, var, alvo = NA) {
  df <- irf[[nm]]; pp <- picos_post[[nm]][, var]
  data.frame(choque = nm, variavel = var,
             modo = pico(df[[var]]), tri_modo = tri_pico(df, var),
             post_mediana = median(pp, na.rm = TRUE),
             post_q05 = quantile(pp, .05, na.rm = TRUE), post_q95 = quantile(pp, .95, na.rm = TRUE),
             modas_bcb = pico(df[[paste0(var, "_bcb")]]), tri_bcb = tri_pico(df, paste0(var, "_bcb")),
             boxe = alvo, row.names = NULL)
}
tab_irf <- rbind(
  linha("juros",  "ipca_acum4",   ALVOS[["juros"]]),
  linha("juros",  "livres_acum4"),
  linha("juros",  "hiato"),
  linha("cambio", "ipca_acum4",   ALVOS[["cambio_ipca"]]),
  linha("cambio", "livres_acum4", ALVOS[["cambio_livres"]]),
  linha("cambio", "adm_acum4",    ALVOS[["cambio_adm"]]),
  linha("hiato",  "ipca_acum4",   ALVOS[["hiato"]]))
cat("\n=== Efeitos máximos (p.p.) ===\n"); imprimir(tab_irf)


# 3. Canais de transmissão: trajetória fixa da Selic
CANAIS <- list(`Expectativa e câmbio livres`    = c(cambio = TRUE,  exp = TRUE),
               `Expectativa livre, câmbio fixo` = c(cambio = FALSE, exp = TRUE),
               `Expectativa e câmbio fixos`     = c(cambio = FALSE, exp = FALSE))
canais <- data.frame(trimestre = 0:(H_IRF - 1))
for (nm in names(CANAIS)) {
  r <- simular_irf(p_modo, "juros", canal_cambio = CANAIS[[nm]][["cambio"]], canal_exp = CANAIS[[nm]][["exp"]])
  canais[[nm]] <- r$ipca_acum4
  cat(sprintf("%-32s efeito máximo %+.3f no trimestre %d | no 4º trimestre %+.3f\n",
              nm, pico(r$ipca_acum4), tri_pico(r, "ipca_acum4"), r$ipca_acum4[5]))
}

#De onde vêm as diferenças
# Troca de blocos: parte das modas do BCB e troca um bloco por vez pelo nosso
BLOCOS <- list(Taylor = c("theta_1", "theta_2", "theta_3"), IS = c("beta1", "beta2", "beta5"),
               Phillips = c("alpha_1L", "alpha_1I", "alpha_2", "alpha_3", "alpha_4"),
               Expectativas = c("phi_1", "phi_2", "phi_3"), Paridade = "delta")
resumo_juros <- function(par, rotulo) {
  r <- simular_irf(par, "juros_regra")
  data.frame(parametros = rotulo, ipca = pico(r$ipca_acum4), tri = tri_pico(r, "ipca_acum4"),
             hiato = pico(r$hiato), selic_6tri = r$selic[7])
}
blocos <- rbind(
  resumo_juros(p_bcb, "Modas do BCB"),
  do.call(rbind, lapply(names(BLOCOS), function(b) {
    pp <- p_bcb; pp[BLOCOS[[b]]] <- p_modo[BLOCOS[[b]]]
    resumo_juros(pp, paste("BCB com o nosso bloco:", b))
  })),
  resumo_juros(p_modo, "Nosso modo"))
cat("\n=== Choque de juros: troca de blocos de parâmetros (seção 5.2) ===\n"); imprimir(blocos)

# Combustíveis: um quarto dos administrados passa a seguir o câmbio
casos <- expand.grid(fracao_cambio = c(0, 0.25), parametros = c("Nosso modo", "Modas do BCB"),
                     stringsAsFactors = FALSE)
combustiveis <- do.call(rbind, lapply(seq_len(nrow(casos)), function(k) {
  par <- if (casos$parametros[k] == "Nosso modo") p_modo else p_bcb
  f   <- casos$fracao_cambio[k]
  rj  <- simular_irf(par, "juros_regra", comb_adm = f)
  rc  <- simular_irf(par, "cambio", comb_adm = f)
  data.frame(casos[k, ], juros_ipca = pico(rj$ipca_acum4), juros_tri = tri_pico(rj, "ipca_acum4"),
             cambio_ipca = pico(rc$ipca_acum4), cambio_adm = pico(rc$adm_acum4))
}))
cat("\n=== Administrados: indexação pura (0) e 25% seguindo o câmbio ===\n"); imprimir(combustiveis)


#graficos
LARANJA <- "#C2571A"; VERDE <- "#3B8E4F"
EIXO_X <- scale_x_continuous(breaks = seq(0, 16, 2))
tema <- theme_minimal(base_size = TAM, base_family = FONTE) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(),
        strip.text = element_text(hjust = 0, size = TAM), plot.margin = margin(4, 8, 4, 4))
linha_zero <- geom_hline(yintercept = 0, colour = CINZA, linewidth = 0.3)

SERIES <- c("Nosso modo", "Nosso modelo nas modas do BCB")
painel_ab <- function(df, va, rot_a, alvo) {
  rot_b  <- "B. IPCA acumulado em quatro trimestres"
  niveis <- c(rot_a, rot_b)
  lng <- function(v, painel) data.frame(trimestre = df$trimestre, painel = painel,
                                        lo = df[[paste0(v, "_q05")]], hi = df[[paste0(v, "_q95")]],
                                        modo = df[[v]], bcb = df[[paste0(v, "_bcb")]])
  d <- rbind(lng(va, rot_a), lng("ipca_acum4", rot_b)); d$painel <- factor(d$painel, niveis)
  dl <- pivot_longer(d, c(modo, bcb), names_to = "serie", values_to = "valor")
  dl$serie <- factor(ifelse(dl$serie == "modo", SERIES[1], SERIES[2]), SERIES)
  rf <- data.frame(painel = factor(rot_b, niveis), y = alvo,
                   rot = paste0("boxe: ", ifelse(alvo > 0, "+", "−"), virgula(abs(alvo))))
  ggplot(d, aes(trimestre)) +
    linha_zero +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = "Intervalo de 90% (posteriori)")) +
    geom_hline(data = rf, aes(yintercept = y), linetype = "dotted", linewidth = 0.5) +
    geom_text(data = rf, aes(x = 0, y = y, label = rot), hjust = 0, vjust = ifelse(alvo > 0, -0.5, 1.5),
              size = TAM / 3.4, family = FONTE) +
    geom_line(data = dl, aes(y = valor, colour = serie, linetype = serie, linewidth = serie)) +
    facet_wrap(~ painel, scales = "free_y") +
    scale_fill_manual(values = adjustcolor(AZUL, 0.18)) +
    scale_colour_manual(values = c(AZUL, CINZA)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_linewidth_manual(values = c(0.9, 0.6)) +
    EIXO_X + scale_y_continuous(labels = virgula) +
    labs(x = "Trimestres", y = "p.p.") + tema +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),
           linewidth = guide_legend(order = 1), fill = guide_legend(order = 2))
}
salvar(painel_ab(irf$juros, "selic", "A. Taxa Selic", ALVOS[["juros"]]), "fig5_1_juros.pdf", 8.5)
salvar(painel_ab(irf$hiato, "hiato", "A. Hiato do produto", ALVOS[["hiato"]]), "fig5_3_hiato.pdf", 8.5)

#Câmbio. A: câmbio acumulado; B: IPCA (com banda), livres e administrados,
cb <- irf$cambio
NIVEIS <- c("IPCA", "Livres", "Administrados")
dB <- rbind(
  data.frame(trimestre = cb$trimestre, serie = "IPCA", valor = cb$ipca_acum4, lo = cb$ipca_acum4_q05, hi = cb$ipca_acum4_q95),
  data.frame(trimestre = cb$trimestre, serie = "Livres", valor = cb$livres_acum4, lo = NA, hi = NA),
  data.frame(trimestre = cb$trimestre, serie = "Administrados", valor = cb$adm_acum4, lo = NA, hi = NA))
dB$serie <- factor(dB$serie, NIVEIS)
refB <- data.frame(serie = factor(NIVEIS, NIVEIS), y = unname(ALVOS[c("cambio_ipca", "cambio_livres", "cambio_adm")]))
gA <- ggplot(data.frame(trimestre = cb$trimestre, valor = cb$cambio_acum), aes(trimestre, valor)) +
  linha_zero + geom_line(colour = "black", linewidth = 0.9) +
  EIXO_X + scale_y_continuous(labels = virgula, limits = c(0, 12)) +
  labs(x = "Trimestres", y = "%", subtitle = "A. Câmbio, variação acumulada") + tema
gB <- ggplot(dB, aes(trimestre, valor, colour = serie)) +
  linha_zero +
  geom_ribbon(data = dB[dB$serie == "IPCA", ], aes(ymin = lo, ymax = hi), fill = adjustcolor(AZUL, 0.18), colour = NA) +
  geom_hline(data = refB, aes(yintercept = y, colour = serie), linetype = "dotted", linewidth = 0.5, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(IPCA = AZUL, Livres = LARANJA, Administrados = VERDE)) +
  EIXO_X + scale_y_continuous(labels = virgula) +
  labs(x = "Trimestres", y = "p.p.", subtitle = "B. Inflação acumulada em quatro trimestres") + tema
salvar(patchwork::wrap_plots(gA, gB, widths = c(1, 1.25)) +
         patchwork::plot_layout(guides = "collect") & theme(legend.position = "top"),
       "fig5_2_cambio.pdf", 8.5)

#Canais de transmissão, no nosso modo
dc <- pivot_longer(canais, -trimestre, names_to = "canal", values_to = "valor")
dc$canal <- factor(dc$canal, names(CANAIS))
g54 <- ggplot(dc, aes(trimestre, valor, colour = canal, linetype = canal)) +
  linha_zero + geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c(AZUL, LARANJA, CINZA)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  EIXO_X + scale_y_continuous(labels = virgula) +
  labs(x = "Trimestres", y = "p.p.") + tema
salvar(g54, "fig5_4_canais.pdf", 8)