library(tidyr)

N_DRAWS  <- 300                                    # sorteios da posteriori para as bandas
H_IRF    <- 17                                     # trimestres 0 a 16, como no boxe
W_LIVRES <- mean(tail(dados_exog$peso_livres, 4))  # peso dos preços livres em 2023
RHO_ICBR <- ar_coefs$icbr_desvio_an$rho            # AR(1) do IC-Br na projeção

pico     <- function(x) x[which.max(abs(x))]                     # efeito máximo (com sinal)
tri_pico <- function(df, v) df$trimestre[which.max(abs(df[[v]]))]
reqm     <- function(a, b) sqrt(mean((a - b)^2))                 # raiz do erro quadrático médio
imprimir <- function(df) { num <- sapply(df, is.numeric); df[num] <- round(df[num], 3); print(df, row.names = FALSE) }

# Respostas do boxe (modelo agregado), do anexo estatístico do RI de jun/2024.
# Convenção do boxe, adotada aqui: o trimestre 0 é o anterior ao choque e o choque
# ocorre no trimestre 1.
ler_irf <- function(nome, col = 2) {
  a <- aba_anexo(nome); a <- a[como_num(a[[1]]) %in% 0:16, ]
  como_num(a[[col]])[order(como_num(a[[1]]))]
}
BOXE <- list(juros  = list(selic = ler_irf("C2 Boxe3 Graf 1A"), ipca_acum4 = ler_irf("C2 Boxe3 Graf 1B")),
             cambio = list(cambio_acum = ler_irf("C2 Boxe3 Graf 2A"), ipca_acum4 = ler_irf("C2 Boxe3 Graf 2B")),
             hiato  = list(hiato = ler_irf("C2 Boxe3 Graf 3A"), ipca_acum4 = ler_irf("C2 Boxe3 Graf 3B")))
BOXE_CANAIS <- data.frame(trimestre = 0:16,
                          `Expectativa e câmbio livres`    = ler_irf("C2 Boxe3 Graf 4B", 2),
                          `Expectativa livre, câmbio fixo` = ler_irf("C2 Boxe3 Graf 4B", 3),
                          `Expectativa e câmbio fixos`     = ler_irf("C2 Boxe3 Graf 4B", 4), check.names = FALSE)
# Livres e administrados só aparecem no texto do boxe, como efeitos máximos
BOXE_TEXTO <- c(cambio_livres = 0.72, cambio_adm = 1.65)

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
# comb_adm passa a seguir o câmbio no mesmo trimestre (ver FRACAO_ADM abaixo).
# canal_cambio = FALSE fixa o câmbio; canal_exp = FALSE fixa as expectativas.
# Saída na convenção do boxe: trimestre 0 antes do choque, choque no trimestre 1.
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
  k <- (t0 - 1):(t0 + H - 2)                    # do trimestre anterior ao choque em diante
  data.frame(trimestre    = 0:(H - 1),
             selic        = i[k],
             hiato        = h[k],
             cambio_acum  = cumsum(de / 4)[k],        # câmbio, % acumulada desde o choque
             ipca_acum4   = acum4(ipca)[k],           # inflação acumulada em quatro trimestres
             livres_acum4 = acum4(piL)[k],
             adm_acum4    = acum4(adm)[k])
}

p_modo <- combinar_par(theta_modo)
p_bcb  <- combinar_par(par_bcb_base)

# ---- Administrados ----
# Versão principal (bandas, figuras e Tabela de efeitos máximos): os administrados
# seguem o IPCA dos quatro trimestres anteriores, a mesma regra da projeção usada na
# estimação (proj_afim).
# O boxe informa que as simulações do BCB incorporam o seu modelo de administrados
# (Brasil, 2017b), em que parte dos itens, como os combustíveis, acompanha o câmbio.
# Para medir quanto isso explica da diferença em relação ao boxe, FRACAO_ADM é a fração
# dos administrados que acompanha o câmbio no mesmo trimestre, calibrada para que, com
# as modas do BCB, o efeito máximo da depreciação sobre os administrados seja o do boxe
# (1,65). Ela só aparece nas comparações, não nos resultados principais.
adm_max   <- function(f) max(simular_irf(p_bcb, "cambio", comb_adm = f)$adm_acum4)
FRACAO_ADM <- uniroot(function(f) adm_max(f) - BOXE_TEXTO[["cambio_adm"]], c(0, 0.3), tol = 1e-6)$root
cat(sprintf("FRACAO_ADM = %.3f (efeito máximo nos administrados, modas do BCB: %.2f)\n",
            FRACAO_ADM, adm_max(FRACAO_ADM)))

#Respostas no nosso modo, com bandas, e nas modas do BCB
CHOQUES <- c(juros = "juros_regra", cambio = "cambio", hiato = "hiato")
VARS    <- c("selic", "hiato", "cambio_acum", "ipca_acum4", "livres_acum4", "adm_acum4")
VERSOES <- list(`Nosso modo`                  = list(par = p_modo, f = 0),
                `Modas do BCB`                = list(par = p_bcb,  f = 0),
                `Nosso modo, adm. no câmbio`   = list(par = p_modo, f = FRACAO_ADM),
                `Modas do BCB, adm. no câmbio` = list(par = p_bcb,  f = FRACAO_ADM))

set.seed(13)
sorteios <- sample(nrow(amostras), N_DRAWS)
irf <- list(); picos_post <- list(); versoes <- list()
for (nm in names(CHOQUES)) {
  versoes[[nm]] <- lapply(VERSOES, function(v) simular_irf(v$par, CHOQUES[[nm]], comb_adm = v$f))
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
    r[[paste0(v, "_boxe")]] <- if (!is.null(BOXE[[nm]][[v]])) BOXE[[nm]][[v]] else NA_real_
  }
  irf[[nm]] <- r
  picos_post[[nm]] <- sapply(VARS, function(v) apply(arr[, , v], 1, pico))
  cat(sprintf("%-7s pronto (%d sorteios sem solução estável)\n", nm, sum(is.na(arr[, 1, 1]))))
}

# Efeitos máximos contra o boxe (trimestres na convenção do boxe)
linha <- function(nm, var, boxe = NA, tri_boxe = NA) {
  df <- irf[[nm]]; pp <- picos_post[[nm]][, var]
  data.frame(choque = nm, variavel = var,
             modo = pico(df[[var]]), tri_modo = tri_pico(df, var),
             post_mediana = median(pp, na.rm = TRUE),
             post_q05 = quantile(pp, .05, na.rm = TRUE), post_q95 = quantile(pp, .95, na.rm = TRUE),
             modas_bcb = pico(df[[paste0(var, "_bcb")]]), tri_bcb = tri_pico(df, paste0(var, "_bcb")),
             boxe = boxe, tri_boxe = tri_boxe, row.names = NULL)
}
max_boxe <- function(nm) { x <- BOXE[[nm]]$ipca_acum4; c(pico(x), which.max(abs(x)) - 1) }
tab_irf <- rbind(
  linha("juros",  "ipca_acum4",   max_boxe("juros")[1],  max_boxe("juros")[2]),
  linha("juros",  "livres_acum4"),
  linha("juros",  "hiato"),
  linha("cambio", "ipca_acum4",   max_boxe("cambio")[1], max_boxe("cambio")[2]),
  linha("cambio", "livres_acum4", BOXE_TEXTO[["cambio_livres"]]),
  linha("cambio", "adm_acum4",    BOXE_TEXTO[["cambio_adm"]]),
  linha("hiato",  "ipca_acum4",   max_boxe("hiato")[1],  max_boxe("hiato")[2]))
cat(sprintf("\n=== Efeitos máximos (p.p.), administrados indexados ===\n"))
imprimir(tab_irf)

# Trajetórias completas contra o anexo: efeito máximo, trimestre e REQM de cada versão
VAR_IPCA <- c(juros = "ipca_acum4", cambio = "ipca_acum4", hiato = "ipca_acum4")
tab_traj <- do.call(rbind, lapply(names(CHOQUES), function(nm) {
  bx <- BOXE[[nm]]$ipca_acum4
  rbind(data.frame(choque = nm, versao = "Boxe (anexo)", maximo = pico(bx), tri = which.max(abs(bx)) - 1, reqm = 0),
        do.call(rbind, lapply(names(VERSOES), function(v) {
          x <- versoes[[nm]][[v]]$ipca_acum4
          data.frame(choque = nm, versao = v, maximo = pico(x), tri = which.max(abs(x)) - 1, reqm = reqm(x, bx))
        })))
}))
cat("\n=== IPCA acumulado em quatro trimestres: trajetórias contra o anexo ===\n"); imprimir(tab_traj)
for (nm in names(CHOQUES)) with(irf[[nm]], cat(sprintf(
  "%-7s trajetória do boxe dentro da nossa faixa de 90%%: %d de 16 trimestres\n",
  nm, sum(ipca_acum4_boxe[-1] >= ipca_acum4_q05[-1] & ipca_acum4_boxe[-1] <= ipca_acum4_q95[-1]))))
cat(sprintf("Selic após o choque de juros, modas do BCB x anexo: REQM %.3f (trimestres 0 a 8: %.3f)\n",
            reqm(versoes$juros$`Modas do BCB`$selic, BOXE$juros$selic),
            reqm(versoes$juros$`Modas do BCB`$selic[1:9], BOXE$juros$selic[1:9])))
cat(sprintf("Hiato após o choque de demanda, modas do BCB x anexo: REQM %.3f\n",
            reqm(versoes$hiato$`Modas do BCB`$hiato, BOXE$hiato$hiato)))


# 3. Canais de transmissão: trajetória fixa da Selic
CANAIS <- list(`Expectativa e câmbio livres`    = c(cambio = TRUE,  exp = TRUE),
               `Expectativa livre, câmbio fixo` = c(cambio = FALSE, exp = TRUE),
               `Expectativa e câmbio fixos`     = c(cambio = FALSE, exp = FALSE))
canais <- data.frame(trimestre = 0:(H_IRF - 1))
tab_canais <- do.call(rbind, lapply(names(CANAIS), function(nm) {
  do.call(rbind, lapply(names(VERSOES), function(v) {
    r <- simular_irf(VERSOES[[v]]$par, "juros", canal_cambio = CANAIS[[nm]][["cambio"]],
                     canal_exp = CANAIS[[nm]][["exp"]], comb_adm = VERSOES[[v]]$f)
    if (v == "Nosso modo") canais[[nm]] <<- r$ipca_acum4
    data.frame(canal = nm, versao = v, maximo = pico(r$ipca_acum4), tri = tri_pico(r, "ipca_acum4"),
               tri4 = r$ipca_acum4[5], tri16 = r$ipca_acum4[17], reqm = reqm(r$ipca_acum4, BOXE_CANAIS[[nm]]))
  }))
}))
cat("\n=== Canais de transmissão (Selic em trajetória fixa) contra o anexo ===\n")
tab_canais <- rbind(data.frame(canal = names(CANAIS), versao = "Boxe (anexo)",
                               maximo = sapply(names(CANAIS), function(nm) pico(BOXE_CANAIS[[nm]])),
                               tri = sapply(names(CANAIS), function(nm) which.max(abs(BOXE_CANAIS[[nm]])) - 1),
                               tri4 = sapply(names(CANAIS), function(nm) BOXE_CANAIS[[nm]][5]),
                               tri16 = sapply(names(CANAIS), function(nm) BOXE_CANAIS[[nm]][17]), reqm = 0),
                    tab_canais)
imprimir(tab_canais[order(match(tab_canais$canal, names(CANAIS)), tab_canais$versao != "Boxe (anexo)"), ])

#De onde vêm as diferenças
# Troca de blocos: parte das modas do BCB e troca um bloco por vez pelo nosso,
# com os administrados indexados e com a fração FRACAO_ADM no câmbio
BLOCOS <- list(Taylor = c("theta_1", "theta_2", "theta_3"), IS = c("beta1", "beta2", "beta5"),
               Phillips = c("alpha_1L", "alpha_1I", "alpha_2", "alpha_3", "alpha_4"),
               Expectativas = c("phi_1", "phi_2", "phi_3"), Paridade = "delta")
resumo_juros <- function(par, rotulo) {
  r0 <- simular_irf(par, "juros_regra"); rf <- simular_irf(par, "juros_regra", comb_adm = FRACAO_ADM)
  data.frame(parametros = rotulo, ipca = pico(r0$ipca_acum4), tri = tri_pico(r0, "ipca_acum4"),
             ipca_adm_cambio = pico(rf$ipca_acum4), tri_adm_cambio = tri_pico(rf, "ipca_acum4"),
             hiato = pico(r0$hiato), selic_t6 = r0$selic[7])
}
blocos <- rbind(
  resumo_juros(p_bcb, "Modas do BCB"),
  do.call(rbind, lapply(names(BLOCOS), function(b) {
    pp <- p_bcb; pp[BLOCOS[[b]]] <- p_modo[BLOCOS[[b]]]
    resumo_juros(pp, paste("BCB com o nosso bloco:", b))
  })),
  resumo_juros(p_modo, "Nosso modo"))
cat(sprintf("\n=== Choque de juros: troca de blocos de parâmetros (boxe: %.2f no trimestre %d) ===\n",
            max_boxe("juros")[1], max_boxe("juros")[2])); imprimir(blocos)

# Fração dos administrados no câmbio: como cada valor aproxima o simulador do boxe
# (modas do BCB). O choque de juros não entra na calibração e serve de validação.
fracoes <- do.call(rbind, lapply(sort(c(0, 0.05, 0.10, FRACAO_ADM, 0.15, 0.20, 0.25)), function(f) {
  rc <- simular_irf(p_bcb, "cambio", comb_adm = f); rj <- simular_irf(p_bcb, "juros_regra", comb_adm = f)
  data.frame(fracao = f, cambio_ipca = pico(rc$ipca_acum4), cambio_livres = pico(rc$livres_acum4),
             cambio_adm = pico(rc$adm_acum4), reqm_cambio = reqm(rc$ipca_acum4, BOXE$cambio$ipca_acum4),
             juros_ipca = pico(rj$ipca_acum4), juros_tri = tri_pico(rj, "ipca_acum4"),
             reqm_juros = reqm(rj$ipca_acum4, BOXE$juros$ipca_acum4))
}))
cat(sprintf("\n=== Administrados no câmbio, modas do BCB (boxe: câmbio %.2f, livres %.2f, adm %.2f; juros %.2f) ===\n",
            max_boxe("cambio")[1], BOXE_TEXTO[["cambio_livres"]], BOXE_TEXTO[["cambio_adm"]], max_boxe("juros")[1]))
imprimir(fracoes)


#graficos
LARANJA <- "#C2571A"; VERDE <- "#3B8E4F"
EIXO_X <- scale_x_continuous(breaks = seq(0, 16, 2))
tema <- theme_minimal(base_size = TAM, base_family = FONTE) +
  theme(panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(),
        strip.text = element_text(hjust = 0, size = TAM), plot.margin = margin(4, 8, 4, 4))
linha_zero <- geom_hline(yintercept = 0, colour = CINZA, linewidth = 0.3)

SERIES <- c("Nosso modo", "Nosso modelo nas modas do BCB", "Boxe (anexo do RI)")
painel_ab <- function(df, va, rot_a) {
  rot_b  <- "B. IPCA acumulado em quatro trimestres"
  niveis <- c(rot_a, rot_b)
  lng <- function(v, painel) data.frame(trimestre = df$trimestre, painel = painel,
                                        lo = df[[paste0(v, "_q05")]], hi = df[[paste0(v, "_q95")]],
                                        modo = df[[v]], bcb = df[[paste0(v, "_bcb")]], boxe = df[[paste0(v, "_boxe")]])
  d <- rbind(lng(va, rot_a), lng("ipca_acum4", rot_b)); d$painel <- factor(d$painel, niveis)
  dl <- pivot_longer(d, c(modo, bcb, boxe), names_to = "serie", values_to = "valor")
  dl$serie <- factor(c(modo = SERIES[1], bcb = SERIES[2], boxe = SERIES[3])[dl$serie], SERIES)
  ggplot(d, aes(trimestre)) +
    linha_zero +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = "Intervalo de 90% (posteriori)")) +
    geom_line(data = dl, aes(y = valor, colour = serie, linetype = serie, linewidth = serie)) +
    facet_wrap(~ painel, scales = "free_y") +
    scale_fill_manual(values = adjustcolor(AZUL, 0.18)) +
    scale_colour_manual(values = c(AZUL, CINZA, "black")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_linewidth_manual(values = c(0.9, 0.6, 0.7)) +
    EIXO_X + scale_y_continuous(labels = virgula) +
    labs(x = "Trimestres", y = "p.p.") + tema +
    theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt")) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),
           linewidth = guide_legend(order = 1), fill = guide_legend(order = 2))
}
salvar(painel_ab(irf$juros, "selic", "A. Taxa Selic"), "fig5_1_juros.pdf", 9)
salvar(painel_ab(irf$hiato, "hiato", "A. Hiato do produto"), "fig5_3_hiato.pdf", 9)

#Câmbio. A: câmbio acumulado; B: IPCA (com banda), livres e administrados.
# Linha pontilhada azul: IPCA do boxe; linhas horizontais: máximos de livres e
# administrados informados no texto do boxe.
cb <- irf$cambio
NIVEIS <- c("IPCA", "Livres", "Administrados")
dB <- rbind(
  data.frame(trimestre = cb$trimestre, serie = "IPCA", valor = cb$ipca_acum4, lo = cb$ipca_acum4_q05, hi = cb$ipca_acum4_q95),
  data.frame(trimestre = cb$trimestre, serie = "Livres", valor = cb$livres_acum4, lo = NA, hi = NA),
  data.frame(trimestre = cb$trimestre, serie = "Administrados", valor = cb$adm_acum4, lo = NA, hi = NA))
dB$serie <- factor(dB$serie, NIVEIS)
refB <- data.frame(serie = factor(NIVEIS[2:3], NIVEIS), y = unname(BOXE_TEXTO[c("cambio_livres", "cambio_adm")]))
gA <- ggplot(data.frame(trimestre = cb$trimestre, valor = cb$cambio_acum), aes(trimestre, valor)) +
  linha_zero + geom_line(colour = "black", linewidth = 0.9) +
  EIXO_X + scale_y_continuous(labels = virgula, limits = c(0, 12)) +
  labs(x = "Trimestres", y = "%", subtitle = "A. Câmbio, variação acumulada") + tema
gB <- ggplot(dB, aes(trimestre, valor, colour = serie)) +
  linha_zero +
  geom_ribbon(data = dB[dB$serie == "IPCA", ], aes(ymin = lo, ymax = hi), fill = adjustcolor(AZUL, 0.18), colour = NA) +
  geom_hline(data = refB, aes(yintercept = y, colour = serie), linetype = "dotted", linewidth = 0.5, show.legend = FALSE) +
  geom_line(data = data.frame(trimestre = cb$trimestre, valor = cb$ipca_acum4_boxe, serie = factor("IPCA", NIVEIS)),
            linetype = "dotted", linewidth = 0.7, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(IPCA = AZUL, Livres = LARANJA, Administrados = VERDE), breaks = NIVEIS) +
  EIXO_X + scale_y_continuous(labels = virgula) +
  labs(x = "Trimestres", y = "p.p.", subtitle = "B. Inflação acumulada em quatro trimestres") + tema
salvar(patchwork::wrap_plots(gA, gB, widths = c(1, 1.25)) +
         patchwork::plot_layout(guides = "collect") & theme(legend.position = "top"),
       "fig5_2_cambio.pdf", 8.5)

#Canais de transmissão, no nosso modo; linhas pontilhadas: boxe (Gráfico 4B do anexo)
dc <- pivot_longer(canais, -trimestre, names_to = "canal", values_to = "valor")
dc$canal <- factor(dc$canal, names(CANAIS))
dcb <- pivot_longer(BOXE_CANAIS, -trimestre, names_to = "canal", values_to = "valor")
dcb$canal <- factor(dcb$canal, names(CANAIS))
g54 <- ggplot(dc, aes(trimestre, valor, colour = canal)) +
  linha_zero +
  geom_line(data = dcb, linetype = "dotted", linewidth = 0.6, show.legend = FALSE) +
  geom_line(aes(linetype = canal), linewidth = 0.9) +
  scale_colour_manual(values = c(AZUL, LARANJA, CINZA)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  EIXO_X + scale_y_continuous(labels = virgula) +
  labs(x = "Trimestres", y = "p.p.") + tema
salvar(g54, "fig5_4_canais.pdf", 8)