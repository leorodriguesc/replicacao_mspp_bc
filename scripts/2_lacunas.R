lp0 <- log_posterior(theta_modo)        # especificação principal, no modo

# PARTE A — De onde vêm as calibrações
# ---- A1 [resíduo] ----
# Critério: no modo, os erros de previsão padronizados de cada equação têm
# desvio-padrão próximo de 1 (entre 0,75 e 1,33) e média próxima de 0 (|média| < 0,25).
# Os desvios foram obtidos reescalando cada um pelo desvio-padrão do seu resíduo
# e reestimando o modo, até o critério ser atendido.
r  <- construir_kfas_lp(combinar_par(theta_modo), dados_exog, dados_obs, s_perfil)
kf <- KFS(r$modelo, filtering = "state", smoothing = "none")
Fm <- if (nrow(kf$F) == nrow(kf$v)) kf$F else t(kf$F)
z  <- kf$v / sqrt(Fm); colnames(z) <- nomes_obs
residuos <- t(sapply(c("pi_L", "selic", "dcambio", "pi_focus"), function(o) {
  x <- z[is.finite(z[, o]), o]
  c(media = mean(x), desvio = sd(x))
}))
cat("A1. Erros de previsão padronizados no modo:\n")
print(round(residuos, 2))
fora <- rownames(residuos)[residuos[, "desvio"] < 0.75 | residuos[, "desvio"] > 1.33 |
                             abs(residuos[, "media"]) > 0.25]
cat("Fora do critério:", if (length(fora)) paste(fora, collapse = ", ") else "nenhuma", "\n\n")

# ---- A2 [grade] ----
# A escala dos termos de crise entra na grade só como informação: é escolha própria
# (ver o comentário em ESPEC_PRINCIPAL, no 1_modelo.R).
# Critério: o valor de maior log-posteriori, reestimando o modo em cada ponto
# (inclusive no valor em uso, para que todos os pontos sejam tratados igual).
# Calibrar uma variância pela log-posteriori equivale a estimá-la por máxima
# verossimilhança numa grade.
GRADES <- list(sigma_eps_h  = c(0.10, 0.15, 0.25, 0.50, 1.00),
               sigma_pi_e   = c(0.15, 0.20, 0.25, 0.30, 0.40),
               escala_ph    = c(0.05, 0.10, 0.25, 0.50, 1.00),
               crise_escala = c(0.5, 1, 2))

if (file.exists("calibracao_grades.csv")) {
  grades <- read.csv("calibracao_grades.csv") 
} else {
  grades <- do.call(rbind, lapply(names(GRADES), function(nm) do.call(rbind, lapply(GRADES[[nm]], function(v) {
    em_uso <- isTRUE(all.equal(v, ESPEC_PRINCIPAL[[nm]]))
    o  <- com_espec(do.call(espec, setNames(list(v), nm)), otimizar(theta_modo))
    lp <- -o$value
    if (em_uso && lp > lp0 + 0.5) cat("  ATENÇÃO: o modo salvo não é o máximo da especificação principal\n")
    cat(sprintf("  %-12s = %5.2f | log-posteriori %9.2f\n", nm, v, lp))
    # β1 e β5 mostram como o choque da IS reparte a persistência do hiato
    data.frame(parametro = nm, valor = v, em_uso = em_uso, logpost = lp,
               beta1 = o$par[["beta1"]], beta5 = o$par[["beta5"]])
  }))))
  write.csv(grades, "calibracao_grades.csv", row.names = FALSE)
}
grades$dif_max <- ave(grades$logpost, grades$parametro, FUN = function(x) x - max(x))
cat("A2. Grades (dif_max = distância ao máximo da grade):\n")
imp <- grades; num <- sapply(imp, is.numeric); imp[num] <- round(imp[num], 2)
print(imp, row.names = FALSE)
for (nm in names(GRADES)) {
  g <- grades[grades$parametro == nm, ]
  cat(sprintf("  %-12s em uso %.2f | maior log-posteriori em %.2f%s\n", nm, ESPEC_PRINCIPAL[[nm]],
              g$valor[which.max(g$logpost)],
              if (which.max(g$logpost) %in% c(1, nrow(g))) " (na borda da grade: estender)" else ""))
}

# ---- A3 Critério do boxe para os choques das taxas neutras ----
# Critério do boxe: o desvio dos choques das taxas neutras é calibrado para que a
# taxa neutra tenha variância semelhante à da média e da mediana das medidas de
# taxa neutra do boxe "Atualização de medidas de taxa de juros real neutra no
# Brasil" (RI jun/2024, Gráfico 1).
# Os dados vêm do anexo estatístico (medianas de cinco grupos de medidas). Os grupos
# não cobrem o mesmo período: até 2005T1 só há o grupo Focus, e o NTN-B começa em
# 2007T1. Por isso o alvo é calculado de duas formas, nos trimestres com os cinco
# grupos e em toda a amostra (com os grupos disponíveis em cada trimestre). No modo.
a3 <- aba_anexo("C2 Boxe2 Graf 1"); a3 <- a3[grepl("^20[0-9]{2}:[1-4]$", a3[[1]]), ]
grupos <- data.frame(trimestre = as.Date(sprintf("%s-%02d-01", substr(a3[[1]], 1, 4),
                                                 3 * como_num(substr(a3[[1]], 6, 6)) - 2)),
                     sapply(a3[2:6], como_num))
names(grupos)[-1] <- c("focus", "ntnb", "externas", "af_bf", "modelos_bc")
datas <- dados_obs$trimestre[efetivo]
k     <- match(datas, grupos$trimestre)
vals  <- as.matrix(grupos[k, -1])
todos <- complete.cases(vals)
dp_alvo <- function(sel) c(media   = sd(rowMeans(vals[sel, , drop = FALSE], na.rm = TRUE)),
                           mediana = sd(apply(vals[sel, , drop = FALSE], 1, median, na.rm = TRUE)))
alvos <- rbind(cinco_grupos = dp_alvo(todos), amostra_toda = dp_alvo(TRUE))
cat("\nA3. Desvio-padrão das medidas do boxe (média e mediana dos grupos):\n"); print(round(alvos, 2))
cat(sprintf("Trimestres com os cinco grupos: %d (%s a %s) | tendência HP sozinha: %.2f nesses trimestres, %.2f na amostra toda\n",
            sum(todos), format(min(datas[todos])), format(max(datas[todos])),
            sd(dados_exog$juro_real_hp[efetivo][todos]), sd(dados_exog$juro_real_hp[efetivo])))
neutra <- t(sapply(c(0.05, 0.10, 0.15, 0.25, 0.40, 0.60, 1.00), function(sg) {
  rr <- com_espec(espec(sigma_rr = sg), {
    a <- KFS(construir_kfas_lp(combinar_par(theta_modo), dados_exog, dados_obs, s_perfil)$modelo,
             smoothing = "state")$alphahat
    dados_exog$juro_real_hp[efetivo] + a[, idx$r_is]
  })
  c(sigma_rr = sg, dp_cinco_grupos = sd(rr[todos]), dp_amostra_toda = sd(rr),
    neutra_2023T2 = rr[datas == TRI_NEUTRA])
}))
print(round(neutra, 2))
cat(sprintf(paste0("Alvo (média de média e mediana): %.2f nos cinco grupos e %.2f na amostra toda | ",
                   "menor desvio na grade: %.2f e %.2f | em uso: %.2f\n"),
            mean(alvos["cinco_grupos", ]), mean(alvos["amostra_toda", ]),
            min(neutra[, "dp_cinco_grupos"]), min(neutra[, "dp_amostra_toda"]), ESPEC_PRINCIPAL$sigma_rr))
cat(sprintf("Taxa neutra do BCB em 2023T2 (modelo agregado): %.1f\n\n", neutra_bcb))

# PARTE B — Variantes: uma alternativa por lacuna (ordem da tabela da seção 3.4)
VARIANTES <- list(
  "Administrados pela expectativa Focus"      = espec(adm_projecao = "focus"),
  "Termos de crise da IS × 2"                 = espec(crise_escala = 2),
  "Choques das taxas neutras = 0,25"          = espec(sigma_rr = 0.25),
  "Variância inicial das taxas neutras = 100" = espec(p1_neutro = 100),
  "Choque das expectativas = 0,15"            = espec(sigma_pi_e = 0.15),
  "Choque da IS = 0,50"                       = espec(sigma_eps_h = 0.50))

# Hiato e taxa neutra suavizados, na especificação ativa
trajetorias <- function(th) {
  r <- construir_kfas_lp(combinar_par(th), dados_exog, dados_obs, s_perfil)
  a <- KFS(r$modelo, smoothing = "state")$alphahat
  list(hiato  = as.numeric(a[, idx$h]),
       neutra = dados_exog$juro_real_hp[r$idx_amostra] + as.numeric(a[, idx$r_is]))
}

# O que se compara em cada modo (unidades do boxe). As três últimas colunas
# comparam com o BCB (anexo estatístico, lido no 1_modelo.R): hiato do Gráfico 2.2.6
# e taxa neutra do modelo agregado em 2023T2.
datas_ef <- dados_obs$trimestre[efetivo]
h_bcb    <- hiato_bcb$bcb[match(datas_ef, hiato_bcb$trimestre)]
resumo <- function(th, lp, ref = NULL) {
  tr <- trajetorias(th)
  tb <- th / escala
  data.frame(
    logpost    = lp,
    modos_IC   = sum(tb >= IC_BCB["inf", ] & tb <= IC_BCB["sup", ]),
    alpha_1I   = tb[["alpha_1I"]],
    phi_2      = tb[["phi_2"]],
    lambda     = 1 - tb[["theta_1"]] - tb[["theta_2"]],
    efeito_LP  = (tb[["beta2"]] / 4) / (1 - tb[["beta1"]]),
    cor_hiato  = if (is.null(ref)) 1 else cor(tr$hiato, ref$hiato),
    dif_neutra = if (is.null(ref)) 0 else mean(tr$neutra - ref$neutra),
    cor_hiato_bcb  = cor(tr$hiato, h_bcb),
    reqm_hiato_bcb = sqrt(mean((tr$hiato - h_bcb)^2)),
    neutra_2023T2  = tr$neutra[datas_ef == TRI_NEUTRA])
}

ref  <- trajetorias(theta_modo)
base <- resumo(theta_modo, lp0)

# Cada variante reestima o modo partindo do modo principal
linhas <- lapply(names(VARIANTES), function(nm) {
  t0 <- Sys.time()
  l <- com_espec(VARIANTES[[nm]], {
    o <- otimizar(theta_modo)
    resumo(o$par, -o$value, ref)
  })
  cat(sprintf("B. %-42s %.1f min\n", nm, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  l
})
tab <- cbind(variante = c("Especificação principal", names(VARIANTES)),
             rbind(base, do.call(rbind, linhas)))
tab$dif_logpost <- tab$logpost - lp0

# PARTE C — Regra de decisão
# Uma escolha pesa se mudar algum resultado:
#   parâmetros: 2 ou mais modos mudam de lado do intervalo do BCB, ou α1I ou φ2
#               andam metade da distância até a moda do BCB
#   dinâmica:   λ ou o efeito de longo prazo andam metade da distância até o BCB
#   hiato e taxa neutra: correlação do hiato < 0,95 ou taxa neutra 0,5 p.p.
#               acima ou abaixo, em média
# Os dados decidem quando a log-posteriori muda 2 pontos ou mais. Exceção: a
# variância inicial muda a verossimilhança das primeiras observações (uma
# variância maior as torna menos previsíveis e reduz a verossimilhança), então a
# log-posteriori não serve para escolhê-la.
alvo_bcb <- c(alpha_1I = 0.38, phi_2 = 0.11, lambda = 0.10, efeito_LP = 0.733)
andou_metade <- function(v) abs(tab[[v]] - base[[v]]) >= abs(base[[v]] - alvo_bcb[[v]]) / 2

tab$pesa_parametros <- abs(tab$modos_IC - base$modos_IC) >= 2 |
  andou_metade("alpha_1I") | andou_metade("phi_2")
tab$pesa_dinamica   <- andou_metade("lambda") | andou_metade("efeito_LP")
tab$pesa_nao_obs    <- tab$cor_hiato < 0.95 | abs(tab$dif_neutra) >= 0.5
tab$dados_decidem   <- abs(tab$dif_logpost) >= 2
tab$dados_decidem[tab$variante == "Variância inicial das taxas neutras = 100"] <- NA

num <- sapply(tab, is.numeric); tab[num] <- round(tab[num], 3)
print(tab[, c("variante", "dif_logpost", "modos_IC", "alpha_1I", "phi_2", "lambda", "efeito_LP",
              "cor_hiato", "dif_neutra", "pesa_parametros", "pesa_dinamica", "pesa_nao_obs",
              "dados_decidem")], row.names = FALSE)

# Comparação com o BCB: alguma variante aproxima o hiato ou a taxa neutra dos do BCB?
# (informativa: não entra na regra de decisão)
cat(sprintf("\nComparação com o BCB (hiato do Gráfico 2.2.6; taxa neutra do modelo agregado em 2023T2 = %.1f):\n",
            neutra_bcb))
print(tab[, c("variante", "cor_hiato_bcb", "reqm_hiato_bcb", "neutra_2023T2")], row.names = FALSE)