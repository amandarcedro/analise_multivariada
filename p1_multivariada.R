# =====================================================================
# ME731 - Projeto 1: MANOVA - aspectos metodologicos e aplicacao
# Dados: Seeds (UCI, id 236) - 210 graos, 3 variedades, 7 atributos
# Pacotes: apenas base R + MASS (MASS ja vem com a instalacao do R)
# Reprodutibilidade: gerador fixado (Mersenne-Twister), set.seed(731) e
# numero FIXO de repeticoes no Monte Carlo => mesmos resultados a cada
# execucao (na mesma versao do R; diferencas minimas de arredondamento
# podem ocorrer entre plataformas).
# Saidas graficas (PDF) na pasta ./figuras. Arquivo 100% ASCII; acentos dos
# graficos via escapes \u.
# =====================================================================
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(731)
fig_dir <- if (dir.create("figuras", showWarnings = FALSE) || dir.exists("figuras")) "figuras" else tempdir()
suppressPackageStartupMessages(library(MASS))

# Formata p-valores: nunca imprime potencias absurdas (ex.: 1e-81)
fp <- function(x) ifelse(x < 2.2e-16, "< 2.2e-16", formatC(x, format = "e", digits = 2))

# ---------------------------------------------------------------------
# 1. DADOS
# ---------------------------------------------------------------------
nomes <- c("area", "perimetro", "compacidade", "comprimento",
           "largura", "assimetria", "sulco", "variedade")

# Fontes em ordem de preferencia (a segunda e o espelho original da UCI)
urls <- "https://raw.githubusercontent.com/Amandacedro/Multivariada_P1/main/seeds_dataset.txt"
dados <- read.table(urls, header = FALSE)
names(dados) <- nomes
if (!is.factor(dados$variedade))
  dados$variedade <- factor(dados$variedade, levels = 1:3,
                            labels = c("Kama", "Rosa", "Canadian"))
dados <- dados[complete.cases(dados), ]
# Checagem de integridade (media conhecida de area no banco original)
stopifnot(nrow(dados) == 210, abs(mean(dados$area) - 14.847) < 0.01,
          all(table(dados$variedade) == 70))
grp  <- droplevels(dados$variedade)
lv   <- levels(grp)
vars <- nomes[1:7]
rot  <- c("\u00e1rea", "per\u00edmetro", "compacidade", "comprimento",
          "largura", "assimetria", "sulco")   # rotulos dos graficos
Y <- as.matrix(dados[, vars])
n <- nrow(Y); p <- ncol(Y); g <- nlevels(grp); s <- min(g - 1, p)
ni <- as.numeric(table(grp))
stopifnot(n >= 30, p >= 4, n > p, is.numeric(Y))
cat("n =", n, " p =", p, " g =", g, "\n"); print(table(grp))

# ---------------------------------------------------------------------
# 2. FUNCOES AUXILIARES (formulas do relatorio)
# ---------------------------------------------------------------------
# Matrizes SSCP: T = H + E (total = entre + dentro)
sscp <- function(Y, grp) {
  E <- Reduce(`+`, lapply(split(as.data.frame(Y), grp), function(x)
    crossprod(scale(as.matrix(x), scale = FALSE))))
  Tt <- crossprod(scale(Y, scale = FALSE))
  list(E = E, H = Tt - E, T = Tt)
}
autovalores <- function(sc, s) {
  Re(eigen(solve(sc$E, sc$H), only.values = TRUE)$values)[1:s]
}
# As 4 estatisticas a partir dos autovalores de E^{-1}H
est4 <- function(Yd, gp) {
  nvec <- as.numeric(table(gp))
  Yc <- scale(Yd, scale = FALSE)
  m  <- rowsum(Yc, gp)                       # somas por grupo (dados centrados)
  H  <- crossprod(m / sqrt(nvec))
  E  <- crossprod(Yc) - H
  lam <- Re(eigen(solve(E, H), only.values = TRUE)$values)[1:s]
  c(Wilks = prod(1 / (1 + lam)), Pillai = sum(lam / (1 + lam)),
    `Hotelling-Lawley` = sum(lam), Roy = lam[1])
}
wilks_det <- function(Yd, gp) { sb <- sscp(Yd, gp); det(sb$E) / det(sb$T) }

mardia <- function(X) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- crossprod(Xc) / n
  D  <- Xc %*% solve(S) %*% t(Xc)
  b1 <- sum(D^3) / n^2
  b2 <- mean(diag(D)^2)
  skew <- n * b1 / 6; gl <- p * (p + 1) * (p + 2) / 6
  kurt <- (b2 - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  c(assim = skew, p_assim = pchisq(skew, gl, lower.tail = FALSE),
    curt = kurt, p_curt = 2 * pnorm(abs(kurt), lower.tail = FALSE))
}

box_m <- function(Y, grp) {
  g <- nlevels(grp); p <- ncol(Y); N <- nrow(Y); ni <- as.numeric(table(grp))
  Si <- lapply(split(as.data.frame(Y), grp), cov)
  Sp <- Reduce(`+`, Map(function(S, k) (k - 1) * S, Si, ni)) / (N - g)
  M  <- (N - g) * log(det(Sp)) -
    sum(mapply(function(S, k) (k - 1) * log(det(S)), Si, ni))
  c1 <- (sum(1 / (ni - 1)) - 1 / (N - g)) *
    (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (g - 1))
  est <- M * (1 - c1); gl <- (g - 1) * p * (p + 1) / 2
  c(M = M, qui2 = est, gl = gl, p_valor = pchisq(est, gl, lower.tail = FALSE))
}

# ---------------------------------------------------------------------
# 3. DESCRITIVA E MULTICOLINEARIDADE
# ---------------------------------------------------------------------
cat("\nMedias por variedade:\n")
print(round(do.call(rbind, lapply(split(as.data.frame(Y), grp), colMeans)), 3))
cat("\nDesvios-padrao por variedade:\n")
print(round(do.call(rbind, lapply(split(as.data.frame(Y), grp),
                                  function(x) apply(x, 2, sd))), 3))
cat("\nCoeficiente de variacao (%) da area por variedade:\n")
print(round(100 * tapply(Y[, "area"], grp, sd) / tapply(Y[, "area"], grp, mean), 1))
R_geral <- cor(Y)
cat("\nCorrelacoes (total):\n"); print(round(R_geral, 3))
cat("Numero de condicao da matriz de correlacao:",
    round(kappa(R_geral, exact = TRUE), 1), "\n")
cat("Autovalores da matriz de correlacao:", round(eigen(R_geral)$values, 4), "\n")
cat("VIF (diagonal da inversa da correlacao):\n")
print(round(setNames(diag(solve(R_geral)), vars), 1))

# Boxplots: 2 x 4 (o 8o painel traz a legenda de cores)
pal <- c("#5568af", "#ec1763", "#cdd629")
pdf(file.path(fig_dir, "boxplots.pdf"), width = 11, height = 6)
par(mfrow = c(2, 4), mar = c(3, 3, 2.5, 1))
for (j in seq_along(vars))
  boxplot(dados[[vars[j]]] ~ grp, main = rot[j], col = adjustcolor(pal, 0.35),
          xlab = "", ylab = "", cex.main = 1.2, cex.axis = 1.1)
plot.new()
legend("center", lv, fill = adjustcolor(pal, 0.35), border = pal, bty = "n",
       cex = 1.4, title = "Variedade")
invisible(dev.off())

pdf(file.path(fig_dir, "dispersao.pdf"), width = 8, height = 8)
pairs(Y, labels = rot, col = pal[grp], pch = 19, cex = 0.5,
      main = "Dispers\u00e3o por variedade (azul = Kama, rosa = Rosa, verde = Canadian)",
      cex.main = 0.9)
invisible(dev.off())

# ---------------------------------------------------------------------
# 4. PRESSUPOSTOS
# ---------------------------------------------------------------------
fit <- manova(Y ~ grp)
res <- residuals(fit)

diagnosticos <- function(Y, grp, rotulo) {
  cat("\n---- Diagnosticos:", rotulo, "----\n")
  cat("Mardia por grupo (assim, p, curt, p):\n")
  print(round(t(sapply(split(as.data.frame(Y), grp), mardia)), 4))
  cat("M de Box:\n"); print(round(box_m(Y, grp), 4))
}
diagnosticos(Y, grp, "escala original")
diagnosticos(log(Y), grp, "escala log (sensibilidade)")

cat("\nShapiro-Wilk dos residuos (p-valores por variavel; Bonferroni = 0,05/7):\n")
print(round(sapply(1:p, function(j) shapiro.test(res[, j])$p.value), 4))

# QQ-plot classico (residuos agrupados) e MCD robusto DENTRO de cada grupo
d2 <- mahalanobis(res, rep(0, p), cov(res))
corte999 <- qchisq(0.999, p)
im <- which.max(d2)
cat("\nD2 classico (residuos agrupados): maior =", round(d2[im], 1),
    "na linha", im, "(", as.character(grp[im]), ") | corte qui2 0,999 =", round(corte999, 1),
    "| linhas acima do corte:", paste(which(d2 > corte999), collapse = ", "), "\n")
q <- qchisq(ppoints(n), p)
pdf(file.path(fig_dir, "mahalanobis.pdf"), width = 6.5, height = 5.5)
plot(q, sort(d2), xlab = "Quantis do qui-quadrado",
     ylab = "D\u00b2 cl\u00e1ssico (res\u00edduos agrupados)",
     main = "QQ-plot das dist\u00e2ncias de Mahalanobis")
abline(0, 1, col = "red"); abline(h = corte999, lty = 3, col = "blue")
text(max(q), max(d2), paste0("linha ", im, " (", as.character(grp[im]), ")"), pos = 2, cex = 0.8)
legend("topleft", c("y = x", "corte de 99,9%"), col = c("red", "blue"), lty = c(1, 3), cex = 0.8)
invisible(dev.off())

# MCD com cobertura de 75% (h = 53 de 70): subconjuntos de ~56% (padrao)
# ficam degenerados com atributos quase colineares (r ~ 0,95-0,99).
pdf(file.path(fig_dir, "mahalanobis_robusto.pdf"), width = 12, height = 4.5)
par(mfrow = c(1, 3))
cat("\nOutliers robustos por grupo (MCD com h = 75%; limiares qui2 0,975 e 0,999):\n")
for (k in 1:g) {
  idx <- which(grp == lv[k]); Xg <- Y[idx, ]
  h <- ceiling(0.75 * length(idx))
  mcdk <- tryCatch(cov.rob(Xg, method = "mcd", quantile.used = h, nsamp = "best"),
                   error = function(e) NULL)
  if (is.null(mcdk)) {
    cat(lv[k], ": MCD falhou (covariancia singular); usando estimativas classicas.\n")
    mcdk <- list(center = colMeans(Xg), cov = cov(Xg))
  }
  d2k  <- mahalanobis(Xg, mcdk$center, mcdk$cov)
  top  <- order(d2k, decreasing = TRUE)[1:5]
  d2c  <- mahalanobis(Xg, colMeans(Xg), cov(Xg))
  cat(lv[k], ": cond. da correlacao classica =", round(kappa(cov2cor(cov(Xg)), exact = TRUE)),
      "| robusta =", round(kappa(cov2cor(mcdk$cov), exact = TRUE)),
      "| D2 classico > 0,975:", sum(d2c > qchisq(0.975, p)), "\n  ")
  cat(lv[k], ": n_g =", length(idx), "| h =", h,
      "| > 0,975:", sum(d2k > qchisq(0.975, p)),
      "| > 0,999:", sum(d2k > qchisq(0.999, p)),
      "| 5 maiores (linha do arquivo, D2):",
      paste0(idx[top], " (", round(d2k[top], 1), ")", collapse = ", "), "\n")
  plot(qchisq(ppoints(length(idx)), p), sort(d2k), log = "xy",
       xlab = "Quantis do qui-quadrado (escala log)", ylab = "D\u00b2 robusto (MCD, escala log)",
       main = paste("QQ-plot robusto (MCD) -", lv[k])); abline(0, 1, col = "red")
}
invisible(dev.off())

# ---------------------------------------------------------------------
# 5. MANOVA - estatisticas a partir dos autovalores de E^{-1}H
# ---------------------------------------------------------------------
cat("\n===== MANOVA =====\n")
for (tst in c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"))
  print(summary(fit, test = tst))
cat("\nANOVAs univariadas:\n"); print(summary.aov(fit))

sc  <- sscp(Y, grp)
lam <- autovalores(sc, s)
manual <- c(Wilks = prod(1 / (1 + lam)), Pillai = sum(lam / (1 + lam)),
            `Hotelling-Lawley` = sum(lam), Roy = lam[1])
usado <- sapply(names(manual), function(t) summary(fit, test = t)$stats[1, t])
cat("\nAutovalores de E^-1 H:", round(lam, 4), "\n")
cat("Correlacoes canonicas ao quadrado:", round(lam / (1 + lam), 4), "\n")
cat("Conferencia manual vs summary():\n")
print(round(rbind(manual = manual, summary = usado), 5))
if (!isTRUE(all.equal(unname(manual), unname(usado), tolerance = 1e-6)))
  warning("Estatisticas manuais diferem de summary(): revisar.")
if (!isTRUE(all.equal(unname(est4(Y, grp)), unname(manual), tolerance = 1e-8)))
  warning("est4() difere do calculo manual: revisar.")

Lambda_full <- unname(manual["Wilks"])
cat("\nTamanho de efeito: eta2 multivariado = 1 - Lambda^(1/s) =",
    round(1 - Lambda_full^(1 / s), 4), "\n")
cat("Tamanho de efeito (Pillai): V/s =", round(unname(manual["Pillai"]) / s, 4), "\n")
cat("Proporcao da separacao por eixo canonico:",
    round(lam / sum(lam), 4), "\n")

# 5.0 Teste de dimensionalidade (Bartlett): H0: apenas k-1 dimensoes
cat("\nTeste de dimensionalidade (Bartlett): H0_k = so k-1 eixos significativos\n")
dimtab <- t(sapply(1:s, function(k) {
  Lk <- prod(1 / (1 + lam[k:s]))
  X2 <- -(n - 1 - (p + g) / 2) * log(Lk)
  gl <- (p - k + 1) * (g - k)
  c(k = k, qui2 = X2, gl = gl, p = pchisq(X2, gl, lower.tail = FALSE))
}))
print(data.frame(k = dimtab[, 1], qui2 = round(dimtab[, 2], 1), gl = dimtab[, 3],
                 p = fp(dimtab[, 4])))

# 5.1 Sensibilidade: MANOVA na escala log
Ylog <- log(Y)
fit_log <- manova(Ylog ~ grp)
est_log <- est4(Ylog, grp)
resumo <- function(e) c(Wilks = unname(e["Wilks"]), Pillai = unname(e["Pillai"]),
                        eta2 = unname(1 - e["Wilks"]^(1 / s)))
cat("\nMANOVA na escala log (sensibilidade):\n")
print(summary(fit_log, test = "Wilks")); print(summary(fit_log, test = "Pillai"))
cat("Comparacao original vs log:\n")
print(round(rbind(original = resumo(manual), log = resumo(est_log)), 4))

# 5.2 Lambda parcial: informacao adicional de cada variavel dada as demais
parcial <- function(Yx, Lfull) {
  tb <- t(sapply(seq_len(p), function(j) {
    Lj <- Lfull / wilks_det(Yx[, -j, drop = FALSE], grp)
    Fj <- (1 - Lj) / Lj * (n - g - p + 1) / (g - 1)
    pj <- pf(Fj, g - 1, n - g - p + 1, lower.tail = FALSE)
    c(Lambda_parcial = Lj, F = Fj, p = pj, p_bonf = min(1, pj * p))
  }))
  rownames(tb) <- vars
  tb[order(tb[, "Lambda_parcial"]), ]
}
fmt_parc <- function(tb) data.frame(Lambda_parcial = round(tb[, 1], 4), F = round(tb[, 2], 2),
                                    p = fp(tb[, 3]), p_bonf = fp(tb[, 4]),
                                    row.names = rownames(tb))
cat("\nLambda parcial por variavel, escala ORIGINAL (menor = mais informativa; gl =",
    g - 1, "e", n - g - p + 1, "):\n")
print(fmt_parc(parcial(Y, Lambda_full)))
cat("\nLambda parcial por variavel, escala LOG (compacidade = 4*pi*A/P^2 vira combinacao linear):\n")
print(fmt_parc(parcial(Ylog, unname(est_log["Wilks"]))))

# 5.3 Sensibilidade: sem os tres graos de Kama (linhas 60-62) mais discrepantes
out_idx <- c(60, 61, 62)
keep <- setdiff(seq_len(n), out_idx)
Y_s <- Y[keep, ]; g_s <- droplevels(grp[keep])
diagnosticos(Y_s, g_s, "sem as linhas 60-62 (Kama)")
fit_s <- manova(Y_s ~ g_s)
cat("\nMANOVA sem as linhas 60-62:\n")
print(summary(fit_s, test = "Wilks")); print(summary(fit_s, test = "Pillai"))
cat("eta2 multivariado sem as linhas 60-62:",
    round(1 - unname(est4(Y_s, g_s)["Wilks"])^(1 / s), 4), "\n")

# ---------------------------------------------------------------------
# 6. CONTRASTES (T2 de Hotelling, pooled), IC SIMULTANEOS, MAHALANOBIS
# ---------------------------------------------------------------------
ybar <- do.call(rbind, lapply(split(as.data.frame(Y), grp), colMeans))
Sp   <- sc$E / (n - g)
gl2  <- n - g - p + 1
contraste <- function(cvec, alfa = 0.05) {
  L  <- drop(t(cvec) %*% ybar)
  f  <- sum(cvec^2 / ni)
  T2 <- drop(t(L) %*% solve(Sp) %*% L) / f
  Fs <- gl2 / (p * (n - g)) * T2
  pv <- pf(Fs, p, gl2, lower.tail = FALSE)
  crit <- (n - g) * p / gl2 * qf(1 - alfa, p, gl2)
  ic <- cbind(estim = L, inf = L - sqrt(crit * diag(Sp) * f),
              sup = L + sqrt(crit * diag(Sp) * f))
  list(T2 = T2, F = Fs, p = pv, ic = ic)
}
C <- list(`Kama - Rosa`      = c(1, -1, 0),
          `Kama - Canadian`  = c(1, 0, -1),
          `Rosa - Canadian`  = c(0, 1, -1),
          `Canadian - media(Kama,Rosa)` = c(-.5, -.5, 1))
stopifnot(identical(lv, c("Kama", "Rosa", "Canadian")))
tab <- t(sapply(C, function(cv) { r <- contraste(cv); c(T2 = r$T2, F = r$F, p = r$p) }))
tab <- cbind(tab, p_holm = p.adjust(tab[, "p"], "holm"))
cat("\nContrastes (p ajustado por Holm):\n")
print(data.frame(T2 = round(tab[, "T2"], 1), F = round(tab[, "F"], 1),
                 p = fp(tab[, "p"]), p_holm = fp(tab[, "p_holm"]), row.names = rownames(tab)))

# IC simultaneos (T2) nas 7 variaveis; Bonferroni nos 3 pares => nivel 0,05/3 por par
cat("\nIC simultaneos (T2) para os 3 pares; cobertura conjunta >= 95% (alfa = 0,05/3 por par):\n")
for (nm in names(C)[1:3]) {
  cat("\n", nm, ":\n"); print(round(contraste(C[[nm]], alfa = 0.05 / 3)$ic, 3))
}

# Distancia de Mahalanobis entre medias: T2 = n_i n_j / (n_i + n_j) * D2
D2m <- matrix(0, g, g, dimnames = list(lv, lv))
for (i in 1:g) for (j in 1:g) {
  d <- ybar[i, ] - ybar[j, ]
  D2m[i, j] <- drop(t(d) %*% solve(Sp) %*% d)
}
cat("\nDistancias de Mahalanobis (D2) entre medias das variedades:\n")
print(round(D2m, 3))
T2_via_D2 <- c(`Kama - Rosa` = ni[1] * ni[2] / (ni[1] + ni[2]) * D2m[1, 2],
               `Kama - Canadian` = ni[1] * ni[3] / (ni[1] + ni[3]) * D2m[1, 3],
               `Rosa - Canadian` = ni[2] * ni[3] / (ni[2] + ni[3]) * D2m[2, 3])
cat("T2 = n_i n_j/(n_i+n_j) * D2 :", round(T2_via_D2, 2), "\n")
if (!isTRUE(all.equal(unname(T2_via_D2), unname(tab[1:3, "T2"]), tolerance = 1e-8)))
  warning("T2 via D2 difere do T2 dos contrastes: revisar.")

# Erro de Bayes (normal, covariancias iguais, priors iguais): Phi(-D/2)
bayes <- c(`Kama - Rosa` = pnorm(-sqrt(D2m[1, 2]) / 2),
           `Kama - Canadian` = pnorm(-sqrt(D2m[1, 3]) / 2),
           `Rosa - Canadian` = pnorm(-sqrt(D2m[2, 3]) / 2))
cat("Erro de Bayes teorico por par (Phi(-D/2), %):", round(100 * bayes, 2), "\n")

# ---------------------------------------------------------------------
# 7. PERMUTACAO (4 estatisticas), PERMUTACAO POR PAR E BOOTSTRAP
#    Obs.: a permutacao pressupoe intercambiabilidade sob H0 (covariancias
#    iguais); com o M de Box rejeitando, e uma aproximacao (ver relatorio).
# ---------------------------------------------------------------------
B <- 2000
obs4  <- est4(Y, grp)
perm4 <- replicate(B, est4(Y, sample(grp)))
tol <- 1e-10
p_perm <- c(Wilks = (1 + sum(perm4["Wilks", ] <= obs4["Wilks"] + tol)) / (B + 1),
            Pillai = (1 + sum(perm4["Pillai", ] >= obs4["Pillai"] - tol)) / (B + 1),
            `Hotelling-Lawley` = (1 + sum(perm4["Hotelling-Lawley", ] >= obs4["Hotelling-Lawley"] - tol)) / (B + 1),
            Roy = (1 + sum(perm4["Roy", ] >= obs4["Roy"] - tol)) / (B + 1))
cat("\nTeste de permutacao (B =", B, "); p minimo atingivel = 1/(B+1) =",
    round(1 / (B + 1), 5), "\n")
print(rbind(observado = round(obs4, 4), p_permutacao = round(p_perm, 5)))

# Permutacao por par: T2 de Hotelling de duas amostras (covariancia agrupada do par)
t2_dois <- function(X1, X2) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  d  <- colMeans(X1) - colMeans(X2)
  S  <- ((n1 - 1) * cov(X1) + (n2 - 1) * cov(X2)) / (n1 + n2 - 2)
  n1 * n2 / (n1 + n2) * drop(t(d) %*% solve(S, d))
}
pares <- list(`Kama - Rosa` = c(1, 2), `Kama - Canadian` = c(1, 3), `Rosa - Canadian` = c(2, 3))
perm_par <- t(sapply(pares, function(pr) {
  sel <- grp %in% lv[pr]; X <- Y[sel, ]; gp <- droplevels(grp[sel])
  obs <- t2_dois(X[gp == lv[pr[1]], ], X[gp == lv[pr[2]], ])
  pm  <- replicate(B, { gp2 <- sample(gp)
  t2_dois(X[gp2 == lv[pr[1]], ], X[gp2 == lv[pr[2]], ]) })
  c(T2 = obs, p_perm = (1 + sum(pm >= obs - tol)) / (B + 1))
}))
perm_par <- cbind(perm_par, p_holm = p.adjust(perm_par[, "p_perm"], "holm"))
cat("\nPermutacao por par (T2 de duas amostras, B =", B, "):\n")
print(round(perm_par, 5))

eta2_boot <- replicate(B, {
  idx <- unlist(lapply(split(seq_len(n), grp), function(i) sample(i, replace = TRUE)))
  sb  <- sscp(Y[idx, ], grp[idx])
  1 - (det(sb$E) / det(sb$T))^(1 / s)
})
cat("IC bootstrap 95% (percentil) para eta2 multivariado:",
    round(quantile(eta2_boot, c(.025, .975)), 4), "\n")

# ---------------------------------------------------------------------
# 8. VARIAVEIS CANONICAS
# ---------------------------------------------------------------------
ev <- eigen(solve(sc$E, sc$H))
V  <- Re(ev$vectors[, 1:2])
rownames(V) <- vars
V  <- sweep(V, 2, sqrt(diag(t(V) %*% Sp %*% V)), "/")   # var intra = 1
colnames(V) <- c("CV1", "CV2")
escores <- scale(Y, scale = FALSE) %*% V
# Ancoragem de sinal (autovetores tem sinal arbitrario): CV1 correlaciona
# positivamente com area; CV2 correlaciona negativamente com assimetria.
sg <- c(sign(cor(Y[, "area"], escores[, 1])), -sign(cor(Y[, "assimetria"], escores[, 2])))
sg[sg == 0] <- 1
V <- sweep(V, 2, sg, "*"); escores <- sweep(escores, 2, sg, "*")
colnames(escores) <- c("CV1", "CV2")
cat("\nCoeficientes de estrutura (cor variavel x eixo canonico):\n")
print(round(cor(Y, escores), 3))
cat("\nCoeficientes canonicos padronizados:\n")
print(round(V * sqrt(diag(Sp)), 3))

elipse <- function(centro, S, nivel = 0.95, k = 200) {
  th <- seq(0, 2 * pi, length.out = k)
  z  <- cbind(cos(th), sin(th)) %*% chol(S) * sqrt(qchisq(nivel, 2))
  sweep(z, 2, centro, "+")
}
pct <- round(100 * lam / sum(lam))
elips <- lapply(1:g, function(k) {
  X <- escores[grp == lv[k], ]; elipse(colMeans(X), cov(X)) })
lim_x <- range(escores[, 1], sapply(elips, function(e) range(e[, 1])))
lim_y <- range(escores[, 2], sapply(elips, function(e) range(e[, 2])))
pdf(file.path(fig_dir, "canonicas.pdf"), width = 7, height = 6)
plot(escores, col = pal[grp], pch = 19, cex = .6, xlim = lim_x, ylim = lim_y,
     xlab = paste0("Vari\u00e1vel can\u00f4nica 1 (", pct[1], "% da separa\u00e7\u00e3o)"),
     ylab = paste0("Vari\u00e1vel can\u00f4nica 2 (", pct[2], "% da separa\u00e7\u00e3o)"),
     main = "Separa\u00e7\u00e3o das variedades")
for (k in 1:g) {
  X <- escores[grp == lv[k], ]
  lines(elips[[k]], col = pal[k])
  points(rbind(colMeans(X)), pch = 4, cex = 2, lwd = 2, col = pal[k])
}
legend("topright", lv, col = pal, pch = 19); invisible(dev.off())

# ---------------------------------------------------------------------
# 9. CLASSIFICACAO (ponte com discriminante) - validacao leave-one-out
# ---------------------------------------------------------------------
dd <- data.frame(Y, variedade = grp)
lda_cv <- lda(variedade ~ ., data = dd, CV = TRUE)
qda_cv <- qda(variedade ~ ., data = dd, CV = TRUE)
lda_lg <- lda(variedade ~ ., data = data.frame(Ylog, variedade = grp), CV = TRUE)
cm <- table(grp, lda_cv$class)
cat("\nLDA (LOO) - matriz de confusao:\n"); print(cm)
cat("Acuracia LDA:", round(mean(lda_cv$class == grp), 4),
    "| QDA:", round(mean(qda_cv$class == grp), 4),
    "| LDA na escala log:", round(mean(lda_lg$class == grp), 4), "\n")
ic_bin <- function(cl) {
  k <- sum(cl == grp)
  c(acertos = k, erros = n - k, acuracia = k / n, binom.test(k, n)$conf.int)
}
tab_acc <- rbind(LDA = ic_bin(lda_cv$class), QDA = ic_bin(qda_cv$class),
                 `LDA (log)` = ic_bin(lda_lg$class))
colnames(tab_acc)[4:5] <- c("IC95_inf", "IC95_sup")
cat("Acuracia LOO com IC binomial exato de 95%:\n"); print(round(tab_acc, 4))

# Erros LOO por par vs erro de Bayes teorico (contagem de erros entre as duas classes)
err_par <- c(`Kama - Rosa` = cm["Kama", "Rosa"] + cm["Rosa", "Kama"],
             `Kama - Canadian` = cm["Kama", "Canadian"] + cm["Canadian", "Kama"],
             `Rosa - Canadian` = cm["Rosa", "Canadian"] + cm["Canadian", "Rosa"])
cat("\nErros LOO por par (contagem):", err_par, "\n")
cat("Erro de Bayes teorico por par (%):", round(100 * bayes, 2), "\n")

# ---------------------------------------------------------------------
# 10. MONTE CARLO: erro tipo I e curva de poder das 4 estatisticas
#     R fixo em 1000 repeticoes por celula (resultados reproduziveis)
#     Roy: alem da aprox. F, poder com TAMANHO CORRIGIDO (corte = quantil
#     95% da nula simulada em cada cenario).
# ---------------------------------------------------------------------
R <- 1000                    # repeticoes por celula (FIXO)
deltas <- c(`H0` = 0, `d=0.05` = 0.05, `d=0.075` = 0.075, `d=0.10` = 0.10)
grande <- colMeans(Y)
mus_de <- function(dl) sweep(sweep(ybar, 2, grande) * dl, 2, grande, "+")
S_g  <- lapply(split(as.data.frame(Y), grp), cov)
ord  <- order(sapply(S_g, det), decreasing = TRUE)
ns_het <- numeric(g); ns_het[ord] <- c(30, 70, 110)  # menor n com maior covariancia

gera <- function(ns, mus, Sigs, gl_t = Inf) {
  Z <- lapply(seq_along(ns), function(k) {
    X <- mvrnorm(ns[k], rep(0, p), Sigs[[k]])
    if (is.finite(gl_t)) {
      w <- sqrt(rchisq(ns[k], gl_t) / gl_t)
      X <- (X / w) * sqrt((gl_t - 2) / gl_t)     # t multivariada, mesma cov
    }
    sweep(X, 2, mus[k, ], "+")
  })
  list(Y = do.call(rbind, Z), g = factor(rep(seq_along(ns), ns)))
}
# Retorna 4 indicadores de rejeicao (aprox. F) + a estatistica de Roy (lambda_1)
rejeita <- function(d, alfa = 0.05) {
  f <- manova(d$Y ~ d$g)
  tst <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
  st <- lapply(tst, function(t) summary(f, test = t)$stats[1, ])
  c(setNames(sapply(st, function(x) as.numeric(x[["Pr(>F)"]] < alfa)), tst),
    RoyStat = st[[4]][["Roy"]])
}
cen <- list(
  `A normal, cov iguais, balanceado`          = list(ns = rep(70, g), Sig = rep(list(Sp), g), gl = Inf),
  `B t(5), cov iguais, balanceado`            = list(ns = rep(70, g), Sig = rep(list(Sp), g), gl = 5),
  `C normal, cov distintas, desbalanceado`    = list(ns = ns_het,     Sig = S_g,              gl = Inf),
  `D normal, cov REAIS, balanceado (observado)` = list(ns = rep(70, g), Sig = S_g,              gl = Inf))

n_cel <- length(cen) * length(deltas)
cat("\nMonte Carlo: R =", R, "repeticoes por celula (", n_cel,
    "celulas); erro padrao em 0,05 =", round(sqrt(0.05 * 0.95 / R), 4), "\n")

mc_raw <- lapply(cen, function(cc)
  lapply(deltas, function(dl) replicate(R, rejeita(gera(cc$ns, mus_de(dl), cc$Sig, cc$gl)))))
mc <- lapply(mc_raw, function(sc_) sapply(sc_, function(m) rowMeans(m[1:4, , drop = FALSE])))
# Roy com tamanho corrigido: corte = quantil 95% de lambda_1 sob H0 do proprio cenario
roy_corr <- lapply(mc_raw, function(sc_) {
  cut <- quantile(sc_[["H0"]]["RoyStat", ], 0.95)
  sapply(sc_, function(m) mean(m["RoyStat", ] > cut))
})
cat("\nTaxa de rejeicao (nivel 0,05). Linhas: H0 = erro tipo I; d = poder\n")
for (nm in names(mc)) {
  cat("\nCenario", nm, "\n")
  tb <- cbind(t(mc[[nm]]), `Roy (corrigido)` = roy_corr[[nm]])
  print(round(tb, 3))
}

tit <- c("A: normal, covari\u00e2ncias iguais, balanceado",
         "B: t(5), covari\u00e2ncias iguais, balanceado",
         "C: normal, covari\u00e2ncias distintas, desbalanceado",
         "D: normal, covari\u00e2ncias reais, balanceado (observado)")
leg_lab <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy (aprox. F, n\u00edvel inflado)",
             "Roy (tamanho corrigido)")
pdf(file.path(fig_dir, "poder.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(2, 0, 0, 0))
for (i in seq_along(mc)) {
  matplot(deltas, cbind(t(mc[[i]]), roy_corr[[i]]), type = "b", pch = 1:5, lty = 1:5,
          col = c(1:4, "purple"), ylim = c(0, 1),
          xlab = "delta (fra\u00e7\u00e3o do efeito observado)",
          ylab = "Taxa de rejei\u00e7\u00e3o", main = tit[i], cex.main = 0.9)
  abline(h = 0.05, lty = 3)
  legend("bottomright", leg_lab, pch = 1:5, lty = 1:5, col = c(1:4, "purple"), cex = 0.65)
}
mtext("Nota: Roy pela aprox. F tem erro tipo I inflado; a curva roxa usa o corte simulado sob H0 de cada cen\u00e1rio.",
      side = 1, outer = TRUE, line = 0.5, cex = 0.8)
invisible(dev.off())

# 10.1 Calibracao de Roy: valor critico da aproximacao F vs. distribuicao nula simulada
df1 <- max(p, g - 1); df2 <- (n - g) - df1 + (g - 1)   # gl usados por summary() para Roy
crit_F <- qf(0.95, df1, df2) * df1 / df2                # corte de lambda1 implicito na aprox. F
R_roy <- 2000
lam1_nula <- replicate(R_roy, {
  d <- gera(rep(70, g), mus_de(0), rep(list(Sp), g))
  unname(est4(d$Y, d$g)["Roy"])
})
cat("\nCalibracao de Roy (cenario A, H0):\n",
    "corte de lambda1 pela aprox. F (5%):", round(crit_F, 4), "\n",
    "quantil 95% simulado de lambda1:   ", round(quantile(lam1_nula, 0.95), 4), "\n",
    "nivel real da regra aprox. F:      ", round(mean(lam1_nula > crit_F), 4), "\n")