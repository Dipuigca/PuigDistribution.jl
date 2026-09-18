# Tests unitarios del port R (port de test/runtests.jl del paquete Julia).
# Uso: Rscript R/tests/run_unit_tests.R

rfiles <- list.files("R", pattern = "\\.R$", full.names = TRUE)
rfiles <- rfiles[!grepl("tests", rfiles)]
for (f in rfiles) source(f)

set.seed(0)
PASS <- 0L; FAIL <- 0L

ok <- function(cond, label) {
  if (isTRUE(cond)) {
    assign("PASS", get("PASS", .GlobalEnv) + 1L, envir = .GlobalEnv)
  } else {
    assign("FAIL", get("FAIL", .GlobalEnv) + 1L, envir = .GlobalEnv)
    cat(sprintf("  FAIL: %s\n", label))
  }
  invisible(NULL)
}
expect_error <- function(expr, label) {
  tries <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  ok(tries, label)
}
expect_near <- function(a, b, label, tol = 1e-9) {
  ok(is.numeric(a) && is.numeric(b) && isTRUE(all.equal(a, b, tolerance = tol)),
     sprintf("%s (a=%.10g b=%.10g)", label, as.numeric(a), as.numeric(b)))
}

# ---- 0. Tipos ----
d <- PuigDistribution(3.0, 2.0, 1.0)
ok(is.PuigDistribution(d), "d es PuigDistribution")
ok(d$lam == 3.0 && d$k == 2.0 && d$T == 1.0, "campos d")
mb <- MB(2.0, 1.0)
ok(is.MB(mb) && mb$k == 2.0 && mb$B == 1.0, "campos MB")

# ---- 1. effective_dimension ----
ok(effective_dimension(matrix(rnorm(50), ncol = 1)) == 1.0, "effd 1D")
set.seed(123)
M_uncorr <- matrix(rnorm(30000), ncol = 3)
ok(effective_dimension(M_uncorr) > 2.9 && effective_dimension(M_uncorr) < 3.1,
   "effd incorreladas ~3")
xb <- rnorm(100)
M_corr <- cbind(xb, xb, xb)
ok(effective_dimension(M_corr) < 1.0001, "effd correladas ~1")
k_exacto <- 3^2 / (3 + 2 * (0.8^2 + 0.6^2 + 0.7^2))
ok(abs(k_exacto - 1.5050167224) < 1e-6, "k formula exacta")

# ---- 2. empirical_survival ----
t <- c(1, 2, 3, 4, 5); grid <- c(0, 1, 3, 5, 6)
s <- empirical_survival(t, grid)
ok(all(c(s[1] == 1, s[2] == 1, abs(s[3] - 0.6) < 1e-12,
         abs(s[4] - 0.2) < 1e-12, s[5] == 0)), "emp_survival valores")

# ---- 3. Métricas ----
S_true <- c(1, 0.8, 0.5, 0.2, 0); xg <- 1:5
st <- calc_all_metrics(xg, S_true, S_true)
ok(abs(st$R2 - 1) < 1e-12 && st$MAE == 0 && st$RMSE == 0 && st$MaxAE == 0 && st$IAE == 0,
   "metricas perfectas")
st2 <- calc_all_metrics(xg, S_true, S_true + 0.05)
ok(abs(st2$MAE - 0.05) < 1e-12 && abs(st2$MaxAE - 0.05) < 1e-12 && st2$R2 < 1,
   "metricas discrepancia")

# ---- 4. Puig_fit (DataFrame / Matrix / Vector) ----
set.seed(42)
N <- 150
X1 <- 150 + rnorm(N) * 20
X2 <- 100 + 0.8 * X1 + rnorm(N) * 10
X3 <- 200 + 0.5 * X1 + rnorm(N) * 15
df <- data.frame(Squat = X1, Bench = X2, Deadlift = X3)

res_dyn <- Puig_fit(df, method = "dynamic")
ok(inherits(res_dyn, "PuigFitResult"), "res es PuigFitResult")
ok(length(res_dyn$indicator) == N, "len indicador")
ok(res_dyn$params$lam > 0 && res_dyn$params$T > 0, "lam,T > 0")
ok(res_dyn$params$k >= 1 && res_dyn$params$k <= 3, "k en [1,3]")
ok(res_dyn$stats$R2 > 0.95, "R2 dynamic > 0.95")

res_fix <- Puig_fit(df, method = "fixed")
ok(abs(res_fix$params$k - 3) < 1e-9 && res_fix$stats$R2 > 0.90, "fixed k=3")
res_num <- Puig_fit(df, method = "numerical")
ok(res_num$params$k >= 1 && res_num$params$k <= 3 && res_num$stats$R2 > 0.95, "numerical")

res_sub <- Puig_fit(df, vars = c("Squat", "Bench"))
ok(res_sub$params$k <= 2 && length(res_sub$indicator) == N, "vars subset")

df2 <- df
df2$Total <- sqrt(X1^2 + X2^2 + X3^2)
res_time <- Puig_fit(df2, vars = c("Squat", "Bench", "Deadlift"), time_col = "Total")
ok(max(abs(res_time$indicator - df2$Total)) < 1e-9, "time_col indicador")

res_mat <- Puig_fit(as.matrix(df), method = "dynamic")
ok(abs(res_mat$params$lam - res_dyn$params$lam) / res_dyn$params$lam < 1e-6 &&
     abs(res_mat$params$k - res_dyn$params$k) < 1e-6, "matrix == dynamic")

res_vec <- Puig_fit(df2$Total)
ok(res_vec$params$lam > 0 && res_vec$params$k >= 1 && res_vec$params$T > 0 &&
     res_vec$stats$R2 > 0.95, "vector 1D")

set.seed(77)
muestra <- Puig_rand(500, 40.0, 2.5, 0.05)
res_1d <- Puig_fit(muestra)
ok(abs(res_1d$params$lam - 40) < 5 && res_1d$params$k >= 1 && res_1d$stats$R2 > 0.98,
   sprintf("1D sintetico lam=%.2f R2=%.4f", res_1d$params$lam, res_1d$stats$R2))
res_1d_kf <- Puig_fit(muestra, k_fixed = 2.5)
ok(abs(res_1d_kf$params$k - 2.5) < 1e-6 && res_1d_kf$stats$R2 > 0.98, "1D k_fixed")

# ---- 5. Desempaquetado ----
set.seed(99)
dfe <- data.frame(A = 10 + rnorm(50), B = 20 + rnorm(50))
res <- Puig_fit(dfe)
p2 <- res$params; in2 <- res$indicator; st2 <- res$stats
ok(is.PuigDistribution(p2) && is.numeric(in2) && inherits(st2, "PuigStats"),
   "desempaquetado list")

# ---- 6. Errores ----
expect_error(Puig_fit(data.frame(A = runif(5), B = runif(5))), "fit n<10")
expect_error(Puig_fit(runif(5)), "fit vec n<10")
expect_error(Puig_fit(df, method = "metodo_inexistente"), "metodo desconocido")
expect_error(Puig_fit(data.frame(A = runif(20), B = runif(20)), method = "metodo_inexistente"),
             "metodo desconocido df")
expect_error(Puig_fit(runif(20), k_fixed = 0.5), "k_fixed < 1")
expect_error(Puig_fit(data.frame(Nombre = LETTERS[1:11])), "sin columnas numeric")

# ---- 7. Interfaz básica ----
ok(Puig_pdf(0, 3, 2, 1) == 0, "pdf(0)=0")
ok(Puig_logpdf(0, 3, 2, 1) == -Inf, "logpdf(0)=-Inf")
ok(Puig_cumulative(0, 3, 2, 1) == 0 && Puig_surviving(0, 3, 2, 1) == 1, "cdf/surv(0)")
ok(abs(Puig_cumulative(100, 3, 2, 1) - 1) < 1e-6, "cdf(100)~1")
ok(abs(Puig_surviving(100, 3, 2, 1)) < 1e-6, "surv(100)~0")
ok(abs(Puig_cumulative(5, 3, 2, 1) + Puig_surviving(5, 3, 2, 1) - 1) < 1e-10,
   "cdf+surv=1")
ok(Puig_quantile(0.01, 3, 2, 1) >= 0, "quantile>=0")
m <- Puig_mean(3, 2, 1); v <- Puig_var(3, 2, 1); st3 <- Puig_std(3, 2, 1)
ok(abs(v - st3^2) < 1e-6 && m > 0 && v > 0, "mean/var/std coherentes")
set.seed(42)
samples <- Puig_rand(100, 3, 2, 1)
ok(all(samples >= 0) && abs(mean(samples) - m) < 1.0, "rand coherente")

# ---- 8. Robustez asymp vs arb ----
for (lamv in c(10, 50, 100, 250, 500, 1000)) {
  for (fr in c(0.7, 0.9, 1.0, 1.1, 1.3)) {
    xv <- fr * lamv
    va <- Puig_pdf(xv, lamv, 2.4, 0.01, method = "asymp")
    vb <- Puig_pdf(xv, lamv, 2.4, 0.01, method = "arb")
    ok(if (vb > 1e-15) abs(va - vb) / vb < 1e-12 else abs(va - vb) < 1e-15,
       sprintf("asymp/arb lam=%g frac=%g", lamv, fr))
  }
}
for (kv in c(1.0001, 1.2, 2.0, 3.0, 5.0, 10.0)) {
  for (fr in c(0.8, 1.0, 1.2)) {
    xv <- fr * 50
    va <- Puig_pdf(xv, 50, kv, 0.05, method = "asymp")
    vb <- Puig_pdf(xv, 50, kv, 0.05, method = "arb")
    ok(if (vb > 1e-15) abs(va - vb) / vb < 1e-12 else abs(va - vb) < 1e-15,
       sprintf("asymp/arb k=%g frac=%g", kv, fr))
  }
}
###ok(Puig_pdf(-5, 100, 2.5, 0.05) == 0, "pdf negativo")  # no: debe fallar
expect_error(Puig_pdf(-5, 100, 2.5, 0.05), "pdf negativo error")

# ---- 9. Momentos ----
t4 <- Puig_moments(2.5, 3.2, 1.4)
v4 <- Puig_moments(2.5, 3.2, 1.4, 4)
v8 <- Puig_moments(2.5, 3.2, 1.4, 8)
ok(length(v4) == 4 && all(abs(v4 - t4) < 1e-12), "mom4 == tuple")
ok(length(v8) == 8, "mom8 len")
mu1 <- v4[1]; mu2 <- v4[2]; mu3 <- v4[3]; mu4 <- v4[4]
lamv <- 2.5; kv <- 3.2; Tv <- 1.4
mu5_rec <- ((10 - 4 + kv) / Tv + lamv^2) * mu3 - 3 * (1 + kv) / Tv^2 * mu1
ok(abs(v8[5] - mu5_rec) < 1e-10, "mu5 recurrencia")
mu6_rec <- ((12 - 4 + kv) / Tv + lamv^2) * mu4 - 4 * (2 + kv) / Tv^2 * mu2
ok(abs(v8[6] - mu6_rec) < 1e-10, "mu6 recurrencia")

Bmb <- 1.0; kmb <- 2.0
vec_mb <- Puig_moments(0, kmb, 1 / Bmb, 8)
for (nn in 1:8) {
  mu_an <- (2 * Bmb)^(nn / 2) * gamma((kmb + nn) / 2) / gamma(kmb / 2)
  ok(abs(vec_mb[nn] - mu_an) < 1e-10, sprintf("MB mom %d", nn))
}
expect_error(Puig_moments(2.5, 3.2, 1.4, 0), "n=0 no valido")
sts <- Puig_stats(2.5, 3.2, 1.4)
ok(abs(sts$mean - v4[1]) < 1e-12 && abs(sts$var - (v4[2] - v4[1]^2)) < 1e-12,
   "stats coherentes")

# ---- 10. Marcum-Q ----
mq <- marcumq_asymp(10, 15, 1)
ok(mq > 0 && mq < 1, "marcumq_asymp(10,15,1) en (0,1)")
sT <- sqrt(1); a_t <- 10 * sT
mq_lo <- marcumq(a_t, 25 * sT, 3.2 / 2)
ok(mq_lo < 1e-30, "cola lejana < 1e-30")
Sext <- Puig_surviving(c(20, 25, 30), 10, 3.2, 1)
ok(all(is.finite(Sext)) && all(Sext >= 0) && !is.unsorted(rev(Sext)), "colas monot")
Sbnd <- Puig_surviving(10 + 4, 10, 3.2, 1)
ok(is.finite(Sbnd) && Sbnd > 0 && Sbnd < 1, "frontera asymp")
ms <- marcumq(1, 3, 1)
ok(is.finite(ms) && ms > 0 && ms < 1, "ab<30 fallback ncx2")
ok(marcumq_asymp(10, 20, 2) > marcumq_asymp(10, 20, 1), "m mayor -> mayor Q")

cat(sprintf("\n%s\nRESULTADOS: %d PASS, %d FAIL\n",
            paste(rep("=", 60), collapse = ""),
            get("PASS", .GlobalEnv), get("FAIL", .GlobalEnv)))
if (get("FAIL", .GlobalEnv) > 0) quit(status = 1)
cat("Todos los tests pasaron (OK)\n")