# Validación del port R contra los golden files de Julia.
# Uso: Rscript R/tests/run_golden_tests.R
# (debe invocarse desde la raíz del repositorio)

rfiles <- list.files("R", pattern = "\\.R$", full.names = TRUE)
rfiles <- rfiles[!grepl("tests", rfiles)]
for (f in rfiles) source(f)

golden <- file.path("golden")

PASS <- 0L
FAIL <- 0L

check <- function(label, py, gl, rel_tol = 1e-9, abs_tol = 1e-14) {
  g <- as.numeric(gl); p <- as.numeric(py)
  if (is.infinite(g) || is.infinite(p)) {
    ok <- isTRUE(all.equal(p, g))
  } else if (g == 0 && p == 0) {
    ok <- TRUE
  } else if (g == 0) {
    ok <- abs(p) < abs_tol
  } else if (abs(g) > 1e-6) {
    ok <- abs(p - g) / abs(g) < rel_tol
  } else {
    ok <- abs(p - g) < abs_tol
  }
  assign("PASS", PASS + ok, envir = .GlobalEnv)
  assign("FAIL", FAIL + !ok, envir = .GlobalEnv)
  if (!ok) {
    cat(sprintf("  FAIL %s: R=% .16g golden=% .16g err=%.2e\n",
                label, p, g, abs(p - g)))
  }
  invisible(NULL)
}

read_golden <- function(name) {
  read.csv(file.path(golden, name), fileEncoding = "UTF-8", check.names = FALSE)
}

cat("=== PDF ===\n")
pdf <- read_golden("golden_pdf.csv")
for (i in seq_len(nrow(pdf))) {
  lam <- pdf[[1]][i]; kv <- pdf[[2]][i]; Tv <- pdf[[3]][i]; xv <- pdf[[4]][i]
  p_asymp <- Puig_pdf(xv, lam, kv, Tv, method = "asymp")
  check(sprintf("pdf_asymp lam=%g k=%g T=%g x=%g", lam, kv, Tv, xv),
        p_asymp, pdf[[5]][i], rel_tol = 1e-10, abs_tol = 1e-15)
  if (lam > 0) {
    p_arb <- Puig_pdf(xv, lam, kv, Tv, method = "arb")
    check(sprintf("pdf_arb lam=%g k=%g T=%g x=%g", lam, kv, Tv, xv),
          p_arb, pdf[[6]][i], abs_tol = 1e-12)
  }
}

cat("=== SURV ===\n")
surv <- read_golden("golden_surv.csv")
for (i in seq_len(nrow(surv))) {
  lam <- surv[[1]][i]; kv <- surv[[2]][i]; Tv <- surv[[3]][i]; xv <- surv[[4]][i]
  p <- Puig_surviving(xv, lam, kv, Tv)
  check(sprintf("surv lam=%g k=%g T=%g x=%g", lam, kv, Tv, xv),
        p, surv[[5]][i], rel_tol = 5e-8, abs_tol = 1e-14)
}

cat("=== CDF ===\n")
cdf <- read_golden("golden_cdf.csv")
for (i in seq_len(nrow(cdf))) {
  lam <- cdf[[1]][i]; kv <- cdf[[2]][i]; Tv <- cdf[[3]][i]; xv <- cdf[[4]][i]
  p <- Puig_cumulative(xv, lam, kv, Tv)
  check(sprintf("cdf lam=%g k=%g T=%g x=%g", lam, kv, Tv, xv),
        p, cdf[[5]][i], rel_tol = 5e-8, abs_tol = 1e-14)
}

cat("=== MOMENTS ===\n")
mom <- read_golden("golden_moments.csv")
for (i in seq_len(nrow(mom))) {
  lam <- mom[[1]][i]; kv <- mom[[2]][i]; Tv <- mom[[3]][i]; n <- as.integer(mom[[4]][i])
  mu <- Puig_moments(lam, kv, Tv, n)
  check(sprintf("moments lam=%g k=%g T=%g n=%d", lam, kv, Tv, n),
        mu[n], mom[[5]][i], rel_tol = 1e-10)
}

cat("=== QUANTILE ===\n")
qg <- read_golden("golden_quantile.csv")
for (i in seq_len(nrow(qg))) {
  lam <- qg[[1]][i]; kv <- qg[[2]][i]; Tv <- qg[[3]][i]; pval <- qg[[4]][i]
  q <- Puig_quantile(pval, lam, kv, Tv)
  check(sprintf("quantile lam=%g k=%g T=%g p=%g", lam, kv, Tv, pval),
        q, qg[[5]][i], rel_tol = 1e-7, abs_tol = 1e-10)
}

cat("=== MB ===\n")
mb <- read_golden("golden_mb.csv")
for (i in seq_len(nrow(mb))) {
  kv <- mb[[1]][i]; Bv <- mb[[2]][i]; xv <- mb[[3]][i]
  check(sprintf("mb_pdf k=%g B=%g x=%g", kv, Bv, xv),
        MB_pdf(xv, kv, Bv), mb[[4]][i], rel_tol = 1e-10)
  check(sprintf("mb_surv k=%g B=%g x=%g", kv, Bv, xv),
        MB_surviving(xv, kv, Bv), mb[[5]][i], rel_tol = 1e-10)
  check(sprintf("mb_cdf k=%g B=%g x=%g", kv, Bv, xv),
        MB_cumulative(xv, kv, Bv), mb[[6]][i], rel_tol = 1e-10)
}
cat("=== MB MOMENTS ===\n")
mbm <- read_golden("golden_mb_moments.csv")
for (i in seq_len(nrow(mbm))) {
  kv <- mbm[[1]][i]; Bv <- mbm[[2]][i]; n <- as.integer(mbm[[3]][i])
  mu <- Puig_moments(0.0, kv, 1.0 / Bv, n)
  check(sprintf("mb_mom k=%g B=%g n=%d", kv, Bv, n),
        mu[n], mbm[[4]][i], rel_tol = 1e-10)
}

cat(sprintf("\n%s\nRESULTADOS: %d PASS, %d FAIL\n",
            paste(rep("=", 60), collapse = ""), PASS, FAIL))
if (FAIL > 0) quit(status = 1)
cat("Todos los tests pasaron (OK)\n")