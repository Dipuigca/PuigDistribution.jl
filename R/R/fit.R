# Métricas de bondad de ajuste + dimensión efectiva + supervivencia empírica.

calc_r2 <- function(S_emp, S_mod) {
  S_emp <- as.numeric(S_emp); S_mod <- as.numeric(S_mod)
  1.0 - sum((S_emp - S_mod)^2) / sum((S_emp - mean(S_emp))^2)
}

calc_mae <- function(S_emp, S_mod) {
  mean(abs(as.numeric(S_emp) - as.numeric(S_mod)))
}

calc_rmse <- function(S_emp, S_mod) {
  sqrt(mean((as.numeric(S_emp) - as.numeric(S_mod))^2))
}

calc_maxae <- function(S_emp, S_mod) {
  max(abs(as.numeric(S_emp) - as.numeric(S_mod)))
}

calc_iae <- function(x, S_emp, S_mod) {
  x <- as.numeric(x)
  diffs <- abs(as.numeric(S_emp) - as.numeric(S_mod))
  dx <- diff(x)
  sum(dx * (diffs[-length(diffs)] + diffs[-1]) / 2.0)
}

#' Todas las métricas (R\u00b2, MAE, RMSE, MaxAE, IAE).
#' @export
calc_all_metrics <- function(x, S_emp, S_mod) {
  new_PuigStats(
    calc_r2(S_emp, S_mod),
    calc_mae(S_emp, S_mod),
    calc_rmse(S_emp, S_mod),
    calc_maxae(S_emp, S_mod),
    calc_iae(x, S_emp, S_mod)
  )
}

#' Dimensión efectiva k = d\u00b2/(d + 2\u00b7\u03a3_{i<j} R_ij\u00b2).
#' @export
effective_dimension <- function(M) {
  M <- as.matrix(M)
  d <- ncol(M)
  if (d <= 1) return(1.0)
  R <- cor(M)
  iu <- upper.tri(R)
  s_off <- sum(R[iu]^2)
  denom <- d + 2.0 * s_off
  if (denom == 0) return(1.0)
  as.numeric(d^2) / denom
}

#' Curva de supervivencia empírica sobre una rejilla.
#' @export
empirical_survival <- function(times, grid) {
  times <- as.numeric(times); grid <- as.numeric(grid)
  n <- length(times)
  if (n == 0) return(rep(0, length(grid)))
  vapply(grid, function(g) sum(times >= g) / n, numeric(1))
}

# ---- Funciones objetivo ----

.err_l1 <- function(B, x_opt, S_emp, lam, k_target) {
  if (B <= 0) return(1e12)
  T_val <- 1.0 / B^2
  S_mod <- tryCatch(as.numeric(.puig_surviving(x_opt, lam, k_target, T_val)),
                    error = function(e) rep(1e6, length(S_emp)))
  if (!all(is.finite(S_mod))) return(1e12)
  sum(abs(S_emp - S_mod))
}

.err_l2_matrix <- function(par, x_opt, S_emp, lam) {
  k_val <- par[1]; B_val <- par[2]
  if (k_val < 1.0 || B_val <= 0.0) return(1e12)
  T_val <- 1.0 / B_val^2
  S_mod <- tryCatch(as.numeric(.puig_surviving(x_opt, lam, k_val, T_val)),
                    error = function(e) rep(1e6, length(S_emp)))
  if (!all(is.finite(S_mod))) return(1e12)
  sum((S_emp - S_mod)^2)
}

._fit_matrix <- function(Data, times = NULL, method = "dynamic") {
  Data <- as.matrix(Data)
  n_obs <- nrow(Data)
  d <- ncol(Data)
  if (n_obs < 10) {
    stop(sprintf("Se requieren al menos 10 observaciones para el ajuste (recibidas: %d)", n_obs))
  }

  t_raw <- if (is.null(times)) sqrt(rowSums(Data^2)) else as.numeric(times)
  if (length(t_raw) != n_obs) stop("times debe tener tantas observaciones como filas de Data")

  valid <- is.finite(t_raw) & (t_raw >= 0)
  t <- t_raw[valid]
  if (length(t) < 10) stop("Se requieren al menos 10 observaciones válidas para el ajuste")

  Dv <- Data[valid, , drop = FALSE]
  col_means <- colMeans(Dv)
  lam <- sqrt(sum(col_means^2))
  k_corr <- effective_dimension(Dv)

  min_t <- min(t); max_t <- max(t)
  x_opt <- seq(min_t, max_t, length.out = 150)
  S_emp_opt <- empirical_survival(t, x_opt)

  sdt <- sd(t)
  B0 <- if (is.finite(sdt) && sdt > 0) sdt else 1.0

  k_est <- k_corr
  B_est <- B0

  if (method %in% c("dynamic", "fixed")) {
    k_target <- if (method == "fixed") as.numeric(d) else k_corr
    res <- optimize(function(B) .err_l1(B, x_opt, S_emp_opt, lam, k_target),
                    interval = c(0.01, 5.0 * B0), tol = 1e-12)
    B_est <- res$minimum
    k_est <- k_target
  } else if (method == "numerical") {
    upper_k <- max(3.0, as.numeric(d))
    init_par <- c(min(max(k_corr, 1.0), upper_k), B0)
    res <- optim(init_par,
                 function(par) .err_l2_matrix(par, x_opt, S_emp_opt, lam),
                 method = "L-BFGS-B",
                 lower = c(1.0, 0.001),
                 upper = c(upper_k, 10.0 * B0),
                 control = list(maxit = 300))
    if (res$convergence == 0 || is.finite(res$value) && res$value < 1e12) {
      k_est <- res$par[1]; B_est <- res$par[2]
    } else {
      k_est <- min(max(k_corr, 1.0), upper_k); B_est <- B0
    }
  } else {
    stop(sprintf("M\u00e9todo desconocido :%s. Opciones v\u00e1lidas: :dynamic, :numerical, :fixed", method))
  }

  T_est <- 1.0 / B_est^2
  x_eval <- seq(min_t, max_t, length.out = 500)
  S_emp_eval <- empirical_survival(t, x_eval)
  S_mod_eval <- as.numeric(.puig_surviving(x_eval, lam, k_est, T_est))
  fit_metrics <- calc_all_metrics(x_eval, S_emp_eval, S_mod_eval)

  new_PuigFitResult(PuigDistribution(lam, k_est, T_est), t, fit_metrics)
}

._fit_dataframe <- function(df, vars = NULL, time_col = NULL, method = "dynamic") {
  if (is.null(vars)) {
    exclude <- if (is.null(time_col)) character(0) else as.character(time_col)
    vars <- colnames(df)[vapply(df, is.numeric, logical(1)) &
                           !(colnames(df) %in% exclude)]
  }
  if (length(vars) == 0) stop("Se requiere al menos una columna num\u00e9rica para el ajuste")
  mat <- as.matrix(df[, vars, drop = FALSE])
  if (ncol(mat) < 2) stop("vars debe seleccionar m\u00e1s de una columna num\u00e9rica")
  times <- if (!is.null(time_col)) as.numeric(df[[time_col]]) else NULL
  ._fit_matrix(mat, times = times, method)
}

.err_2d_fixed <- function(p, x_opt, S_emp, k_target) {
  lam_val <- exp(p[1])
  B_val <- exp(p[2])
  T_val <- 1.0 / B_val^2
  S_mod <- tryCatch(as.numeric(.puig_surviving(x_opt, lam_val, k_target, T_val)),
                    error = function(e) rep(1e6, length(S_emp)))
  if (!all(is.finite(S_mod))) return(1e12)
  sum((S_emp - S_mod)^2)
}

.err_3d <- function(p, x_opt, S_emp) {
  lam_val <- exp(p[1])
  k_val <- 1.0 + exp(p[2])
  if (k_val > 10.0) return(1e12)
  B_val <- exp(p[3])
  T_val <- 1.0 / B_val^2
  S_mod <- tryCatch(as.numeric(.puig_surviving(x_opt, lam_val, k_val, T_val)),
                    error = function(e) rep(1e6, length(S_emp)))
  if (!all(is.finite(S_mod))) return(1e12)
  sum((S_emp - S_mod)^2)
}

._fit_vector <- function(Data, method = "dynamic", k_fixed = NULL) {
  Data <- as.numeric(Data)
  valid <- is.finite(Data) & (Data >= 0)
  t <- Data[valid]
  n_obs <- length(t)
  if (n_obs < 10) {
    stop(sprintf("Se requieren al menos 10 observaciones para el ajuste (recibidas: %d)", n_obs))
  }

  min_t <- min(t); max_t <- max(t)
  x_opt <- seq(min_t, max_t, length.out = 60)
  S_emp_opt <- empirical_survival(t, x_opt)

  m1 <- mean(t)
  m2 <- mean(t^2)
  s <- sd(t)
  B0 <- if (is.finite(s) && s > 0) s else 1.0
  lam0 <- sqrt(max(0.0, m1^2 - s^2))
  k0 <- min(max((m2 - lam0^2) / B0^2, 1.0), 10.0)

  is_fixed_k <- method == "fixed" || !is.null(k_fixed)
  k_target <- if (!is.null(k_fixed)) as.numeric(k_fixed) else 1.0
  if (k_target < 1.0) {
    stop(sprintf("k_fixed debe ser >= 1.0 (recibido: %g)", k_target))
  }

  k_est <- 1.0; lam_est <- lam0; B_est <- B0

  if (is_fixed_k) {
    k_est <- k_target
    init_p <- c(log(max(1e-3, lam0)), log(max(1e-3, B0)))
    res <- optim(init_p,
                 function(p) .err_2d_fixed(p, x_opt, S_emp_opt, k_target),
                 method = "Nelder-Mead",
                 control = list(maxit = 500, reltol = 1e-10))
    p_opt <- res$par
    lam_est <- exp(p_opt[1]); B_est <- exp(p_opt[2])
  } else if (method %in% c("dynamic", "numerical")) {
    init_p <- c(log(max(1e-3, lam0)), log(max(1e-3, k0 - 1.0)),
                log(max(1e-3, B0)))
    res <- optim(init_p,
                 function(p) .err_3d(p, x_opt, S_emp_opt),
                 method = "Nelder-Mead",
                 control = list(maxit = 600, reltol = 1e-10))
    p_opt <- res$par
    lam_est <- exp(p_opt[1])
    k_est <- 1.0 + exp(p_opt[2])
    B_est <- exp(p_opt[3])
  } else {
    stop(sprintf("M\u00e9todo desconocido :%s. Opciones v\u00e1lidas: :dynamic, :numerical, :fixed", method))
  }

  T_est <- 1.0 / B_est^2
  x_eval <- seq(min_t, max_t, length.out = 500)
  S_emp_eval <- empirical_survival(t, x_eval)
  S_mod_eval <- as.numeric(.puig_surviving(x_eval, lam_est, k_est, T_est))
  fit_metrics <- calc_all_metrics(x_eval, S_emp_eval, S_mod_eval)

  new_PuigFitResult(PuigDistribution(lam_est, k_est, T_est), t, fit_metrics)
}

#' Ajusta datos a la distribución de Puig.
#'
#' Acepta `data.frame`, `matrix` (filas = observaciones) o vector numérico 1D.
#' Métodos: `"dynamic"`, `"fixed"`, `"numerical"`; `k_fixed` para vectores.
#' @export
Puig_fit <- function(Data, vars = NULL, time_col = NULL, times = NULL,
                     method = "dynamic", k_fixed = NULL) {
  if (is.data.frame(Data)) {
    if (!is.null(k_fixed)) stop("k_fixed solo aplica al ajuste de vectores 1D")
    return(._fit_dataframe(Data, vars = vars, time_col = time_col, method))
  }
  if (is.matrix(Data)) {
    if (!is.null(k_fixed)) stop("k_fixed solo aplica al ajuste de vectores 1D")
    return(._fit_matrix(Data, times = times, method))
  }
  if (is.numeric(Data) && is.null(dim(Data))) {
    return(._fit_vector(Data, method, k_fixed))
  }
  stop("Entrada no soportada: use data.frame, matrix o vector num\u00e9rico")
}