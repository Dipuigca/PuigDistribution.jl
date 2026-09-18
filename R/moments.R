# Momentos y estadísticos teóricos (Laguerre fraccionaria vía Rmpfr).

._puig_mean <- function(lam, k, T) {
  ls <- lam * sqrt(T)
  sqrt(pi / (2.0 * T)) * laguerre_real(0.5, k / 2.0 - 1.0, -ls^2 / 2.0)
}

._puig_moment3 <- function(lam, k, T) {
  ls <- lam * sqrt(T)
  L <- laguerre_real(1.5, k / 2.0 - 1.0, -ls^2 / 2.0)
  (3.0 / T) * sqrt(pi / (2.0 * T)) * L
}

._puig_var <- function(lam, k, T) {
  mu2 <- k / T + lam^2
  mu2 - ._puig_mean(lam, k, T)^2
}

._puig_std <- function(lam, k, T) sqrt(._puig_var(lam, k, T))

._puig_skewness <- function(lam, k, T) {
  mu1 <- ._puig_mean(lam, k, T)
  mu2 <- ._puig_var(lam, k, T) + mu1^2
  mu3 <- ._puig_moment3(lam, k, T)
  sigma3 <- ._puig_std(lam, k, T)^3
  (mu3 - 3.0 * mu1 * mu2 + 2.0 * mu1^3) / sigma3
}

._puig_kurtosis <- function(lam, k, T) {
  # Kurtosis NO restada (γ2 = E[(X-μ)^4]/σ^4)
  mu1 <- ._puig_mean(lam, k, T)
  mu2 <- k / T + lam^2
  mu4 <- mu2^2 + 2.0 * k / T^2 + 4.0 * lam^2 / T
  mu3 <- ._puig_moment3(lam, k, T)
  sigma4 <- ._puig_var(lam, k, T)^2
  (mu4 - 4.0 * mu1 * mu3 + 6.0 * mu1^2 * mu2 - 3.0 * mu1^4) / sigma4
}

._puig_moments_raw <- function(lam, k, T, n) {
  n <- as.integer(n)
  if (n < 1) stop(sprintf("n debe ser \u2265 1, recibido n=%d", n))
  mu <- numeric(n)
  mu[1] <- ._puig_mean(lam, k, T)
  if (n == 1) return(mu)
  mu[2] <- k / T + lam^2
  if (n == 2) return(mu)
  mu[3] <- ._puig_moment3(lam, k, T)
  if (n == 3) return(mu)
  mu[4] <- mu[2]^2 + 2.0 * k / T^2 + 4.0 * lam^2 / T
  if (n == 4) return(mu)

  invT <- 1.0 / T
  invT2 <- invT * invT
  lam2 <- lam * lam
  for (i in 5:n) {
    order_v <- i
    mu[i] <- ((2.0 * order_v - 4.0 + k) * invT + lam2) * mu[i - 2] -
      (order_v - 2.0) * (order_v + k - 4.0) * invT2 * mu[i - 4]
  }
  mu
}

._puig_moments <- function(lam, k, T, n = 4) {
  mu <- ._puig_moments_raw(lam, k, T, as.integer(n))
  if (as.integer(n) == 4) {
    names(mu) <- c("\u03bc1", "\u03bc2", "\u03bc3", "\u03bc4")
  }
  mu
}

._puig_stats <- function(lam, k, T) {
  mu <- ._puig_moments_raw(lam, k, T, 4)
  mu1 <- mu[1]; mu2 <- mu[2]; mu3 <- mu[3]; mu4 <- mu[4]
  sig2 <- mu2 - mu1^2
  sig <- sqrt(sig2)
  sk <- (mu3 - 3.0 * mu1 * mu2 + 2.0 * mu1^3) / sig^3
  kur <- (mu4 - 4.0 * mu1 * mu3 + 6.0 * mu1^2 * mu2 - 3.0 * mu1^4) / sig2^2
  list(mean = mu1, var = sig2, sig = sig,
       skewness = sk, kurtosis = kur)
}

# ---- API pública ----

#' Media teórica E[X] = \u221a(\u03c0/(2T)) L_{1/2}^{(k/2-1)}(-\u03bb\u00b2T/2).
#' @export
Puig_mean <- function(lam, k = NULL, T = NULL) {
  if (is.MB(lam)) return(._puig_mean(0, lam$k, 1.0 / lam$B))
  if (is.PuigDistribution(lam)) return(._puig_mean(lam$lam, lam$k, lam$T))
  ._puig_mean(lam, k, T)
}

#' Varianza teórica.
#' @export
Puig_var <- function(lam, k = NULL, T = NULL) {
  if (is.MB(lam)) return(._puig_var(0, lam$k, 1.0 / lam$B))
  if (is.PuigDistribution(lam)) return(._puig_var(lam$lam, lam$k, lam$T))
  ._puig_var(lam, k, T)
}

#' Desviación estándar teórica.
#' @export
Puig_std <- function(lam, k = NULL, T = NULL) {
  if (is.MB(lam)) return(._puig_std(0, lam$k, 1.0 / lam$B))
  if (is.PuigDistribution(lam)) return(._puig_std(lam$lam, lam$k, lam$T))
  ._puig_std(lam, k, T)
}

#' Asimetría estandarizada \u03b3\u2081.
#' @export
Puig_skewness <- function(lam, k = NULL, T = NULL) {
  if (is.MB(lam)) return(._puig_skewness(0, lam$k, 1.0 / lam$B))
  if (is.PuigDistribution(lam)) return(._puig_skewness(lam$lam, lam$k, lam$T))
  ._puig_skewness(lam, k, T)
}

#' Kurtosis (no restada) \u03b3\u2082.
#' @export
Puig_kurtosis <- function(lam, k = NULL, T = NULL) {
  if (is.MB(lam)) return(._puig_kurtosis(0, lam$k, 1.0 / lam$B))
  if (is.PuigDistribution(lam)) return(._puig_kurtosis(lam$lam, lam$k, lam$T))
  ._puig_kurtosis(lam, k, T)
}

#' Momentos raw \u03bc\u2081..\u03bc\u2099. Con n=4 devuelve vector con nombre.
#' @export
Puig_moments <- function(lam, k = NULL, T = NULL, n = NULL) {
  if (is.MB(lam)) {
    params <- lam
    n_val <- if (!is.null(k)) as.integer(k) else if (!is.null(n)) as.integer(n) else 4L
    return(._puig_moments(0, params$k, 1.0 / params$B, n_val))
  }
  if (is.PuigDistribution(lam)) {
    params <- lam
    n_val <- if (!is.null(k)) as.integer(k) else if (!is.null(n)) as.integer(n) else 4L
    return(._puig_moments(params$lam, params$k, params$T, n_val))
  }
  if (is.null(k) || is.null(T)) stop("Se requieren lam, k, T")
  n_val <- if (!is.null(n)) as.integer(n) else 4L
  ._puig_moments(lam, k, T, n_val)
}

#' Estadísticos: mean, var, sig, skewness, kurtosis.
#' @export
Puig_stats <- function(lam, k = NULL, T = NULL) {
  if (is.MB(lam)) return(._puig_stats(0, lam$k, 1.0 / lam$B))
  if (is.PuigDistribution(lam)) return(._puig_stats(lam$lam, lam$k, lam$T))
  ._puig_stats(lam, k, T)
}