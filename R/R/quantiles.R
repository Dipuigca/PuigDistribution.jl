# Cuantiles, intervalos, entropía y muestreo.

.puig_quantile_scalar <- function(p, lam, k, T, tol = 1e-10, maxit = 200) {
  if (p == 0) return(0.0)
  if (p == 1) return(Inf)
  lo <- 0.0
  hi <- max(lam, 1.0) * 2.0
  while (.puig_cumulative(hi, lam, k, T) < p) hi <- hi * 2.0
  for (i in seq_len(as.integer(maxit))) {
    mid <- 0.5 * (lo + hi)
    if (.puig_cumulative(mid, lam, k, T) < p) {
      lo <- mid
    } else {
      hi <- mid
    }
    if ((hi - lo) < tol) break
  }
  0.5 * (lo + hi)
}

.puig_quantile <- function(p, lam, k, T, tol = 1e-10, maxit = 200) {
  if (length(p) == 1L) return(.puig_quantile_scalar(p, lam, k, T, tol, maxit))
  vapply(p, function(pi) .puig_quantile_scalar(pi, lam, k, T, tol, maxit),
         numeric(1))
}

#' Cuantil \eqn{Q(p)} tal que \eqn{F(Q(p)) = p}.
#'
#' Búsqueda por bisección sobre la acumulada; vectorizado sobre \code{p}.
#'
#' @param p Probabilidades en \eqn{[0, 1]} (escalar o vector).
#' @param lam Parámetro \eqn{\lambda \ge 0}, o un objeto \code{PuigDistribution}
#'   o \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @param tol Tolerancia de longitud del intervalo de bisección.
#' @param maxit Número máximo de iteraciones.
#' @return Cuantil \eqn{Q(p)} (escalar o vector). \code{Q(0) = 0},
#'   \code{Q(1) = +\infty}.
#' @export
Puig_quantile <- function(p, lam, k = NULL, T = NULL, tol = 1e-10, maxit = 200) {
  p <- as.numeric(p)
  if (any(p < 0 | p > 1)) stop("p debe estar en [0, 1]")
  if (is.PuigDistribution(lam)) {
    return(.puig_quantile(p, lam$lam, lam$k, lam$T, tol, maxit))
  }
  if (is.MB(lam)) {
    return(.puig_quantile(p, 0, lam$k, 1.0 / lam$B, tol, maxit))
  }
  .puig_quantile(p, lam, k, T, tol, maxit)
}

#' Intervalo central al 98\%: \eqn{(P_1, P_{99})}.
#'
#' @param lam Parámetro \eqn{\lambda \ge 0}, o un objeto \code{PuigDistribution}
#'   o \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Vector con nombres \code{lo}/\code{hi}.
#' @export
Puig_ci98 <- function(lam, k = NULL, T = NULL) {
  q <- Puig_quantile(c(0.01, 0.99), lam, k, T)
  c(lo = q[1], hi = q[2])
}

#' Intervalo central al 99.8\%: \eqn{(P_{0.1}, P_{99.9})}.
#'
#' @param lam Parámetro \eqn{\lambda \ge 0}, o un objeto \code{PuigDistribution}
#'   o \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Vector con nombres \code{lo}/\code{hi}.
#' @export
Puig_ci996 <- function(lam, k = NULL, T = NULL) {
  q <- Puig_quantile(c(0.001, 0.999), lam, k, T)
  c(lo = q[1], hi = q[2])
}

.puig_entropy <- function(lam, k, T, N = 400) {
  ci <- Puig_quantile(c(0.001, 0.999), lam, k, T)
  init <- max(ci[1] / 10.0, 1e-6)
  fin <- ci[2]
  x <- seq(init, fin, length.out = as.integer(N))
  f <- .puig_pdf(x, lam, k, T)
  lf <- .puig_logpdf(x, lam, k, T)
  Fv <- -f * lf
  h <- (fin - init) / (as.integer(N) - 1)
  Fv[1] <- Fv[1] / 2.0
  Fv[length(Fv)] <- Fv[length(Fv)] / 2.0
  sum(Fv * h)
}

#' Entropía diferencial \eqn{H(X) = -\int f(x)\log f(x)\,dx}.
#'
#' Integración numérica por la regla del trapecio sobre el intervalo
#' \eqn{[P_{0.1}, P_{99.9}]}.
#'
#' @param dist Objeto \code{PuigDistribution} o \code{MB}.
#' @param N Número de puntos de la rejilla de integración (por defecto 400).
#' @return Entropía diferencial (escalar, en nats).
#' @export
Puig_entropy <- function(dist, N = 400) {
  if (is.MB(dist)) return(.puig_entropy(0, dist$k, 1.0 / dist$B, N))
  if (is.PuigDistribution(dist)) return(.puig_entropy(dist$lam, dist$k, dist$T, N))
  stop("Puig_entropy requiere un objeto PuigDistribution o MB")
}

.puig_rand <- function(n, lam, k, T) {
  n <- as.integer(n)
  if (lam == 0) {
    W <- rgamma(n, shape = k / 2.0, scale = 2.0 / T)
    return(sqrt(W / T))
  }
  W <- suppressWarnings(rchisq(n, df = k, ncp = lam^2 * T))
  sqrt(W / T)
}

#' Muestreo pseudoaleatorio exacto.
#'
#' Representación exacta: \eqn{X = \sqrt{W/T}} con \eqn{W \sim \chi^2(k,
#' \lambda^2 T)}; para \eqn{\lambda = 0}, \eqn{W \sim \operatorname{Gamma}(k/2,
#' 2/T)}.
#'
#' @param n Número de muestras.
#' @param lam Parámetro \eqn{\lambda \ge 0}, o un objeto \code{PuigDistribution}
#'   o \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Vector numérico de longitud \code{n} con \eqn{X \ge 0}.
#' @export
Puig_rand <- function(n, lam, k = NULL, T = NULL) {
  n <- as.integer(n)
  if (is.PuigDistribution(lam)) return(.puig_rand(n, lam$lam, lam$k, lam$T))
  if (is.MB(lam)) return(.puig_rand(n, 0, lam$k, 1.0 / lam$B))
  .puig_rand(n, lam, k, T)
}