# Wrappers estilo "Distribuciones" (convención R d/p/q/r).

#' Densidad de la distribución de Puig (convención \code{d/p/q/r}).
#'
#' @param x Puntos de evaluación (\eqn{x \ge 0}); puede ser vector.
#' @param lam Parámetro \eqn{\lambda \ge 0}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @param log Si \code{TRUE}, devuelve \code{log(f(x))}.
#' @return Densidad (o log-densidad).
#' @export
dpuig <- function(x, lam, k, T, log = FALSE) {
  f <- Puig_pdf(x, lam, k, T)
  if (log) log(f) else f
}

#' Acumulada (o cola superior) de la distribución de Puig.
#'
#' @param q Puntos de evaluación (\eqn{q \ge 0}); puede ser vector.
#' @param lam Parámetro \eqn{\lambda \ge 0}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @param lower.tail Si \code{TRUE} (por defecto) devuelve \eqn{F(q)};
#'   si \code{FALSE} devuelve la supervivencia \eqn{S(q)}.
#' @return Probabilidad en \eqn{[0, 1]}.
#' @export
ppuig <- function(q, lam, k, T, lower.tail = TRUE) {
  if (lower.tail) Puig_cumulative(q, lam, k, T)
  else Puig_surviving(q, lam, k, T)
}

#' Cuantil de la distribución de Puig.
#'
#' @param p Probabilidades en \eqn{[0, 1]} (escalar o vector).
#' @param lam Parámetro \eqn{\lambda \ge 0}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Cuantil \eqn{Q(p)}.
#' @export
qpuig <- function(p, lam, k, T) {
  Puig_quantile(p, lam, k, T)
}

#' Muestreo aleatorio de la distribución de Puig.
#'
#' @param n Número de muestras.
#' @param lam Parámetro \eqn{\lambda \ge 0}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Vector de longitud \code{n} con \eqn{X \ge 0}.
#' @export
rpuig <- function(n, lam, k, T) {
  Puig_rand(n, lam, k, T)
}