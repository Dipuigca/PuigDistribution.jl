# ---- Constructors ----

new_puig <- function(lam, k, T) {
  structure(list(lam = as.numeric(lam), k = as.numeric(k), T = as.numeric(T)),
            class = "PuigDistribution")
}

#' Distribución de Puig: chi no central generalizada de dimensión real k.
#'
#' @param lam Norma del vector de medias (\eqn{\lambda \ge 0}).
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Parámetro de escala / precisión (\eqn{T > 0}, \eqn{T = 1/\sigma^2}).
#' @return Objeto de clase \code{PuigDistribution}: lista con campos
#'   \code{lam}, \code{k} y \code{T}.
#' @export
PuigDistribution <- function(lam, k, T) {
  lam <- as.numeric(lam); k <- as.numeric(k); T <- as.numeric(T)
  if (lam < 0) stop("\u03bb debe ser \u2265 0")
  if (k < 1) stop("k debe ser \u2265 1")
  if (T <= 0) stop("T debe ser > 0")
  new_puig(lam, k, T)
}

print.PuigDistribution <- function(x, ...) {
  cat(sprintf("PuigDistribution(\u03bb=%g, k=%g, T=%g)\n", x$lam, x$k, x$T))
  invisible(x)
}

new_mb <- function(k, B) {
  structure(list(k = as.numeric(k), B = as.numeric(B)), class = "MB")
}

#' Distribución límite Maxwell–Boltzmann (\eqn{\lambda = 0}).
#'
#' Caso límite de la distribución de Puig con \eqn{\lambda = 0}; parametrizada
#' por \eqn{k} (dimensión) y \eqn{B = 1/T} (escala / dispersión térmica).
#'
#' @param k Dimensión / grados de libertad reales (\eqn{k \ge 1}).
#' @param B Parámetro de escala (\eqn{B > 0}).
#' @return Objeto de clase \code{MB}: lista con campos \code{k} y \code{B}.
#' @export
MB <- function(k, B) {
  k <- as.numeric(k); B <- as.numeric(B)
  if (k < 1) stop("k debe ser \u2265 1")
  if (B <= 0) stop("B debe ser > 0")
  new_mb(k, B)
}

print.MB <- function(x, ...) {
  cat(sprintf("MB(k=%g, B=%g)\n", x$k, x$B))
  invisible(x)
}

is.PuigDistribution <- function(x) inherits(x, "PuigDistribution")
is.MB <- function(x) inherits(x, "MB")


# ---- Fitted result / stats containers ----

new_PuigStats <- function(R2, MAE, RMSE, MaxAE, IAE) {
  structure(list(R2 = R2, MAE = MAE, RMSE = RMSE, MaxAE = MaxAE, IAE = IAE),
            class = "PuigStats")
}

new_PuigFitResult <- function(params, indicator, stats) {
  structure(list(params = params, indicator = indicator, stats = stats),
            class = "PuigFitResult")
}

print.PuigStats <- function(x, ...) {
  cat(sprintf("PuigStats(R\u00b2=%.5f, MAE=%.5f, RMSE=%.5f, MaxAE=%.5f, IAE=%.5f)\n",
              x$R2, x$MAE, x$RMSE, x$MaxAE, x$IAE))
  invisible(x)
}

print.PuigFitResult <- function(x, ...) {
  cat(sprintf("PuigFitResult(\u03bb=%.4f, k=%.4f, T=%.5f)\n",
              x$params$lam, x$params$k, x$params$T))
  cat(sprintf("  indicador: %d obs\n", length(x$indicator)))
  print(x$stats)
  invisible(x)
}

#' Desempaquetado tipo tupla: p, ind, st = Puig_fit(...)
#' @export
`[[.PuigFitResult` <- function(x, i) {
  unclass(x)[[i]]
}
