# Wrappers estilo "Distribuciones" (convención R d/p/q/r).

#' Densidad (wrapper de Puig_pdf).
#' @export
dpuig <- function(x, lam, k, T, log = FALSE) {
  f <- Puig_pdf(x, lam, k, T)
  if (log) log(f) else f
}

#' Acumulada con lower.tail (wrapper de Puig_cumulative / Puig_surviving).
#' @export
ppuig <- function(q, lam, k, T, lower.tail = TRUE) {
  if (lower.tail) Puig_cumulative(q, lam, k, T)
  else Puig_surviving(q, lam, k, T)
}

#' Cuantil (wrapper de Puig_quantile).
#' @export
qpuig <- function(p, lam, k, T) {
  Puig_quantile(p, lam, k, T)
}

#' Muestreo aleatorio (wrapper de Puig_rand).
#' @export
rpuig <- function(n, lam, k, T) {
  Puig_rand(n, lam, k, T)
}