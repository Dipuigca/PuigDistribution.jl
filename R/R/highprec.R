# Momentos μ_1 y μ_3 de la distribución de Puig, que requieren las funciones
# de Laguerre generalizadas de orden fraccionario L_{1/2}^{(α)} y L_{3/2}^{(α)}.
#
# Evaluación numéricamente robusta por ramas (misma estrategia que el núcleo
# Julia corregido y que el port Python):
#   L_n^(α)(x) = Γ(n+α+1)/(Γ(n+1)Γ(α+1)) · M(-n, α+1, x)   (x ≤ 0)
# con y = -x ≥ 0:
#   - y < 20 : serie directa M(-n, b, -y)      (cancelaciones moderadas; doble).
#   - y ≥ 20 : expansión asintótica de la rama decreciente
#              e^{-y}·M(b+n, b, y) ≈ Γ(b)/Γ(a)·y^{a-b}·Σₖ (1-a)ₖ(b-a)ₖ/k!·y^{-k}.
#
# La serie directa de Kummer con argumento positivo enorme (como en la versión
# original de DistributionsPuig.jl) falla incluso en alta precisión.

#' Serie de Kummer convergente M(a, b, z) (doble precisión).
.kummer_M_stable_float <- function(a, b, z, maxterms = 15000) {
  a <- as.numeric(a); b <- as.numeric(b); z <- as.numeric(z)
  term <- 1.0
  total <- 1.0
  k <- 0L
  while (k < maxterms && is.finite(total)) {
    term <- term * (a + k) / (b + k) * (z / (k + 1))
    total <- total + term
    if (abs(term) < abs(total) * 1e-15) break
    k <- k + 1L
  }
  total
}

#' Rama decreciente: e^{-y}·M(a, b, y) por expansión asintótica.
.kummer_1f1_decay <- function(a, b, y, maxterms = 80) {
  a <- as.numeric(a); b <- as.numeric(b); y <- as.numeric(y)
  .poch <- function(x, k) {
    r <- 1.0
    if (k > 0) for (j in seq_len(k)) r <- r * (x + (j - 1))
    r
  }
  total <- 0.0
  fact <- 1.0
  for (k in 0:maxterms) {
    if (k > 0) fact <- fact * k
    c <- .poch(1 - a, k) * .poch(b - a, k) / fact
    term <- c * y^(-k)
    total <- total + term
    if (k >= 6 && abs(term) < abs(total) * 1e-15) break
  }
  exp(lgamma(b) - lgamma(a)) * y^(a - b) * total
}

#' Serie hipergeométrica confluente de Kummer M(a, b, z).
#'
#' Serie directa \eqn{M(a,b,z) = \sum_k (a)_k z^k / ((b)_k k!)} en doble
#' precisión, convergente para \eqn{z} moderado (uso típico: rama
#' \eqn{y < 20} de \code{laguerre_real}).
#'
#' @param a Primer parámetro.
#' @param b Segundo parámetro.
#' @param z Argumento (real).
#' @param maxterms Máximo número de términos antes de abandonar (por defecto
#'   15000).
#' @return Valor de \eqn{M(a,b,z)} en doble precisión.
#' @export
kummer_M_stable <- function(a, b, z, maxterms = 15000) {
  .kummer_M_stable_float(a, b, z, maxterms = maxterms)
}

#' Función de Laguerre generalizada \eqn{L_n^{(\alpha)}(x)} para \eqn{n},
#' \eqn{\alpha} reales y \eqn{x \le 0}.
#'
#' Evaluación numéricamente robusta por ramas (mismo esquema que el núcleo
#' Julia corregido): \eqn{L_n^{(\alpha)}(x) = \Gamma(n+\alpha+1)/(\Gamma(n+1)
#' \Gamma(\alpha+1)) \cdot M(-n, \alpha+1, x)} con \eqn{y = -x \ge 0}: serie
#' directa para \eqn{y < 20} y expansión asintótica de la rama decreciente
#' para \eqn{y \ge 20}, evitando la cancelación catastrófica de la serie de
#' Kummer con argumento positivo enorme (momentos con \eqn{\lambda} grande).
#'
#' @param n Orden real de la Laguerre.
#' @param alpha Parámetro \eqn{\alpha}.
#' @param x Argumento (\eqn{x \le 0}).
#' @param precBits Parámetro de compatibilidad (no utilizado en R; la
#'   evaluación interna es en doble precisión).
#' @return Valor de \eqn{L_n^{(\alpha)}(x)}.
#' @export
laguerre_real <- function(n, alpha, x, precBits = 384) {
  n <- as.numeric(n); alpha <- as.numeric(alpha); x <- as.numeric(x)
  coef <- gamma(n + alpha + 1) / (gamma(n + 1) * gamma(alpha + 1))
  y <- -x
  b <- alpha + 1
  if (y < 20) {
    f <- exp(x) * .kummer_M_stable_float(b + n, b, -x)
  } else {
    f <- .kummer_1f1_decay(b + n, b, y)
  }
  coef * f
}