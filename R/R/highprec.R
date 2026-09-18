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

#' Función de Laguerre generalizada L_n^(α)(x) para n, α reales y x ≤ 0.
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