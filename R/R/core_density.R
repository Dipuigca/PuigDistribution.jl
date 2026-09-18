# Núcleo: Marcum-Q y Bessel escalada + densidad.

# e^{-z}·I_nu(z): besselI escalada de R para z < 200 (precisa) y expansión
# asintótica de 5 términos (DLMF 10.40.1) para z >= 200. R devuelve 0 para
# z muy grande con expon.scaled, por lo que la asintótica es obligatoria.
.bessel_ive <- function(nu, z) {
  if (z >= 200.0) besselix_asymp(nu, z) else besselI(z, nu, expon.scaled = TRUE)
}

#' Bessel I modificada escalada asintótica: \eqn{e^{-z} I_\nu(z)}.
#'
#' Expansión asintótica de 5 términos (DLMF 10.40.1). Precisión relativa
#' < 1e-12 para \eqn{z \ge 100} y < 1e-15 para \eqn{z \ge 200}. R devuelve 0
#' para \eqn{z} muy grande con \code{expon.scaled = TRUE}, por lo que la
#' rama asintótica es obligatoria en ese régimen.
#'
#' @param nu Orden real del Bessel (\eqn{\nu}).
#' @param z Argumento (puede ser escalar; se vectoriza con \code{vapply}).
#' @return \code{e^{-z} I_\nu(z)} como vector numérico.
#' @export
besselix_asymp <- function(nu, z) {
  nu <- as.numeric(nu); z <- as.numeric(z)
  mu <- 4.0 * nu * nu
  res <- 1.0
  t1 <- -(mu - 1.0) / (8.0 * z)
  res <- res + t1
  t2 <- -t1 * (mu - 9.0) / (16.0 * z)
  res <- res + t2
  t3 <- -t2 * (mu - 25.0) / (24.0 * z)
  res <- res + t3
  t4 <- -t3 * (mu - 49.0) / (32.0 * z)
  res <- res + t4
  t5 <- -t4 * (mu - 81.0) / (40.0 * z)
  res <- res + t5
  res / sqrt(2.0 * pi * z)
}

#' Marcum-Q asintótica de cola superior lejana.
#'
#' Expansión asintótica de 2 términos (Cantrell 1986) válida para
#' \eqn{b - a \ge 4} y \eqn{ab \ge 30}. Evita cancelaciones catastróficas
#' y falsos suelos de ruido numérico en probabilidades extremas (~1e-50).
#'
#' @param a Parámetro \eqn{a} (p. ej. \eqn{\lambda\sqrt{T}}).
#' @param b Argumento \eqn{b} (p. ej. \eqn{x\sqrt{T}}).
#' @param m Orden de la Marcum-Q (p. ej. \eqn{k/2}).
#' @return Probabilidad de la cola superior en \eqn{[0, 1]}.
#' @export
marcumq_asymp <- function(a, b, m = 1.0) {
  a <- as.numeric(a); b <- as.numeric(b); m <- as.numeric(m)
  if (b - a < 4.0) stop("marcumq_asymp requiere b - a \u2265 4.0")
  if (a * b < 30.0) stop("marcumq_asymp requiere ab \u2265 30.0")

  diff_ab <- b - a
  arg_exp <- -0.5 * diff_ab * diff_ab
  if (arg_exp < -740.0) return(0.0)

  rho <- b / a
  prefactor <- rho^(m - 0.5) * exp(arg_exp) / (sqrt(2.0 * pi) * diff_ab)
  num1 <- b + a + (4.0 * (m - 0.5)^2 - 1.0) / (4.0 * diff_ab)
  denom1 <- 2.0 * a * b * diff_ab
  corr <- 1.0 - num1 / denom1
  val <- prefactor * corr
  if (!is.finite(val)) return(0.0)
  max(val, 0.0)
}

.marcumq_scalar <- function(a, b, m) {
  if (b - a >= 4.0 && a * b >= 30.0) {
    return(marcumq_asymp(a, b, m))
  }
  suppressWarnings(pchisq(b * b, 2.0 * m, ncp = a * a, lower.tail = FALSE))
}

#' Función Marcum-Q generalizada \eqn{Q_m(a,b)}.
#'
#' \deqn{Q_m(a, b) = \operatorname{ccdf}\left(\text{NoncentralChisq}(2m, a^2), b^2\right)}
#' Para la cola superior lejana (\eqn{b - a \ge 4} y \eqn{ab \ge 30}) conmuta
#' automáticamente a \code{\link{marcumq_asymp}}.
#'
#' @param a Parámetro \eqn{a} (p. ej. \eqn{\lambda\sqrt{T}}).
#' @param b Argumento \eqn{b} (p. ej. \eqn{x\sqrt{T}}); puede ser vector.
#' @param m Orden real de la Marcum-Q (p. ej. \eqn{k/2}).
#' @return Probabilidad de la cola superior (escalar o vector como \code{b}).
#' @export
marcumq <- function(a, b, m = 1.0) {
  a <- as.numeric(a); m <- as.numeric(m)
  if (length(b) > 1L) {
    vapply(as.numeric(b), function(bb) .marcumq_scalar(a, bb, m), numeric(1))
  } else {
    .marcumq_scalar(a, as.numeric(b), m)
  }
}

# ---- MB (caso límite λ = 0, B = 1/T) ----

.mb_pdf <- function(x, k, B) {
  k <- as.numeric(k); B <- as.numeric(B)
  if (any(x < 0)) stop("Dominio de la función es [0, \u221e)")
  kh <- k / 2.0
  C <- 2.0^(1.0 - kh) * B^(-kh) / gamma(kh)
  x^(k - 1.0) * exp(-(x^2) / (2.0 * B)) * C
}

.mb_surviving <- function(x, k, B) {
  k <- as.numeric(k); B <- as.numeric(B)
  if (any(x < 0)) stop("Dominio de la función es [0, \u221e)")
  pgamma(x^2 / (2.0 * B), k / 2.0, lower.tail = FALSE)
}

.mb_cumulative <- function(x, k, B) {
  1.0 - .mb_surviving(x, k, B)
}

#' PDF del caso límite Maxwell–Boltzmann (\eqn{\lambda = 0}).
#'
#' \deqn{f(x; k, B) = 2^{1-k/2} B^{-k/2} / \Gamma(k/2) \cdot x^{k-1}
#'   \exp(-x^2 / (2B))}
#'
#' @param x Vector de puntos de evaluación (\eqn{x \ge 0}).
#' @param k Dimensión / grados de libertad reales (\eqn{k \ge 1}).
#' @param B Parámetro de escala (\eqn{B > 0}, con \eqn{B = 1/T}).
#' @return Densidad \eqn{f(x)} (misma longitud que \code{x}).
#' @export
MB_pdf <- function(x, k, B) {
  xout <- .mb_pdf(as.numeric(x), k, B)
  if (length(x) == 1L) as.numeric(xout) else as.numeric(xout)
}

#' Supervivencia del caso límite Maxwell–Boltzmann.
#'
#' \deqn{S(x) = \Gamma(k/2, x^2/(2B)) / \Gamma(k/2)}
#' (gamma incompleta superior regularizada).
#'
#' @param x Vector de puntos de evaluación (\eqn{x \ge 0}).
#' @param k Dimensión / grados de libertad reales (\eqn{k \ge 1}).
#' @param B Parámetro de escala (\eqn{B > 0}).
#' @return Supervivencia \eqn{S(x)} en \eqn{[0, 1]}.
#' @export
MB_surviving <- function(x, k, B) .mb_surviving(as.numeric(x), k, B)

#' Acumulada del caso límite Maxwell–Boltzmann.
#'
#' \deqn{F(x) = 1 - S(x)}
#'
#' @param x Vector de puntos de evaluación (\eqn{x \ge 0}).
#' @param k Dimensión / grados de libertad reales (\eqn{k \ge 1}).
#' @param B Parámetro de escala (\eqn{B > 0}).
#' @return Acumulada \eqn{F(x)} en \eqn{[0, 1]}.
#' @export
MB_cumulative <- function(x, k, B) .mb_cumulative(as.numeric(x), k, B)

# ---- Núcleo de la PDF de Puig ----

.puig_pdf <- function(x, lam, k, T, method = "asymp") {
  lam <- as.numeric(lam); k <- as.numeric(k); T <- as.numeric(T)
  if (lam < 0) stop("\u03bb debe ser \u2265 0")
  if (k < 1) stop("k debe ser \u2265 1")
  if (T <= 0) stop("T debe ser > 0")
  if (any(x < 0)) stop("Dominio de la función es [0, \u221e)")

  if (lam == 0) {
    return(as.numeric(.mb_pdf(x, k, 1.0 / T)))
  }

  nu <- k / 2.0 - 1.0
  half_T <- 0.5 * T
  pow_k_half <- k / 2.0
  lam_denom <- lam^(pow_k_half - 1.0)
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (xi > 0.0) {
      arg_exp <- -half_T * (xi - lam)^2
      if (arg_exp >= -740.0) {
        expo <- exp(arg_exp)
        if (expo > 0.0) {
          z <- xi * lam * T
          bix <- .bessel_ive(nu, z)
          factor <- T * (xi^pow_k_half) / lam_denom
          val <- factor * expo * bix
          out[i] <- if (is.finite(val)) val else 0.0
        }
      }
    }
  }
  out
}

# Referencia "exacta" doble precisión: usa el Bessel exacto de R para
# z < 200 y la expansión asintótica de 5 términos para z >= 200 (R devuelve 0
# con besselI escalada en z muy grande).
.puig_pdf_arb <- function(x, lam, k, T) {
  lam <- as.numeric(lam); k <- as.numeric(k); T <- as.numeric(T)
  if (lam == 0) return(as.numeric(.mb_pdf(x, k, 1.0 / T)))
  nu <- k / 2.0 - 1.0
  half_T <- 0.5 * T
  pow_k_half <- k / 2.0
  lam_denom <- lam^(pow_k_half - 1.0)
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (xi > 0.0) {
      arg_exp <- -half_T * (xi - lam)^2
      if (arg_exp >= -740.0) {
        expo <- exp(arg_exp)
        if (expo > 0.0) {
          z <- xi * lam * T
          bix <- .bessel_ive(nu, z)
          out[i] <- T * (xi^pow_k_half) / lam_denom * expo * bix
        }
      }
    }
  }
  out
}

.puig_logpdf <- function(x, lam, k, T) {
  lam <- as.numeric(lam); k <- as.numeric(k); T <- as.numeric(T)
  if (any(x < 0)) stop("Dominio de la función es [0, \u221e)")
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (xi == 0) {
      out[i] <- -Inf
    } else if (lam == 0) {
      B <- 1.0 / T
      kh <- k / 2.0
      out[i] <- (1.0 - kh) * log(2.0) - kh * log(B) - lgamma(kh) +
        (k - 1.0) * log(xi) - xi^2 / (2.0 * B)
    } else {
      nu <- k / 2.0 - 1.0
      z <- xi * lam * T
      out[i] <- log(T) + (k / 2.0) * log(xi) - (k / 2.0 - 1.0) * log(lam) -
        T / 2.0 * (xi - lam)^2 + log(besselI(z, nu, expon.scaled = TRUE))
    }
  }
  out
}

# ---- API pública (acepta objetos o parámetros escalares) ----

#' Densidad de la distribución de Puig.
#'
#' \deqn{f_P(x; \lambda, k, T) = T x^{k/2} / \lambda^{k/2-1}
#'   \exp(-T/2(x^2+\lambda^2)) I_{k/2-1}(x\lambda T)}
#' con \eqn{\lambda \ge 0}, \eqn{k \ge 1}, \eqn{T > 0}.
#'
#' @param x Puntos de evaluación (\eqn{x \ge 0}); puede ser vector.
#' @param lam Parámetro de ubicación de medias (\eqn{\lambda \ge 0}), objeto
#'   \code{PuigDistribution} o \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}); ignorado si \code{lam}
#'   es un objeto.
#' @param T Escala / precisión (\eqn{T > 0}, \eqn{T = 1/\sigma^2}); ignorado si
#'   \code{lam} es un objeto.
#' @param method \code{"asymp"} (por defecto): Bessel escalada + expansión
#'   asintótica para \eqn{z \ge 200}; \code{"arb"}: referencia doble con Bessel
#'   exacta para \eqn{z < 200}.
#' @return Densidad \eqn{f(x)} (escalar o vector).
#' @export
Puig_pdf <- function(x, lam, k, T, method = "asymp") {
  stopifnot(length(lam) == 1 || (is.PuigDistribution(lam)) || is.MB(lam))
  x <- as.numeric(x)
  if (is.PuigDistribution(lam)) {
    return(.puig_pdf(x, lam$lam, lam$k, lam$T, method = method))
  }
  if (is.MB(lam)) {
    return(.mb_pdf(x, lam$k, 1.0 / lam$B))
  }
  if (method == "arb") {
    return(.puig_pdf_arb(x, lam, k, T))
  }
  .puig_pdf(x, lam, k, T, method = method)
}

#' Log-densidad de la distribución de Puig.
#'
#' \eqn{\log f_P(x)} numéricamente estable mediante Bessel escalada
#' (\code{expon.scaled = TRUE}); \code{-Inf} en \eqn{x = 0}.
#'
#' @param x Puntos de evaluación (\eqn{x \ge 0}); puede ser vector.
#' @param lam Parámetro \eqn{\lambda \ge 0}, objeto \code{PuigDistribution} o
#'   \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Log-densidad (escalar o vector).
#' @export
Puig_logpdf <- function(x, lam, k, T) {
  x <- as.numeric(x)
  if (is.PuigDistribution(lam)) return(.puig_logpdf(x, lam$lam, lam$k, lam$T))
  if (is.MB(lam)) return(.puig_logpdf(x, 0, lam$k, 1.0 / lam$B))
  .puig_logpdf(x, lam, k, T)
}

.puig_surviving <- function(x, lam, k, T) {
  lam <- as.numeric(lam); k <- as.numeric(k); T <- as.numeric(T)
  if (any(x < 0)) stop("Dominio de la función es [0, \u221e)")
  if (lam == 0) {
    return(as.numeric(.mb_surviving(x, k, 1.0 / T)))
  }
  order_m <- k / 2.0
  sT <- sqrt(T)
  a <- lam * sT
  b <- x * sT
  out <- marcumq(a, b, order_m)
  for (i in seq_along(out)) {
    if (is.nan(out[i])) {
      stop(sprintf("marcumq devolvi\u00f3 NaN en x=%g", x[i]))
    }
    if (out[i] <= 0.0 && x[i] < lam) out[i] <- 1.0
  }
  out
}

#' Supervivencia \eqn{S(x) = P(X \ge x)}.
#'
#' \deqn{S(x) = Q_{k/2}(\lambda\sqrt{T}, x\sqrt{T})}
#' mediante la Marcum-Q generalizada (conmuta a la asintótica en colas lejanas).
#'
#' @param x Puntos de evaluación (\eqn{x \ge 0}); puede ser vector.
#' @param lam Parámetro \eqn{\lambda \ge 0}, objeto \code{PuigDistribution} o
#'   \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Supervivencia en \eqn{[0, 1]} (escalar o vector).
#' @export
Puig_surviving <- function(x, lam, k, T) {
  x <- as.numeric(x)
  if (is.PuigDistribution(lam)) return(.puig_surviving(x, lam$lam, lam$k, lam$T))
  if (is.MB(lam)) return(.puig_surviving(x, 0, lam$k, 1.0 / lam$B))
  .puig_surviving(x, lam, k, T)
}

.puig_cumulative <- function(x, lam, k, T) {
  1.0 - .puig_surviving(x, lam, k, T)
}

#' Acumulada \eqn{F(x) = P(X \le x)}.
#'
#' Complementaria de la supervivencia: \eqn{F(x) = 1 - S(x)}.
#'
#' @param x Puntos de evaluación (\eqn{x \ge 0}); puede ser vector.
#' @param lam Parámetro \eqn{\lambda \ge 0}, objeto \code{PuigDistribution} o
#'   \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Acumulada en \eqn{[0, 1]} (escalar o vector).
#' @export
Puig_cumulative <- function(x, lam, k, T) {
  x <- as.numeric(x)
  if (is.PuigDistribution(lam)) return(.puig_cumulative(x, lam$lam, lam$k, lam$T))
  if (is.MB(lam)) return(.puig_cumulative(x, 0, lam$k, 1.0 / lam$B))
  .puig_cumulative(x, lam, k, T)
}

.puig_hazard <- function(x, lam, k, T) {
  f <- .puig_pdf(x, lam, k, T)
  S <- .puig_surviving(x, lam, k, T)
  f / S
}

#' Función de riesgo \eqn{h(x) = f(x)/S(x)}.
#'
#' Tasa de fallo instantánea de la distribución de Puig para \eqn{x} donde
#' \eqn{S(x) > 0}.
#'
#' @param x Puntos de evaluación (\eqn{x \ge 0}); puede ser vector.
#' @param lam Parámetro \eqn{\lambda \ge 0}, objeto \code{PuigDistribution} o
#'   \code{MB}.
#' @param k Dimensión efectiva continua (\eqn{k \ge 1}).
#' @param T Escala / precisión (\eqn{T > 0}).
#' @return Riesgo \eqn{h(x)} (escalar o vector).
#' @export
Puig_hazard <- function(x, lam, k, T) {
  x <- as.numeric(x)
  if (is.PuigDistribution(lam)) return(.puig_hazard(x, lam$lam, lam$k, lam$T))
  if (is.MB(lam)) return(.puig_hazard(x, 0, lam$k, 1.0 / lam$B))
  .puig_hazard(x, lam, k, T)
}