# Núcleo: Marcum-Q y Bessel escalada + densidad.

# e^{-z}·I_nu(z): besselI escalada de R para z < 200 (precisa) y expansión
# asintótica de 5 términos (DLMF 10.40.1) para z >= 200. R devuelve 0 para
# z muy grande con expon.scaled, por lo que la asintótica es obligatoria.
.bessel_ive <- function(nu, z) {
  if (z >= 200.0) besselix_asymp(nu, z) else besselI(z, nu, expon.scaled = TRUE)
}

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

#' Función Marcum-Q generalizada Q_m(a,b). Vectorizada sobre b.
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

MB_pdf <- function(x, k, B) {
  xout <- .mb_pdf(as.numeric(x), k, B)
  if (length(x) == 1L) as.numeric(xout) else as.numeric(xout)
}

MB_surviving <- function(x, k, B) .mb_surviving(as.numeric(x), k, B)
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

#' Supervivencia S(x) = P(X \u2265 x).
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

#' Acumulada F(x) = P(X \u2264 x).
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

#' Función de riesgo h(x) = f(x)/S(x).
#' @export
Puig_hazard <- function(x, lam, k, T) {
  x <- as.numeric(x)
  if (is.PuigDistribution(lam)) return(.puig_hazard(x, lam$lam, lam$k, lam$T))
  if (is.MB(lam)) return(.puig_hazard(x, 0, lam$k, 1.0 / lam$B))
  .puig_hazard(x, lam, k, T)
}