# Visualización con ggplot2. Devuelve (una lista de) objetos ggplot.

.puig_xgrid <- function(lam, k, T, npoints = 400, pad_frac = 0.05) {
  ci <- Puig_quantile(c(0.001, 0.999), lam, k, T)
  margen <- pad_frac * (ci[2] - ci[1])
  seq(max(0.0, ci[1] - margen), ci[2], length.out = as.integer(npoints))
}

.puig_plot_dataframe <- function(lam, k, T) {
  st <- Puig_stats(lam, k, T)
  ci98 <- Puig_quantile(c(0.01, 0.99), lam, k, T)
  qs <- Puig_quantile(c(0.25, 0.5, 0.75), lam, k, T)
  sprintf("\u03bb = %.4f   k = %.4f   T = %.4f\n\n
          media    = %.4f\n
          varianza = %.4f\n
          \u03c3        = %.4f\n
          skewness = %.4f\n
          kurtosis = %.4f\n\n
          P1  = %.4f   P99 = %.4f\n\n
          Q1 (25%%) = %.4f\n
          Q2 (50%%) = %.4f\n
          Q3 (75%%) = %.4f",
          lam, k, T,
          st$mean, st$var, st$sig, st$skewness, st$kurtosis,
          ci98[1], ci98[2],
          qs[1], qs[2], qs[3])
}

.puig_ggplot <- function(o, lam, k, T) {
  x <- .puig_xgrid(lam, k, T)
  y <- switch(o,
    pdf = Puig_pdf(x, lam, k, T),
    cdf = Puig_cumulative(x, lam, k, T),
    surv = Puig_surviving(x, lam, k, T),
    data = NULL
  )
  if (o == "data") {
    data.frame(estadisticos = .puig_plot_dataframe(lam, k, T))
  } else {
    dfv <- data.frame(x = x, y = y)
    ggplot2::ggplot(dfv, ggplot2::aes(.data[["x"]], .data[["y"]])) +
      ggplot2::geom_line() +
      ggplot2::labs(title = switch(o,
        pdf = "Densidad (pdf)",
        cdf = "Acumulada (cdf)",
        surv = "Supervivencia"),
        x = "x", y = switch(o, pdf = "f(x)", cdf = "F(x)", surv = "S(x)"))
  }
}

#' Representaciones gráficas de la distribución de Puig (ggplot2).
#'
#' `options` puede ser `"all"` (lista de 4 ggplots), o una de `"pdf"`,
#' `"cdf"`, `"surv"`, `"data"`. Devuelve una lista de objetos ggplot.
#' @export
Puig_plot <- function(dist, options = "all") {
  if (is.MB(dist)) {
    lam <- 0.0; k <- dist$k; T <- 1.0 / dist$B
  } else if (is.PuigDistribution(dist)) {
    lam <- dist$lam; k <- dist$k; T <- dist$T
  } else {
    stop("dist debe ser un objeto PuigDistribution o MB")
  }

  if (is.character(options) && length(options) == 1) {
    if (options != "all" && !options %in% c("pdf", "cdf", "surv", "data")) {
      stop(sprintf('options debe ser "all", "pdf", "cdf", "surv" o "data"; recibido: %s', options))
    }
    if (options == "all") options <- c("pdf", "cdf", "surv", "data")
  }
  if (length(options) == 0) stop("options no puede estar vacío")
  bad <- setdiff(options, c("pdf", "cdf", "surv", "data"))
  if (length(bad) > 0) {
    stop(sprintf("variable de opción desconocida: %s", paste(bad, collapse = ", ")))
  }

  out <- lapply(options, function(o) .puig_ggplot(o, lam, k, T))
  if (length(out) == 1) out[[1]] else out
}