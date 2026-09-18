"""Visualización de la distribución de Puig (matplotlib).

Port del módulo ``src/puig_plot.jl``. ``Puig_plot`` acepta las opciones
``"all"`` (panel 2×2), ``"pdf"``, ``"cdf"``, ``"surv"``, ``"data"`` o tuplas
de opciones.
"""

import numpy as np

from .core import PuigDistribution, Puig_pdf, Puig_surviving, Puig_cumulative
from .core import Puig_stats, Puig_quantile, Puig_ci98, Puig_ci996

__all__ = ["Puig_plot"]


def _puig_xgrid(params, npoints=400, pad_frac=0.05):
    P1, P99 = Puig_ci996(params.lam, params.k, params.T)
    margen = pad_frac * (P99 - P1)
    return np.linspace(max(0.0, P1 - margen), P99, int(npoints))


def _plot_pdf(ax, params):
    x = _puig_xgrid(params)
    ax.plot(x, Puig_pdf(x, params.lam, params.k, params.T))
    ax.set_xlabel("x")
    ax.set_ylabel("f(x)")
    ax.set_title("Densidad (pdf)")


def _plot_cdf(ax, params):
    x = _puig_xgrid(params)
    ax.plot(x, Puig_cumulative(x, params.lam, params.k, params.T))
    ax.set_xlabel("x")
    ax.set_ylabel("F(x)")
    ax.set_title("Acumulada (cdf)")


def _plot_surv(ax, params):
    x = _puig_xgrid(params)
    ax.plot(x, Puig_surviving(x, params.lam, params.k, params.T))
    ax.set_xlabel("x")
    ax.set_ylabel("S(x)")
    ax.set_title("Supervivencia")


def _plot_data(ax, params):
    st = Puig_stats(params.lam, params.k, params.T)
    limi, limsu = Puig_ci98(params.lam, params.k, params.T)
    q1, q2, q3 = Puig_quantile([0.25, 0.5, 0.75], params.lam, params.k, params.T)
    texto = (
        f"λ = {params.lam:.4f}   k = {params.k:.4f}   T = {params.T:.4f}\n\n"
        f"media    = {st['mean']:.4f}\n"
        f"varianza = {st['var']:.4f}\n"
        f"σ        = {st['sig']:.4f}\n"
        f"skewness = {st['skewness']:.4f}\n"
        f"kurtosis = {st['kurtosis']:.4f}\n\n"
        f"P1  = {limi:.4f}   P99 = {limsu:.4f}\n\n"
        f"Q1 (25%) = {q1:.4f}\n"
        f"Q2 (50%) = {q2:.4f}\n"
        f"Q3 (75%) = {q3:.4f}"
    )
    ax.axis("off")
    ax.set_title("Parámetros y estadísticos")
    ax.text(0.0, 0.5, texto, va="center", ha="left", fontsize=9)


_PUIG_PLOT_FUNS = {
    "pdf": _plot_pdf,
    "cdf": _plot_cdf,
    "surv": _plot_surv,
    "data": _plot_data,
}


def Puig_plot(dist, options="all"):
    """Genera representaciones gráficas de la distribución de Puig.

    Parámetros
    ----------
    dist : PuigDistribution o MB
    options : str o tupla de str
        ``"all"`` (panel 2×2), ``"pdf"``, ``"cdf"``, ``"surv"``, ``"data"`` o
        una tupla combinando varias de ellas (p. ej. ``("pdf", "surv")``).

    Devuelve un ``matplotlib.figure.Figure``.
    """
    import matplotlib.pyplot as plt

    if hasattr(dist, "B"):  # MB → equivalente con T = 1/B
        dist = PuigDistribution(0.0, float(dist.k), 1.0 / float(dist.B))

    if isinstance(options, str):
        opts = ["pdf", "cdf", "surv", "data"] if options == "all" else [options]
        if options == "all":
            fig, axs = plt.subplots(2, 2, figsize=(9, 7))
            for ax, o in zip(axs.ravel(), opts):
                _PUIG_PLOT_FUNS[o](ax, dist)
            fig.tight_layout()
            return fig
    else:
        opts = list(options)
        if len(opts) == 0:
            raise ValueError("options no puede estar vacío")

    for o in opts:
        if o not in _PUIG_PLOT_FUNS:
            raise ValueError(f'options debe ser "all", "pdf", "cdf", "surv" o "data"; recibido {o!r}')

    n = len(opts)
    ncols = min(n, 3)
    nrows = int(np.ceil(n / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows))
    axs = np.atleast_1d(axs).ravel()
    for ax, o in zip(axs, opts):
        _PUIG_PLOT_FUNS[o](ax, dist)
    for ax in axs[len(opts):]:
        ax.set_axis_off()
    fig.tight_layout()
    return fig