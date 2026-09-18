"""puigdist — Distribución de Puig.

Port a Python del paquete Julia ``DistributionsPuig.jl``: distribución chi no
central generalizada con dimensión real continua k (k ≥ 1).

Módulos:
- :mod:`puigdist.core`    núcleo matemático (pdf, logpdf, surv/cdf, momentos,
                          cuantiles, muestreo, Marcum-Q) + tipo ``PuigDistribution``.
- :mod:`puigdist.mb`      caso límite Maxwell-Boltzmann (λ = 0).
- :mod:`puigdist.highprec` método de referencia :arb (mpmath).
- :mod:`puigdist.dist`    integración con ``scipy.stats.rv_continuous``.
- :mod:`puigdist.fit`     ajuste de datos (DataFrame/Matrix/Vector).
- :mod:`puigdist.plot`    visualización (matplotlib).
"""

from .core import (
    PuigDistribution,
    Puig_pdf, Puig_logpdf, Puig_surviving, Puig_cumulative,
    Puig_mean, Puig_var, Puig_std, Puig_skewness, Puig_kurtosis,
    Puig_moments, Puig_moments_raw, Puig_stats, Puig_hazard,
    Puig_quantile, Puig_ci98, Puig_ci996, Puig_entropy, Puig_rand,
    marcumq, marcumq_asymp, besselix_asymp,
)
from .mb import MB, MB_pdf, MB_surviving, MB_cumulative
from .fit import (
    Puig_fit, PuigStats, PuigFitResult,
    effective_dimension, empirical_survival,
    calc_r2, calc_mae, calc_rmse, calc_maxae, calc_iae, calc_all_metrics,
)
from .plot import Puig_plot

__version__ = "0.1.0"

__all__ = [
    "PuigDistribution",
    "Puig_pdf", "Puig_logpdf", "Puig_surviving", "Puig_cumulative",
    "Puig_mean", "Puig_var", "Puig_std", "Puig_skewness", "Puig_kurtosis",
    "Puig_moments", "Puig_moments_raw", "Puig_stats", "Puig_hazard",
    "Puig_quantile", "Puig_ci98", "Puig_ci996", "Puig_entropy", "Puig_rand",
    "marcumq", "marcumq_asymp", "besselix_asymp",
    "MB", "MB_pdf", "MB_surviving", "MB_cumulative",
    "Puig_fit", "PuigStats", "PuigFitResult",
    "effective_dimension", "empirical_survival",
    "calc_r2", "calc_mae", "calc_rmse", "calc_maxae", "calc_iae",
    "calc_all_metrics",
    "Puig_plot",
]