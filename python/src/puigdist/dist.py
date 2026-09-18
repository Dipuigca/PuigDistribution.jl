"""Integración de la distribución de Puig con ``scipy.stats.rv_continuous``.

Proporciona el generador ``puig``, un subclase de ``rv_continuous`` con
parámetros de forma ``lam, k, T``::

    from puigdist.dist import puig
    dist = puig(3.5, 2.4, 0.25)     # distribución "congelada" estilo scipy
    dist.pdf([...]), dist.cdf([...]), dist.sf([...]), dist.ppf(0.95)
    dist.mean(), dist.var(), dist.std(), dist.skew(), dist.kurtosis()
    dist.rvs(size=1000)

Convención scipy: ``kurtosis()`` devuelve la kurtosis **excesiva**
(γ₂ − 3); la kurtosis no restada γ₂ se obtiene con ``kurtosis() + 3`` o con
``Puig_kurtosis`` de :mod:`puigdist.core`.
"""

import numpy as np
from scipy.stats import rv_continuous

from .core import (
    PuigDistribution as _PuigDistribution,
    Puig_pdf, Puig_logpdf, Puig_cumulative, Puig_surviving,
    Puig_quantile, Puig_mean, Puig_var,
    Puig_skewness, Puig_kurtosis, Puig_entropy,
)
from .core import Puig_rand

__all__ = ["puig_gen", "puig"]


def _scalar_param(p, name):
    """rv_continuous pasa los parámetros de forma como arrays broadcast de x.

    Para una distribución congelada de parámetros escalares, todos los
    elementos son idénticos; se extrae el escalar.
    """
    a = np.asarray(p, dtype=float)
    if a.size == 1:
        return float(a.reshape(-1)[0])
    first = float(a.reshape(-1)[0])
    if np.allclose(a, first):
        return first
    raise ValueError(f"Parámetro '{name}' debe ser escalar en la integración scipy")


class puig_gen(rv_continuous):
    """Generador scipy de la distribución de Puig (λ, k, T)."""

    def _argcheck(self, lam, k, T):
        """Validez de parámetros: λ ≥ 0, k ≥ 1, T > 0 (hook scipy)."""
        return (lam >= 0) & (k >= 1) & (T > 0)

    def _pdf(self, x, lam, k, T):
        """Densidad (hook scipy; delega en Puig_pdf)."""
        return np.atleast_1d(
            Puig_pdf(x, _scalar_param(lam, "lam"), _scalar_param(k, "k"),
                     _scalar_param(T, "T"))
        )

    def _logpdf(self, x, lam, k, T):
        """Log-densidad (hook scipy; delega en Puig_logpdf)."""
        return np.atleast_1d(
            Puig_logpdf(x, _scalar_param(lam, "lam"), _scalar_param(k, "k"),
                        _scalar_param(T, "T"))
        )

    def _cdf(self, x, lam, k, T):
        """Acumulada (hook scipy; delega en Puig_cumulative)."""
        return np.atleast_1d(
            Puig_cumulative(x, _scalar_param(lam, "lam"), _scalar_param(k, "k"),
                            _scalar_param(T, "T"))
        )

    def _sf(self, x, lam, k, T):
        """Supervivencia (hook scipy; delega en Puig_surviving)."""
        return np.atleast_1d(
            Puig_surviving(x, _scalar_param(lam, "lam"), _scalar_param(k, "k"),
                           _scalar_param(T, "T"))
        )

    def _ppf(self, q, lam, k, T):
        """Cuantil (hook scipy; delega en Puig_quantile)."""
        qa = np.asarray(q, float).ravel()
        out = Puig_quantile(qa, _scalar_param(lam, "lam"), _scalar_param(k, "k"),
                            _scalar_param(T, "T"))
        return np.atleast_1d(out)

    def _rvs(self, lam, k, T, size=None, random_state=None):
        """Muestreo aleatorio (hook scipy; delega en Puig_rand)."""
        if size is None:
            size = 1
        return np.asarray(
            Puig_rand(size, _scalar_param(lam, "lam"), _scalar_param(k, "k"),
                      _scalar_param(T, "T"), rng=random_state),
            float,
        )

    def _mean(self, lam, k, T):
        """Media (hook scipy; delega en Puig_mean)."""
        return float(Puig_mean(_scalar_param(lam, "lam"), _scalar_param(k, "k"),
                               _scalar_param(T, "T")))

    def _var(self, lam, k, T):
        """Varianza (hook scipy; delega en Puig_var)."""
        return float(Puig_var(_scalar_param(lam, "lam"), _scalar_param(k, "k"),
                              _scalar_param(T, "T")))

    def _skew(self, lam, k, T):
        """Asimetría (hook scipy; delega en Puig_skewness)."""
        return float(Puig_skewness(_scalar_param(lam, "lam"), _scalar_param(k, "k"),
                                   _scalar_param(T, "T")))

    def _kurtosis(self, lam, k, T):
        """Kurtosis excesiva γ₂ − 3 (hook scipy)."""
        # Convención scipy: kurtosis EXCESIVA (γ₂ − 3)
        return float(Puig_kurtosis(_scalar_param(lam, "lam"), _scalar_param(k, "k"),
                                   _scalar_param(T, "T"))
                     - 3.0)

    def _entropy(self, lam, k, T, *args):
        """Entropía diferencial (hook scipy; delega en Puig_entropy)."""
        d = _PuigDistribution(_scalar_param(lam, "lam"), _scalar_param(k, "k"),
                              _scalar_param(T, "T"))
        return float(Puig_entropy(d))


puig = puig_gen(name="puig", shapes="lam, k, T", a=0.0, b=np.inf)