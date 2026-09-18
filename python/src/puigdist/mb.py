"""Caso límite Maxwell-Boltzmann (λ → 0) de la distribución de Puig.

Corresponde a ``src/Maxwell_Boltzmann.jl`` del paquete Julia
``DistributionsPuig.jl``.

La densidad es la de una chi generalizada con dimensión real ``k``::

    f(x; k, B) = 2^{1-k/2} B^{-k/2} / Γ(k/2) · x^{k-1} · exp(-x²/(2B))
    S(x) = Γ(k/2, x²/(2B)) / Γ(k/2)   (gamma incompleta superior, regularizada)
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.special import gammaincc, gammaln, gamma


@dataclass(frozen=True)
class MB:
    """Distribución límite Maxwell-Boltzmann / chi generalizada (λ = 0).

    Parámetros
    ----------
    k : float
        Dimensión / grados de libertad reales (k ≥ 1).
    B : float
        Parámetro de escala o dispersión térmica (B > 0, con B = 1/T).
    """

    k: float
    B: float

    def __post_init__(self):
        """Valida los parámetros (k ≥ 1, B > 0)."""
        if not (self.k >= 1):
            raise ValueError("k debe ser ≥ 1")
        if not (self.B > 0):
            raise ValueError("B debe ser > 0")

    def pdf(self, x):
        """Densidad f(x; k, B)."""
        return MB_pdf(x, self.k, self.B)

    def survival(self, x):
        """Supervivencia S(x) = P(X ≥ x)."""
        return MB_surviving(x, self.k, self.B)

    def cumulative(self, x):
        """Acumulada F(x) = P(X ≤ x)."""
        return MB_cumulative(x, self.k, self.B)

    def mean(self):
        """Media teórica E[X] (caso límite λ = 0)."""
        from .core import Puig_mean

        return Puig_mean(0.0, self.k, 1.0 / self.B)

    def var(self):
        """Varianza teórica Var(X) (caso límite λ = 0)."""
        from .core import Puig_var

        return Puig_var(0.0, self.k, 1.0 / self.B)

    def std(self):
        """Desviación estándar teórica σ (caso límite λ = 0)."""
        from .core import Puig_std

        return Puig_std(0.0, self.k, 1.0 / self.B)

    def moments(self, n=4):
        """Momentos raw μ₁, …, μₙ (caso límite λ = 0)."""
        from .core import Puig_moments

        return Puig_moments(0.0, self.k, 1.0 / self.B, n)


def _check(k, B):
    """Valida y convierte (k ≥ 1, B > 0) a float."""
    if not (B > 0):
        raise ValueError("B debe ser > 0")
    if not (k >= 1):
        raise ValueError("k debe ser ≥ 1")
    return float(k), float(B)


def MB_pdf(x, k, B):
    """PDF del caso límite Maxwell-Boltzmann, vectorizada sobre ``x``."""
    k, B = _check(k, B)
    xa = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any(xa < 0):
        raise ValueError("Dominio de la función es [0, ∞)")
    kh = k / 2.0
    E = np.exp(-(xa ** 2) / (2.0 * B))
    G = gamma(kh)
    C = 2.0 ** (1.0 - kh) * B ** (-kh) / G
    xt = xa ** (k - 1.0)
    out = xt * E * C
    return float(out[0]) if np.ndim(x) == 0 else out


def MB_surviving(x, k, B):
    """Supervivencia S(x) = Γ(k/2, x²/(2B)) / Γ(k/2), vectorizada."""
    k, B = _check(k, B)
    xa = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any(xa < 0):
        raise ValueError("Dominio de la función es [0, ∞)")
    kh = k / 2.0
    P = xa ** 2 / (2.0 * B)
    out = gammaincc(kh, P)
    return float(out[0]) if np.ndim(x) == 0 else out


def MB_cumulative(x, k, B):
    """Acumulada F(x) = 1 - S(x), vectorizada."""
    out = 1.0 - MB_surviving(x, k, B)
    return out


__all__ = ["MB", "MB_pdf", "MB_surviving", "MB_cumulative"]