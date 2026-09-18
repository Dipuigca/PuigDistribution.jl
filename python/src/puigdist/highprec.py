"""Método de referencia de alta precisión (:arb).

Replica ``puig_pdf_arb`` de ``src/puig_pdf.jl``, que evalúa la PDF de la
distribución de Puig con ArbNumerics a 128 bits (~38 dígitos). Aquí se usa
``mpmath`` con ~40 dígitos.
"""

import mpmath as mp


def puig_pdf_arb(x, lam, k, T, dps=40):
    """PDF de alta precisión: T·x^(k/2)/λ^(k/2-1)·e^{-T/2(x²+λ²)}·I_{k/2-1}(xλT)."""
    if x <= 0.0:
        return 0.0
    mp.mp.dps = max(dps, 30)
    nu = mp.mpf(k) / 2 - 1
    lb = mp.mpf(lam)
    Ta = mp.mpf(T)
    ka = mp.mpf(k)
    xf = mp.mpf(x)

    pref = Ta * xf ** (ka / 2) / lb ** (ka / 2 - 1)
    fac = mp.exp(-Ta / 2 * (xf ** 2 + lb ** 2))
    bessel = mp.besseli(nu, xf * lb * Ta)
    return float(pref * fac * bessel)


__all__ = ["puig_pdf_arb"]