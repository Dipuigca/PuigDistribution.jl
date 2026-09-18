"""Funciones hipergeométricas de alta precisión usadas por los momentos.

La media μ₁ y el tercer momento μ₃ de la distribución de Puig requieren las
funciones de Laguerre generalizadas de orden fraccionario L_{1/2}^{(α)} y
L_{3/2}^{(α)}, que el paquete Julia evalúa con ArbNumerics (128 bits) mediante
la función de Kummer M(a, b, z).

Aquí se replican con mpmath (dps ≈ 40), como referencia de alta precisión.
"""

import mpmath as mp


def kummer_M_stable(a, b, z, maxterms=15000):
    """Serie de Kummer convergente M(a, b, z) con criterio de parada relativo.

    Corresponde a ``kummer_M_stable`` de ``src/puig_pdf.jl`` (ArbNumerics).
    Válida para |z| moderado; para |z| grande usar :func:`kummer_1f1_decay`.
    """
    a = mp.mpf(a)
    b = mp.mpf(b)
    z = mp.mpf(z)
    term = mp.mpf(1)
    total = mp.mpf(1)
    k = 0
    while k < maxterms:
        term *= (a + k) / (b + k) * (z / (k + 1))
        total += term
        if abs(term) < abs(total) * mp.mpf(10) ** (-30):
            break
        k += 1
    return total


def _pochhammer(x, k):
    """Símbolo de Pochhammer (x)ₖ = x(x+1)···(x+k−1)."""
    out = mp.mpf(1)
    for j in range(k):
        out *= x + j
    return out


def kummer_1f1_decay(a, b, y, maxterms=80):
    """Rama decreciente: e^{-y}·M(a, b, y) ≈ Γ(b)/Γ(a)·y^{a-b}·Σₖ (1-a)ₖ(b-a)ₖ/k!·y^{-k}.

    Expansión asintótica estable para y positivo grande (necesario cuando el
    argumento de la serie directa de Kummer es enorme y pierde precisión).
    """
    a = mp.mpf(a)
    b = mp.mpf(b)
    y = mp.mpf(y)
    total = mp.mpf(0)
    fact = mp.mpf(1)
    for k in range(maxterms + 1):
        if k > 0:
            fact *= k
        c = _pochhammer(1 - a, k) * _pochhammer(b - a, k) / fact
        term = c * y ** (-k)
        total += term
        if k >= 6 and abs(term) < abs(total) * mp.mpf(10) ** (-30):
            break
    return mp.gamma(b) / mp.gamma(a) * y ** (a - b) * total


def laguerre_real(n, alpha, x, dps=40):
    """Función de Laguerre generalizada L_n^(α)(x) para n, α reales y x ≤ 0.

    Usa la transformación de Kummer::

        M(-n, α+1, x) = exp(x) · M(α+1+n, α+1, -x)

    con dos ramas numéricamente robustas (nótese y = -x ≥ 0):

    - ``y < 20``: serie directa M(-n, β, -y) (cancelaciones moderadas,
      precisa incluso en doble precisión; aquí con mpmath ≈ exacta).
    - ``y ≥ 20``: expansión asintótica de la rama decreciente
      e^{-y}·M(α+1+n, α+1, y) (la serie directa de Kummer con argumento
      positivo enorme falla incluso en alta precisión).

    Como en ``laguerre_real`` de ``src/puig_pdf.jl``.
    """
    mp.mp.dps = max(dps, 30)
    n = mp.mpf(n)
    alpha = mp.mpf(alpha)
    x = mp.mpf(x)

    coef = mp.gamma(n + alpha + 1) / (
        mp.gamma(n + 1) * mp.gamma(alpha + 1)
    )
    y = -x
    b = alpha + 1
    if y < 20:
        Mval = mp.exp(x) * kummer_M_stable(b + n, b, -x)
    else:
        Mval = kummer_1f1_decay(b + n, b, y)
    return float(coef * Mval)


__all__ = ["kummer_M_stable", "kummer_1f1_decay", "laguerre_real"]