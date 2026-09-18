"""Núcleo matemático de la distribución de Puig.

Implementación en Python (NumPy/SciPy) del paquete Julia
``DistributionsPuig.jl`` · núcleo analítico ``src/puig_pdf.jl``.

Definición (distribución chi no central generalizada, dimensión real k)::

    f_P(x; λ, k, T) = T·x^{k/2}/λ^{k/2-1}·exp(-T/2·(x²+λ²))·I_{k/2-1}(x·λ·T)

con λ ≥ 0 (norma del vector de medias), k ≥ 1 (dimensión efectiva continua)
y T > 0 (precisión / escala inversa, T = 1/σ²).

La supervivencia / acumulada se calculan mediante la función Marcum-Q
generalizada de orden m = k/2:

    S(x) = Q_{k/2}(λ√T, x√T) = ccdf(NoncentralChisq(2m, a²), b²)

con conmutación automática a la expansión asintótica ``marcumq_asymp`` en la
cola superior lejana (b - a ≥ 4 y ab ≥ 30).
"""

import numpy as np
from scipy.special import ive
from scipy.stats import ncx2

from .hypergeom import laguerre_real
from .highprec import puig_pdf_arb
from .mb import MB, MB_pdf, MB_surviving, MB_cumulative

__all__ = [
    "PuigDistribution",
    "Puig_pdf", "Puig_logpdf", "Puig_surviving", "Puig_cumulative",
    "Puig_mean", "Puig_var", "Puig_std", "Puig_skewness", "Puig_kurtosis",
    "Puig_moments", "Puig_moments_raw", "Puig_stats", "Puig_hazard",
    "Puig_quantile", "Puig_ci98", "Puig_ci996", "Puig_entropy", "Puig_rand",
    "marcumq", "marcumq_asymp", "besselix_asymp",
    "MB", "MB_pdf", "MB_surviving", "MB_cumulative",
]


# ---------------------------------------------------------------------------
# Funciones auxiliares: Marcum-Q y Bessel escalada asintótica
# ---------------------------------------------------------------------------

def besselix_asymp(nu, z):
    """e^{-z}·I_ν(z) mediante expansión asintótica de 5 términos (DLMF 10.40.1).

    Precisión relativa < 1e-12 para z ≥ 100 y < 1e-15 para z ≥ 200.
    Equivale a ``besselix_asymp`` de ``src/puig_pdf.jl``.
    """
    mu = 4.0 * nu * nu
    res = 1.0
    t1 = -(mu - 1.0) / (8.0 * z)
    res += t1
    t2 = -t1 * (mu - 9.0) / (16.0 * z)
    res += t2
    t3 = -t2 * (mu - 25.0) / (24.0 * z)
    res += t3
    t4 = -t3 * (mu - 49.0) / (32.0 * z)
    res += t4
    t5 = -t4 * (mu - 81.0) / (40.0 * z)
    res += t5
    return res / np.sqrt(2.0 * np.pi * z)


def marcumq_asymp(a, b, m=1.0):
    """Q_m(a, b) por expansión asintótica de 2 términos (Cantrell 1986).

    Requiere b - a ≥ 4.0 y ab ≥ 30.0 (cola superior lejana). Evita
    cancelaciones catastróficas y falsos suelos de ruido en probabilidades
    extremas (~1e-50). Equivale a ``marcumq_asymp`` de Julia.
    """
    a = float(a)
    b = float(b)
    m = float(m)
    if (b - a) < 4.0:
        raise ValueError("marcumq_asymp requiere b - a ≥ 4.0")
    if a * b < 30.0:
        raise ValueError("marcumq_asymp requiere ab ≥ 30.0")

    diff = b - a
    arg_exp = -0.5 * diff * diff
    if arg_exp < -740.0:
        return 0.0

    rho = b / a
    prefactor = rho ** (m - 0.5) * np.exp(arg_exp) / (np.sqrt(2.0 * np.pi) * diff)
    num1 = b + a + (4.0 * (m - 0.5) ** 2 - 1.0) / (4.0 * diff)
    denom1 = 2.0 * a * b * diff
    corr = 1.0 - num1 / denom1
    val = prefactor * corr
    return float(max(val, 0.0)) if np.isfinite(val) else 0.0


def _marcumq_scalar(a, b, m):
    """Marcum-Q escalar: asintótica en cola lejana o ncx2.sf en otro caso."""
    if (b - a) >= 4.0 and a * b >= 30.0:
        return marcumq_asymp(a, b, m)
    return float(ncx2.sf(b * b, 2.0 * m, a * a))


def marcumq(a, b, m=1.0):
    """Función Marcum-Q generalizada Q_m(a, b) = ccdf(NoncentralChisq(2m, a²), b²).

    - ``b`` puede ser escalar o array.
    - Para la cola superior lejana (b - a ≥ 4 y ab ≥ 30) conmuta a
      ``marcumq_asymp``.

    Equivale a ``marcumq`` de ``src/puig_pdf.jl``.
    """
    a = float(a)
    m = float(m)
    b_arr = np.asarray(b, dtype=float)
    scalar = b_arr.ndim == 0
    out = np.zeros_like(b_arr, dtype=float)
    flat = out.ravel()
    bflat = b_arr.ravel()
    for i in range(bflat.size):
        flat[i] = _marcumq_scalar(a, bflat[i], m)
    return float(out) if scalar else out


# ---------------------------------------------------------------------------
# Tipo principal
# ---------------------------------------------------------------------------

class PuigDistribution:
    """Distribución de Puig: chi no central generalizada de dimensión real k.

    Parámetros
    ----------
    lam : float
        Norma del vector de medias de las componentes (λ ≥ 0).
    k : float
        Dimensión efectiva continua (k ≥ 1).
    T : float
        Parámetro de escala / precisión (T > 0, T = 1/σ²).

    Compatible con las funciones ``Puig_*`` del módulo; la variante integrada
    con ``scipy.stats.rv_continuous`` está en :mod:`puigdist.dist`.
    """

    def __init__(self, lam, k, T):
        """Construye una distribución de Puig.

        Parámetros
        ----------
        lam : float
            Norma del vector de medias (λ ≥ 0).
        k : float
            Dimensión efectiva continua (k ≥ 1).
        T : float
            Precisión / escala inversa (T > 0, T = 1/σ²).

        Lanza ``ValueError`` si los parámetros quedan fuera de su dominio.
        """
        lam, k, T = float(lam), float(k), float(T)
        if not lam >= 0:
            raise ValueError("λ debe ser ≥ 0")
        if not k >= 1:
            raise ValueError("k debe ser ≥ 1")
        if not T > 0:
            raise ValueError("T debe ser > 0")
        self.lam = lam
        self.k = k
        self.T = T

    def __iter__(self):
        """Desempaquetado idiomático: ``lam, k, T = dist``."""
        yield self.lam
        yield self.k
        yield self.T

    def __repr__(self):
        """Representación compacta de la distribución."""
        return (f"PuigDistribution(λ={self.lam:g}, k={self.k:g}, "
                f"T={self.T:g})")

    @property
    def lambda_(self):
        """Alias del parámetro λ (evita colisión con la palabra clave ``lambda``)."""
        return self.lam

    def pdf(self, x, method="asymp"):
        """Densidad f_P(x).

        ``method='asymp'`` (por defecto) o ``'arb'`` (referencia mpmath).
        """
        return Puig_pdf(x, self.lam, self.k, self.T, method=method)

    def logpdf(self, x):
        """Log-densidad log f_P(x), estable numéricamente."""
        return Puig_logpdf(x, self.lam, self.k, self.T)

    def survival(self, x):
        """Función de supervivencia S(x) = P(X ≥ x)."""
        return Puig_surviving(x, self.lam, self.k, self.T)

    def cumulative(self, x):
        """Función de distribución F(x) = P(X ≤ x)."""
        return Puig_cumulative(x, self.lam, self.k, self.T)

    def quantile(self, p, tol=1e-10, maxit=200):
        """Cuantil Q(p) tal que F(Q(p)) = p, por bisección."""
        return Puig_quantile(p, self.lam, self.k, self.T, tol=tol, maxit=maxit)

    def mean(self):
        """Media teórica E[X]."""
        return Puig_mean(self.lam, self.k, self.T)

    def var(self):
        """Varianza teórica Var(X)."""
        return Puig_var(self.lam, self.k, self.T)

    def std(self):
        """Desviación estándar teórica σ."""
        return Puig_std(self.lam, self.k, self.T)

    def skewness(self):
        """Asimetría estandarizada γ₁."""
        return Puig_skewness(self.lam, self.k, self.T)

    def kurtosis(self):
        """Kurtosis (4º momento estandarizado, no restado) γ₂."""
        return Puig_kurtosis(self.lam, self.k, self.T)

    def moments(self, n=4):
        """Momentos raw μ₁, …, μₙ (los 4 primeros como tupla con n=4)."""
        return Puig_moments(self.lam, self.k, self.T, n)

    def stats(self):
        """Diccionario con mean, var, sig, skewness, kurtosis."""
        return Puig_stats(self.lam, self.k, self.T)

    def hazard(self, x):
        """Función de riesgo h(x) = f(x) / S(x)."""
        return Puig_hazard(x, self.lam, self.k, self.T)

    def rvs(self, size=1, rng=None):
        """Muestreo pseudoaleatorio exacto de tamaño ``size``.

        ``rng`` puede ser un ``numpy.random.Generator``/``RandomState``.
        """
        return Puig_rand(size, self.lam, self.k, self.T, rng=rng)

    def entropy(self, N=400):
        """Entropía diferencial H(X) = -∫ f(x)·log f(x) dx (trapecio)."""
        return Puig_entropy(self, N=N)


# ---------------------------------------------------------------------------
# PDF / Log-PDF
# ---------------------------------------------------------------------------

def _check_params(lam, k, T):
    """Valida y convierte los parámetros (λ≥0, k≥1, T>0)."""
    if not lam >= 0:
        raise ValueError("λ debe ser ≥ 0")
    if not k >= 1:
        raise ValueError("k debe ser ≥ 1")
    if not T > 0:
        raise ValueError("T debe ser > 0")
    return float(lam), float(k), float(T)


def Puig_pdf(x, lam, k, T, method="asymp"):
    """Función de densidad de probabilidad f_P(x).

    - ``method='asymp'`` (por defecto): besselix (escalada) para xλT < 200 y
      expansión asintótica de 5 términos para xλT ≥ 200 (Float64 puro).
    - ``method='arb'``: referencia de alta precisión (mpmath, ~40 dígitos).
    """
    lam, k, T = _check_params(lam, k, T)
    xa = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any(xa < 0):
        raise ValueError("Dominio de la función es [0, ∞)")
    scalar = np.ndim(x) == 0

    if lam == 0:
        out = MB_pdf(xa, k, 1.0 / T)
        return float(out[0]) if scalar else out

    if method == "arb":
        out = np.array([puig_pdf_arb(xi, lam, k, T) for xi in xa])
        return float(out[0]) if scalar else out

    if method != "asymp":
        raise ValueError("method debe ser 'asymp' o 'arb'")

    nu = k / 2.0 - 1.0
    half_T = 0.5 * T
    pow_k_half = k / 2.0
    lam_denom = lam ** (pow_k_half - 1.0)
    out = np.zeros_like(xa)
    for i, xi in enumerate(xa):
        if xi > 0.0:
            arg_exp = -half_T * (xi - lam) ** 2
            if arg_exp >= -740.0:
                expo = np.exp(arg_exp)
                if expo > 0.0:
                    z = xi * lam * T
                    bix = besselix_asymp(nu, z) if z >= 200.0 else ive(nu, z)
                    factor = T * (xi ** pow_k_half) / lam_denom
                    val = factor * expo * bix
                    out[i] = val if np.isfinite(val) else 0.0
    return float(out[0]) if scalar else out


def Puig_logpdf(x, lam, k, T):
    """Log-densidad log f_P(x) numéricamente estable (besselix escalada)."""
    lam, k, T = _check_params(lam, k, T)
    xa = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any(xa < 0):
        raise ValueError("Dominio de la función es [0, ∞)")
    scalar = np.ndim(x) == 0

    out = np.empty_like(xa)
    for i, xi in enumerate(xa):
        if xi == 0:
            out[i] = -np.inf
        elif lam == 0:
            B = 1.0 / T
            kh = k / 2.0
            out[i] = ((1.0 - kh) * np.log(2.0) - kh * np.log(B)
                      - _loggamma(kh) + (k - 1.0) * np.log(xi)
                      - xi ** 2 / (2.0 * B))
        else:
            nu = k / 2.0 - 1.0
            z = xi * lam * T
            out[i] = (np.log(T) + (k / 2.0) * np.log(xi)
                      - (k / 2.0 - 1.0) * np.log(lam)
                      - T / 2.0 * (xi - lam) ** 2 + np.log(ive(nu, z)))
    return float(out[0]) if scalar else out


def _loggamma(x):
    """Log-gamma (wrapper de scipy.special.gammaln)."""
    from scipy.special import gammaln

    return gammaln(x)


# ---------------------------------------------------------------------------
# Supervivencia / Acumulada (Marcum-Q)
# ---------------------------------------------------------------------------

def Puig_surviving(x, lam, k, T):
    """Supervivencia S(x) = P(X ≥ x) = Q_{k/2}(λ√T, x√T)."""
    lam, k, T = _check_params(lam, k, T)
    xa = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any(xa < 0):
        raise ValueError("Dominio de la función es [0, ∞)")
    scalar = np.ndim(x) == 0

    if lam == 0:
        out = MB_surviving(xa, k, 1.0 / T)
    else:
        order = k / 2.0
        sT = np.sqrt(T)
        a = lam * sT
        b = xa * sT
        out = marcumq(a, b, order)
        # Réplica del ajuste de robustez de Julia
        for i in range(len(out)):
            if np.isnan(out[i]):
                raise ValueError(f"marcumq devolvió NaN en x={xa[i]}, a={a}, m={order}")
            if out[i] <= 0.0 and xa[i] < lam:
                out[i] = 1.0
    return float(out[0]) if scalar else out


def Puig_cumulative(x, lam, k, T):
    """Acumulada F(x) = P(X ≤ x) = 1 - S(x)."""
    lam, k, T = _check_params(lam, k, T)
    xa = np.atleast_1d(np.asarray(x, dtype=float))
    if np.any(xa < 0):
        raise ValueError("Dominio de la función es [0, ∞)")
    scalar = np.ndim(x) == 0
    out = np.array(1.0 - np.atleast_1d(Puig_surviving(xa, lam, k, T)))
    return float(out[0]) if scalar else out


# ---------------------------------------------------------------------------
# Momentos y estadísticos teóricos
# ---------------------------------------------------------------------------

def _mean_formula(lam, k, T):
    """Fórmula de μ₁ = √(π/(2T))·L_{1/2}^{(k/2-1)}(-λ²T/2)."""
    ls = lam * np.sqrt(T)
    return np.sqrt(np.pi / (2.0 * T)) * laguerre_real(0.5, k / 2.0 - 1.0, -ls ** 2 / 2.0)


def _moment3_formula(lam, k, T):
    """Fórmula de μ₃ = (3/T)·√(π/(2T))·L_{3/2}^{(k/2-1)}(-λ²T/2)."""
    ls = lam * np.sqrt(T)
    L = laguerre_real(1.5, k / 2.0 - 1.0, -ls ** 2 / 2.0)
    return (3.0 / T) * np.sqrt(np.pi / (2.0 * T)) * L


def Puig_mean(lam, k, T):
    """Media teórica E[X] = √(π/(2T))·L_{1/2}^{(k/2-1)}(-λ²T/2)."""
    lam, k, T = _check_params(lam, k, T)
    return _mean_formula(lam, k, T)


def Puig_var(lam, k, T):
    """Varianza teórica E[X²] − (E[X])² con E[X²] = k/T + λ²."""
    lam, k, T = _check_params(lam, k, T)
    mu2 = k / T + lam ** 2
    return mu2 - Puig_mean(lam, k, T) ** 2


def Puig_std(lam, k, T):
    """Desviación estándar σ = √Var(X)."""
    lam, k, T = _check_params(lam, k, T)
    return np.sqrt(Puig_var(lam, k, T))


def Puig_skewness(lam, k, T):
    """Asimetría estandarizada γ₁ = E[(X-μ)³]/σ³."""
    lam, k, T = _check_params(lam, k, T)
    mu1 = Puig_mean(lam, k, T)
    mu2 = Puig_var(lam, k, T) + mu1 ** 2
    mu3 = _moment3_formula(lam, k, T)
    sigma3 = Puig_std(lam, k, T) ** 3
    return (mu3 - 3.0 * mu1 * mu2 + 2.0 * mu1 ** 3) / sigma3


def Puig_kurtosis(lam, k, T):
    """Kurtosis (4º momento estandarizado no restado) γ₂ = E[(X-μ)⁴]/σ⁴."""
    lam, k, T = _check_params(lam, k, T)
    mu1 = Puig_mean(lam, k, T)
    mu2 = k / T + lam ** 2
    mu4 = mu2 ** 2 + 2.0 * k / T ** 2 + 4.0 * lam ** 2 / T
    mu3 = _moment3_formula(lam, k, T)
    sigma4 = Puig_var(lam, k, T) ** 2
    return (mu4 - 4.0 * mu1 * mu3 + 6.0 * mu1 ** 2 * mu2 - 3.0 * mu1 ** 4) / sigma4


def Puig_moments_raw(lam, k, T, n):
    """Los n primeros momentos raw μ₁,…,μₙ (recurrencia de tres términos O(n)).

    Semillas: μ₁ (Laguerre), μ₂ = k/T + λ², μ₃ (Laguerre), μ₄ = μ₂² + 2k/T² + 4λ²/T.

    Recurrencia (n ≥ 5)::

        μₙ = ((2n-4+k)/T + λ²)·μₙ₋₂ − (n-2)(n+k-4)/T²·μₙ₋₄
    """
    lam, k, T = _check_params(lam, k, T)
    n = int(n)
    if n < 1:
        raise ValueError(f"n debe ser ≥ 1, recibido n={n}")

    mu = np.zeros(n)
    mu[0] = Puig_mean(lam, k, T)
    if n == 1:
        return mu
    mu[1] = k / T + lam ** 2
    if n == 2:
        return mu
    mu[2] = _moment3_formula(lam, k, T)
    if n == 3:
        return mu
    mu[3] = mu[1] ** 2 + 2.0 * k / T ** 2 + 4.0 * lam ** 2 / T
    if n == 4:
        return mu

    invT = 1.0 / T
    invT2 = invT * invT
    lam2 = lam * lam
    for idx in range(4, n):
        order = idx + 1
        mu[idx] = (((2.0 * order - 4.0 + k) * invT + lam2) * mu[idx - 2]
                   - (order - 2.0) * (order + k - 4.0) * invT2 * mu[idx - 4])
    return mu


def Puig_moments(lam, k, T, n=4):
    """Momentos raw. Con ``n=4`` devuelve tupla (μ₁,μ₂,μ₃,μ₄); con ``n``
    genérico devuelve ``np.ndarray`` [μ₁,…,μₙ]."""
    lam, k, T = _check_params(lam, k, T)
    mu = Puig_moments_raw(lam, k, T, int(n))
    if int(n) == 4:
        return tuple(mu)
    return mu


def Puig_stats(lam, k, T):
    """NamedTuple-like con mean, var, sig, skewness, kurtosis."""
    lam, k, T = _check_params(lam, k, T)
    mu = Puig_moments_raw(lam, k, T, 4)
    mu1, mu2, mu3, mu4 = mu
    sig2 = mu2 - mu1 ** 2
    sig = np.sqrt(sig2)
    sk = (mu3 - 3.0 * mu1 * mu2 + 2.0 * mu1 ** 3) / sig ** 3
    kur = (mu4 - 4.0 * mu1 * mu3 + 6.0 * mu1 ** 2 * mu2 - 3.0 * mu1 ** 4) / sig2 ** 2
    return {"mean": mu1, "var": sig2, "sig": sig, "skewness": sk, "kurtosis": kur}


def Puig_hazard(x, lam, k, T):
    """Función de riesgo h(x) = f(x) / S(x)."""
    f = np.atleast_1d(Puig_pdf(x, lam, k, T))
    S = np.atleast_1d(Puig_surviving(x, lam, k, T))
    with np.errstate(divide="ignore", invalid="ignore"):
        out = f / S
    return float(out[0]) if np.ndim(x) == 0 else out


# ---------------------------------------------------------------------------
# Cuantiles e intervalos
# ---------------------------------------------------------------------------

def Puig_quantile(p, lam, k, T, tol=1e-10, maxit=200):
    """Cuantil Q(p) tal que F(Q(p)) = p, por bisección."""
    lam, k, T = _check_params(lam, k, T)
    p_arr = np.atleast_1d(np.asarray(p, dtype=float))
    if np.any((p_arr < 0) | (p_arr > 1)):
        raise ValueError("p debe estar en [0, 1]")

    out = np.empty_like(p_arr)
    for i, pi in enumerate(p_arr):
        if pi == 0:
            out[i] = 0.0
        elif pi == 1:
            out[i] = np.inf
        else:
            lo, hi = 0.0, max(lam, 1.0) * 2.0
            while Puig_cumulative([hi], lam, k, T)[0] < pi:
                hi *= 2.0
            for _ in range(int(maxit)):
                mid = 0.5 * (lo + hi)
                if Puig_cumulative([mid], lam, k, T)[0] < pi:
                    lo = mid
                else:
                    hi = mid
                if (hi - lo) < tol:
                    break
            out[i] = 0.5 * (lo + hi)
    return float(out[0]) if np.ndim(p) == 0 else out


def Puig_ci98(lam, k, T):
    """Intervalo central al 98%: (P1, P99)."""
    q = Puig_quantile([0.01, 0.99], lam, k, T)
    return (q[0], q[1])


def Puig_ci996(lam, k, T):
    """Intervalo central al 99.8%: (P0.1, P99.9)."""
    q = Puig_quantile([0.001, 0.999], lam, k, T)
    return (q[0], q[1])


def Puig_entropy(dist, N=400):
    """Entropía diferencial H(X) = -∫ f(x)·log f(x) dx (trapecio)."""
    P1, P99 = Puig_ci996(dist.lam, dist.k, dist.T)
    init = max(P1 / 10.0, 1e-6)
    fin = P99
    x = np.linspace(init, fin, int(N))
    f = Puig_pdf(x, dist.lam, dist.k, dist.T)
    lf = Puig_logpdf(x, dist.lam, dist.k, dist.T)
    F = -f * lf
    h = (fin - init) / (int(N) - 1)
    FH = F * h
    FH[0] /= 2.0
    FH[-1] /= 2.0
    return float(np.sum(FH))


# ---------------------------------------------------------------------------
# Muestreo exacto
# ---------------------------------------------------------------------------

def Puig_rand(size, lam, k, T, rng=None):
    """Muestreo pseudoaleatorio exacto.

    Representación exacta: X = √(W/T) con W ~ NoncentralChisq(k, λ²T); para
    λ = 0, X = √(W/T) con W ~ Gamma(k/2, 2/T).
    """
    lam, k, T = _check_params(lam, k, T)
    size = tuple(size) if hasattr(size, "__iter__") else (int(size),)
    if lam == 0:
        W = rng.gamma(shape=k / 2.0, scale=2.0 / T, size=size) if rng is not None \
            else np.random.gamma(shape=k / 2.0, scale=2.0 / T, size=size)
        return np.sqrt(W / T)
    W = ncx2.rvs(df=k, nc=lam ** 2 * T, size=size, random_state=rng)
    return np.sqrt(W / T)


def Puig_rand_sample(params, size, rng=None):
    """Compatibilidad: muestreo desde un objeto ``PuigDistribution``."""
    return Puig_rand(size, params.lam, params.k, params.T, rng=rng)