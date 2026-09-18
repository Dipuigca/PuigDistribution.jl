"""Ajuste de datos a la distribución de Puig.

Port del módulo ``src/puig_fit.jl``. Proporciona:

- ``effective_dimension``  : dimensión efectiva k = d²/‖R‖_F² (O(d²), sin autovalores).
- ``empirical_survival``   : curva de supervivencia empírica sobre una rejilla.
- métricas ``calc_*``      : R², MAE, RMSE, MaxAE, IAE.
- ``Puig_fit``             : ajuste desde DataFrame (pandas), Matrix (numpy) o
  Vector 1D, con los métodos ``:dynamic``, ``:fixed`` y ``:numerical``.
"""

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .core import PuigDistribution, Puig_surviving

__all__ = [
    "PuigStats", "PuigFitResult", "Puig_fit",
    "effective_dimension", "empirical_survival",
    "calc_r2", "calc_mae", "calc_rmse", "calc_maxae", "calc_iae",
    "calc_all_metrics",
]


# ---------------------------------------------------------------------------
# Métricas de bondad de ajuste (sobre la supervivencia)
# ---------------------------------------------------------------------------

class PuigStats:
    """Métricas de bondad de ajuste frente a la supervivencia empírica.

    Atributos: ``R2`` (coeficiente de determinación), ``MAE``, ``RMSE``,
    ``MaxAE``, ``IAE`` (error absoluto integrado, regla del trapecio).
    """

    def __init__(self, R2, MAE, RMSE, MaxAE, IAE):
        self.R2 = float(R2)
        self.MAE = float(MAE)
        self.RMSE = float(RMSE)
        self.MaxAE = float(MaxAE)
        self.IAE = float(IAE)

    def __iter__(self):
        yield self.R2
        yield self.MAE
        yield self.RMSE
        yield self.MaxAE
        yield self.IAE

    def __repr__(self):
        return (f"PuigStats(R²={self.R2:.5f}, MAE={self.MAE:.5f}, "
                f"RMSE={self.RMSE:.5f}, MaxAE={self.MaxAE:.5f}, IAE={self.IAE:.5f})")


def calc_r2(S_emp, S_mod):
    """Coeficiente de determinación R² = 1 - Σ(S_emp−S_mod)² / Σ(S_emp−S̄_emp)²."""
    S_emp = np.asarray(S_emp, float)
    S_mod = np.asarray(S_mod, float)
    denom = np.sum((S_emp - np.mean(S_emp)) ** 2)
    return float(1.0 - np.sum((S_emp - S_mod) ** 2) / denom)


def calc_mae(S_emp, S_mod):
    """Error absoluto medio MAE = mean(|S_emp − S_mod|)."""
    return float(np.mean(np.abs(np.asarray(S_emp, float) - np.asarray(S_mod, float))))


def calc_rmse(S_emp, S_mod):
    """Raíz del error cuadrático medio RMSE = sqrt(mean((S_emp − S_mod)²))."""
    return float(np.sqrt(np.mean((np.asarray(S_emp, float) - np.asarray(S_mod, float)) ** 2)))


def calc_maxae(S_emp, S_mod):
    """Error absoluto máximo MaxAE = max(|S_emp − S_mod|)."""
    return float(np.max(np.abs(np.asarray(S_emp, float) - np.asarray(S_mod, float))))


def calc_iae(x, S_emp, S_mod):
    """Error absoluto integrado IAE = ∫|S_emp − S_mod| dx (regla del trapecio)."""
    x = np.asarray(x, float)
    diffs = np.abs(np.asarray(S_emp, float) - np.asarray(S_mod, float))
    dx = np.diff(x)
    return float(np.sum(dx * (diffs[:-1] + diffs[1:]) / 2.0))


def calc_all_metrics(x, S_emp, S_mod):
    """Calcula R², MAE, RMSE, MaxAE e IAE a la vez y devuelve un ``PuigStats``."""
    return PuigStats(
        calc_r2(S_emp, S_mod),
        calc_mae(S_emp, S_mod),
        calc_rmse(S_emp, S_mod),
        calc_maxae(S_emp, S_mod),
        calc_iae(x, S_emp, S_mod),
    )


# ---------------------------------------------------------------------------
# Dimensión efectiva y supervivencia empírica
# ---------------------------------------------------------------------------

def effective_dimension(M):
    """Dimensión efectiva k = d²/(d + 2·Σ_{i<j} R_ij²) en O(d²) sin autovalores."""
    M = np.asarray(M, float)
    _, d = M.shape
    if d <= 1:
        return 1.0
    R = np.corrcoef(M, rowvar=False)
    iu = np.triu_indices(d, k=1)
    s_off = np.sum(R[iu] ** 2)
    denom = d + 2.0 * s_off
    return 1.0 if denom == 0.0 else float(d) ** 2 / denom


def empirical_survival(times, grid):
    """S_emp(x) = #{t_i ≥ x} / N evaluada sobre la rejilla ``grid``."""
    times = np.asarray(times, float)
    grid = np.asarray(grid, float)
    n = times.size
    if n == 0:
        return np.zeros_like(grid)
    return np.array([np.count_nonzero(times >= g) / n for g in grid])


# ---------------------------------------------------------------------------
# Contenedor de resultado
# ---------------------------------------------------------------------------

class PuigFitResult:
    """Resultado del ajuste: ``params`` (PuigDistribution), ``indicator``
    (vector de normas/indicador) y ``stats`` (PuigStats).

    Soporta acceso por atributos y desempaquetado de tres elementos::

        params, indicator, stats = Puig_fit(...)
    """

    def __init__(self, params, indicator, stats):
        self.params = params
        self.indicator = np.asarray(indicator, float)
        self.stats = stats

    def __iter__(self):
        yield self.params
        yield self.indicator
        yield self.stats

    def __getitem__(self, i):
        return (self.params, self.indicator, self.stats)[i]

    def __len__(self):
        return 3

    def __repr__(self):
        lo, hi = self.indicator.min(), self.indicator.max()
        return (f"PuigFitResult(λ={self.params.lam:.4f}, k={self.params.k:.4f}, "
                f"T={self.params.T:.5f})\n"
                f"  indicador: {self.indicator.size} obs (rango [{lo:.3f}, {hi:.3f}])\n"
                f"  métricas: {self.stats!r}")


# ---------------------------------------------------------------------------
# Funciones objetivo de la optimización
# ---------------------------------------------------------------------------

def _err_l1(B, x_opt, S_emp, lam, k_target):
    """Error L1 (matriz) entre S empírica y modelo con (λ, k_target, T=1/B²)."""
    if B <= 0:
        return 1e12
    T_val = 1.0 / B ** 2
    try:
        S_mod = Puig_surviving(x_opt, lam, k_target, T_val)
        return float(np.sum(np.abs(S_emp - np.asarray(S_mod, float))))
    except Exception:
        return 1e12


def _err_l2_matrix(par, x_opt, S_emp, lam):
    """Error L2 (matriz) en la parametrización (k, B) con λ fijo."""
    k_val, B_val = par[0], par[1]
    if k_val < 1.0 or B_val <= 0.0:
        return 1e12
    T_val = 1.0 / B_val ** 2
    try:
        S_mod = Puig_surviving(x_opt, lam, k_val, T_val)
        S_mod = np.asarray(S_mod, float)
        if not np.all(np.isfinite(S_mod)):
            return 1e12
        return float(np.sum((S_emp - S_mod) ** 2))
    except Exception:
        return 1e12


# ---------------------------------------------------------------------------
# Ajuste desde Matrix
# ---------------------------------------------------------------------------

def _fit_matrix(Data, times=None, method="dynamic"):
    """Ajuste desde matriz N×d (filas = observaciones). Ver ``Puig_fit``."""
    Data = np.asarray(Data, float)
    n_obs, d = Data.shape
    if n_obs < 10:
        raise ValueError(f"Se requieren al menos 10 observaciones para el ajuste (recibidas: {n_obs})")

    if times is None:
        t_raw = np.sqrt(np.sum(Data ** 2, axis=1))
    else:
        t_raw = np.asarray(times, float)

    if t_raw.shape[0] != n_obs:
        raise ValueError("times debe tener tantas observaciones como filas de Data")

    valid = np.isfinite(t_raw) & (t_raw >= 0)
    t = t_raw[valid]
    if t.size < 10:
        raise ValueError("Se requieren al menos 10 observaciones válidas para el ajuste")

    col_means = np.mean(Data[valid, :], axis=0)
    lam = np.sqrt(np.sum(col_means ** 2))
    k_corr = effective_dimension(Data[valid, :])

    min_t, max_t = t.min(), t.max()
    x_opt = np.linspace(min_t, max_t, 150)
    S_emp_opt = empirical_survival(t, x_opt)

    sdt = np.std(t, ddof=1)
    B0 = float(sdt) if (np.isfinite(sdt) and sdt > 0) else 1.0

    k_est, B_est = k_corr, B0

    if method in ("dynamic", "fixed"):
        k_target = float(d) if method == "fixed" else k_corr
        res = minimize_scalar(
            lambda B: _err_l1(B, x_opt, S_emp_opt, lam, k_target),
            bounds=(0.01, 5.0 * B0),
            method="bounded",
            options={"xatol": 1e-12},
        )
        B_est = res.x
        k_est = k_target
    elif method == "numerical":
        lower = np.array([1.0, 0.001])
        upper = np.array([max(3.0, float(d)), 10.0 * B0])
        init_par = np.array([np.clip(k_corr, 1.0, upper[0]), B0])
        res = minimize(
            lambda p: _err_l2_matrix(p, x_opt, S_emp_opt, lam),
            init_par,
            method="L-BFGS-B",
            bounds=list(zip(lower, upper)),
        )
        if res.success:
            k_est, B_est = res.x[0], res.x[1]
        else:
            k_est = np.clip(k_corr, 1.0, float(d))
            B_est = B0
    else:
        raise ValueError("Método desconocido :%s. Opciones válidas: :dynamic, :numerical, :fixed" % method)

    T_est = 1.0 / B_est ** 2
    x_eval = np.linspace(min_t, max_t, 500)
    S_emp_eval = empirical_survival(t, x_eval)
    S_mod_eval = np.asarray(Puig_surviving(x_eval, lam, k_est, T_est), float)
    fit_metrics = calc_all_metrics(x_eval, S_emp_eval, S_mod_eval)

    params = PuigDistribution(float(lam), float(k_est), float(T_est))
    return PuigFitResult(params, t, fit_metrics)


# ---------------------------------------------------------------------------
# Ajuste desde DataFrame (pandas)
# ---------------------------------------------------------------------------

def _fit_dataframe(df, vars=None, time_col=None, method="dynamic"):
    """Ajuste desde DataFrame pandas. Ver ``Puig_fit``."""
    import pandas as pd

    var_names = vars
    if var_names is None:
        exclude = set() if time_col is None else {time_col}
        var_names = [
            c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]) and c not in exclude
        ]
    if not var_names:
        raise ValueError("Se requiere al menos una columna numérica para el ajuste")

    mat = df[var_names].to_numpy(dtype=float)
    if mat.ndim != 2:
        raise ValueError("vars debe seleccionar más de una columna numérica")

    times = df[time_col].to_numpy(dtype=float) if time_col is not None else None
    return _fit_matrix(mat, times=times, method=method)


# ---------------------------------------------------------------------------
# Ajuste desde Vector 1D
# ---------------------------------------------------------------------------

def _err_2d_fixed(p, x_opt, S_emp, k_target):
    """Error L2 (vector 1D) en (log λ, log B) con k fijo."""
    lam_val = np.exp(p[0])
    B_val = np.exp(p[1])
    T_val = 1.0 / B_val ** 2
    try:
        S_mod = Puig_surviving(x_opt, lam_val, k_target, T_val)
        S_mod = np.asarray(S_mod, float)
        if not np.all(np.isfinite(S_mod)):
            return 1e12
        return float(np.sum((S_emp - S_mod) ** 2))
    except Exception:
        return 1e12


def _err_3d(p, x_opt, S_emp):
    """Error L2 (vector 1D) en (log λ, log(k−1), log B)."""
    lam_val = np.exp(p[0])
    k_val = 1.0 + np.exp(p[1])
    B_val = np.exp(p[2])
    if k_val > 10.0:
        return 1e12
    T_val = 1.0 / B_val ** 2
    try:
        S_mod = Puig_surviving(x_opt, lam_val, k_val, T_val)
        S_mod = np.asarray(S_mod, float)
        if not np.all(np.isfinite(S_mod)):
            return 1e12
        return float(np.sum((S_emp - S_mod) ** 2))
    except Exception:
        return 1e12


def _fit_vector(Data, method="dynamic", k_fixed=None):
    """Ajuste desde vector 1D de normas observadas. Ver ``Puig_fit``."""
    Data = np.asarray(Data, float)
    valid = np.isfinite(Data) & (Data >= 0)
    t = Data[valid]
    n_obs = t.size
    if n_obs < 10:
        raise ValueError(f"Se requieren al menos 10 observaciones para el ajuste (recibidas: {n_obs})")

    min_t, max_t = t.min(), t.max()
    x_opt = np.linspace(min_t, max_t, 60)
    S_emp_opt = empirical_survival(t, x_opt)

    m1 = np.mean(t)
    m2 = np.mean(t ** 2)
    s = np.std(t, ddof=1)
    B0 = float(s) if (np.isfinite(s) and s > 0) else 1.0
    lam0 = np.sqrt(max(0.0, m1 ** 2 - s ** 2))
    k0 = float(np.clip((m2 - lam0 ** 2) / B0 ** 2, 1.0, 10.0))

    is_fixed_k = (method == "fixed") or (k_fixed is not None)
    k_target = float(k_fixed) if k_fixed is not None else 1.0
    if k_target < 1.0:
        raise ValueError(f"k_fixed debe ser >= 1.0 (recibido: {k_target})")

    k_est, lam_est, B_est = 1.0, lam0, B0

    if is_fixed_k:
        k_est = k_target
        init_p = np.array([np.log(max(1e-3, lam0)), np.log(max(1e-3, B0))])
        res = minimize(
            lambda p: _err_2d_fixed(p, x_opt, S_emp_opt, k_target),
            init_p,
            method="Nelder-Mead",
            options={"maxiter": 500, "xatol": 1e-8, "fatol": 1e-8},
        )
        p_opt = res.x
        lam_est = np.exp(p_opt[0])
        B_est = np.exp(p_opt[1])
    elif method in ("dynamic", "numerical"):
        init_p = np.array([np.log(max(1e-3, lam0)), np.log(max(1e-3, k0 - 1.0)),
                           np.log(max(1e-3, B0))])
        res = minimize(
            lambda p: _err_3d(p, x_opt, S_emp_opt),
            init_p,
            method="Nelder-Mead",
            options={"maxiter": 600, "xatol": 1e-8, "fatol": 1e-8},
        )
        p_opt = res.x
        lam_est = np.exp(p_opt[0])
        k_est = 1.0 + np.exp(p_opt[1])
        B_est = np.exp(p_opt[2])
    else:
        raise ValueError("Método desconocido :%s. Opciones válidas: :dynamic, :numerical, :fixed" % method)

    T_est = 1.0 / B_est ** 2
    x_eval = np.linspace(min_t, max_t, 500)
    S_emp_eval = empirical_survival(t, x_eval)
    S_mod_eval = np.asarray(Puig_surviving(x_eval, lam_est, k_est, T_est), float)
    fit_metrics = calc_all_metrics(x_eval, S_emp_eval, S_mod_eval)

    params = PuigDistribution(float(lam_est), float(k_est), float(T_est))
    return PuigFitResult(params, t, fit_metrics)


# ---------------------------------------------------------------------------
# API principal
# ---------------------------------------------------------------------------

def Puig_fit(Data, vars=None, time_col=None, times=None, method="dynamic", k_fixed=None):
    """Ajusta datos a la distribución de Puig.

    Entradas aceptadas:
    - ``pandas.DataFrame``: columnas numéricas (o las de ``vars``); el indicador
      se calcula como ``√(Σ vars²)`` o se toma de ``time_col``.
    - ``numpy.ndarray`` 2D (filas = observaciones, columnas = variables).
    - ``numpy.ndarray`` 1D / lista (serie de normas observadas).

    Métodos (matriz/DataFrame): ``:dynamic`` (k por correlación), ``:fixed``
    (k = d), ``:numerical`` (k y B conjuntos, L-BFGS-B).
    Métodos (vector 1D): ``:dynamic``/``:numerical`` (3D conjunta Nelder-Mead)
    o ``:fixed`` con ``k_fixed``.
    """
    import pandas as pd

    if isinstance(Data, pd.DataFrame):
        if k_fixed is not None:
            raise ValueError("k_fixed solo aplica al ajuste de vectores 1D")
        return _fit_dataframe(Data, vars=vars, time_col=time_col, method=method)

    arr = np.asarray(Data, float)
    if arr.ndim == 2:
        if k_fixed is not None:
            raise ValueError("k_fixed solo aplica al ajuste de vectores 1D")
        return _fit_matrix(arr, times=times, method=method)
    if arr.ndim == 1:
        return _fit_vector(arr, method=method, k_fixed=k_fixed)
    raise ValueError(f"Entrada no soportada (revisa dimensiones): shape={arr.shape}")