"""Tests unitarios del núcleo, momentos, Marcum-Q y ajuste.

Port de la suite ``test/runtests.jl`` del paquete Julia.
"""

import math

import numpy as np
import pandas as pd
import pytest
from scipy.special import gamma

from puigdist.core import (
    PuigDistribution, MB, MB_pdf, MB_surviving, MB_cumulative,
    Puig_pdf, Puig_logpdf, Puig_surviving, Puig_cumulative,
    Puig_mean, Puig_var, Puig_std, Puig_skewness, Puig_kurtosis,
    Puig_moments, Puig_stats, Puig_quantile, Puig_rand,
    marcumq, marcumq_asymp,
)
from puigdist.fit import (
    effective_dimension, empirical_survival,
    calc_all_metrics, PuigStats, Puig_fit, PuigFitResult,
)


# ============================================================
# 0. Módulo y tipos
# ============================================================
def test_modulo_tipos():
    d = PuigDistribution(3.0, 2.0, 1.0)
    assert isinstance(d, PuigDistribution)
    assert d.lam == 3.0
    assert d.k == 2.0
    assert d.T == 1.0

    mb = MB(2.0, 1.0)
    assert isinstance(mb, MB)
    assert mb.k == 2.0
    assert mb.B == 1.0


# ============================================================
# 1. effective_dimension
# ============================================================
def test_effective_dimension():
    rng = np.random.RandomState(0)
    M1 = rng.randn(50, 1)
    assert effective_dimension(M1) == pytest.approx(1.0)

    rng = np.random.RandomState(123)
    M_uncorr = rng.randn(10000, 3)
    assert effective_dimension(M_uncorr) == pytest.approx(3.0, abs=0.1)

    x_base = rng.randn(100)
    M_corr = np.column_stack([x_base, x_base, x_base])
    assert effective_dimension(M_corr) == pytest.approx(1.0, abs=1e-5)

    d = 3
    s_off = 0.8 ** 2 + 0.6 ** 2 + 0.7 ** 2
    k_exacto = float(d) ** 2 / (d + 2.0 * s_off)
    assert k_exacto == pytest.approx(1.5050167224, abs=1e-6)


# ============================================================
# 2. empirical_survival
# ============================================================
def test_empirical_survival():
    t = [1.0, 2.0, 3.0, 4.0, 5.0]
    grid = [0.0, 1.0, 3.0, 5.0, 6.0]
    s = empirical_survival(t, grid)
    assert s[0] == 1.0
    assert s[1] == 1.0
    assert s[2] == pytest.approx(0.6)
    assert s[3] == pytest.approx(0.2)
    assert s[4] == 0.0
    assert np.all(np.diff(s) <= 0.0)


# ============================================================
# 3. Métricas y PuigStats
# ============================================================
def test_metricas():
    S_true = [1.0, 0.8, 0.5, 0.2, 0.0]
    x_grid = [1.0, 2.0, 3.0, 4.0, 5.0]
    st_perf = calc_all_metrics(x_grid, S_true, S_true)
    assert isinstance(st_perf, PuigStats)
    assert st_perf.R2 == pytest.approx(1.0)
    assert st_perf.MAE == pytest.approx(0.0)
    assert st_perf.RMSE == pytest.approx(0.0)
    assert st_perf.MaxAE == pytest.approx(0.0)
    assert st_perf.IAE == pytest.approx(0.0)

    S_approx = np.array(S_true) + 0.05
    st_err = calc_all_metrics(x_grid, S_true, S_approx)
    assert st_err.MAE == pytest.approx(0.05)
    assert st_err.MaxAE == pytest.approx(0.05)
    assert st_err.RMSE == pytest.approx(0.05)
    assert st_err.R2 < 1.0


# ============================================================
# 4. Puig_fit sobre DataFrame / Matrix / Vector
# ============================================================
def test_puig_fit_dataframe():
    rng = np.random.RandomState(42)
    N = 150
    X1 = 150.0 + rng.randn(N) * 20.0
    X2 = 100.0 + 0.8 * X1 + rng.randn(N) * 10.0
    X3 = 200.0 + 0.5 * X1 + rng.randn(N) * 15.0
    df = pd.DataFrame({"Squat": X1, "Bench": X2, "Deadlift": X3})

    res_dyn = Puig_fit(df, method="dynamic")
    assert isinstance(res_dyn, PuigFitResult)
    assert isinstance(res_dyn.params, PuigDistribution)
    assert isinstance(res_dyn.indicator, np.ndarray)
    assert isinstance(res_dyn.stats, PuigStats)
    assert len(res_dyn.indicator) == N
    assert res_dyn.params.lam > 0
    assert 1.0 <= res_dyn.params.k <= 3.0
    assert res_dyn.params.T > 0
    assert res_dyn.stats.R2 > 0.95

    res_fix = Puig_fit(df, method="fixed")
    assert res_fix.params.k == pytest.approx(3.0, abs=1e-12)
    assert res_fix.stats.R2 > 0.90

    res_num = Puig_fit(df, method="numerical")
    assert 1.0 <= res_num.params.k <= 3.0
    assert res_num.stats.R2 > 0.95

    res_sub = Puig_fit(df, vars=["Squat", "Bench"])
    assert res_sub.params.k <= 2.0
    assert len(res_sub.indicator) == N

    df_with_time = df.copy()
    df_with_time["Total"] = np.sqrt(df.Squat ** 2 + df.Bench ** 2 + df.Deadlift ** 2)
    res_time = Puig_fit(df_with_time, vars=["Squat", "Bench", "Deadlift"], time_col="Total")
    assert np.array_equal(res_time.indicator, np.asarray(df_with_time["Total"], float))

    mat = df.to_numpy(dtype=float)
    res_mat = Puig_fit(mat, method="dynamic")
    assert res_mat.params.lam == pytest.approx(res_dyn.params.lam, rel=1e-6)
    assert res_mat.params.k == pytest.approx(res_dyn.params.k, rel=1e-6)
    assert res_mat.params.T == pytest.approx(res_dyn.params.T, rel=1e-6)

    res_vec = Puig_fit(np.asarray(df_with_time["Total"], float))
    assert isinstance(res_vec, PuigFitResult)
    assert len(res_vec.indicator) == N
    assert res_vec.params.lam > 0
    assert res_vec.params.k >= 1.0
    assert res_vec.params.T > 0
    assert res_vec.stats.R2 > 0.95


def test_puig_fit_1d_sintetico():
    rng = np.random.RandomState(77)
    dist_sint = PuigDistribution(40.0, 2.5, 0.05)
    muestra = Puig_rand(500, 40.0, 2.5, 0.05, rng=rng)

    res_opt = Puig_fit(muestra)
    assert res_opt.params.lam == pytest.approx(40.0, abs=5.0)
    assert res_opt.params.k >= 1.0
    assert res_opt.stats.R2 > 0.98

    res_kfix = Puig_fit(muestra, k_fixed=2.5)
    assert res_kfix.params.k == pytest.approx(2.5, abs=1e-9)
    assert res_kfix.stats.R2 > 0.98


# ============================================================
# 5. Desempaquetado de PuigFitResult
# ============================================================
def test_desempaquetado():
    rng = np.random.RandomState(99)
    df_test = pd.DataFrame({
        "A": 10.0 + rng.randn(50),
        "B": 20.0 + rng.randn(50),
    })
    res = Puig_fit(df_test)

    p, ind, st = res
    assert p is res.params
    assert ind is res.indicator
    assert st is res.stats
    assert res[0] is res.params
    assert res[1] is res.indicator
    assert res[2] is res.stats
    assert len(res) == 3


# ============================================================
# 6. Manejo de errores y casos límite
# ============================================================
def test_errores():
    with pytest.raises(ValueError):
        Puig_fit(pd.DataFrame({"A": np.random.rand(5), "B": np.random.rand(5)}))
    with pytest.raises(ValueError):
        Puig_fit(np.random.rand(5))
    with pytest.raises(ValueError):
        Puig_fit(pd.DataFrame({"A": np.random.rand(20), "B": np.random.rand(20)}),
                 method="metodo_inexistente")
    with pytest.raises(ValueError):
        Puig_fit(np.random.rand(20), method="metodo_inexistente")
    with pytest.raises(ValueError):
        Puig_fit(np.random.rand(20), k_fixed=0.5)
    df_non_num = pd.DataFrame({"Nombre": list("ABCDEFGHIJK")})
    with pytest.raises(ValueError):
        Puig_fit(df_non_num)


# ============================================================
# 7. Interfaz tipo Distributions (equivalente scipy)
# ============================================================
def test_interfaz_base():
    d = PuigDistribution(3.0, 2.0, 1.0)

    assert isinstance(d.pdf(1.0), float)
    assert d.pdf(0.0) == 0.0
    assert isinstance(d.logpdf(1.0), float)
    assert d.logpdf(0.0) == -math.inf

    assert d.cumulative(0.0) == 0.0
    assert d.survival(0.0) == 1.0
    assert d.cumulative(100.0) == pytest.approx(1.0, abs=1e-6)
    assert d.survival(100.0) == pytest.approx(0.0, abs=1e-6)
    assert d.cumulative(5.0) + d.survival(5.0) == pytest.approx(1.0, abs=1e-10)

    assert isinstance(d.quantile(0.5), float)
    assert d.quantile(0.0) == 0.0
    assert d.quantile(1.0) == math.inf

    m = d.mean()
    v = d.var()
    s = d.std()
    assert v == pytest.approx(s ** 2, abs=1e-6)
    assert m > 0
    assert v > 0

    assert isinstance(d.skewness(), float)
    assert isinstance(d.kurtosis(), float)

    rng = np.random.RandomState(42)
    samples = Puig_rand(100, 3.0, 2.0, 1.0, rng=rng)
    assert samples.shape[0] == 100
    assert np.all(samples >= 0)
    assert np.mean(samples) == pytest.approx(m, abs=1.0)

    lam, k, T = d
    assert lam == 3.0
    assert k == 2.0
    assert T == 1.0


def test_interfaz_scipy():
    from puigdist.dist import puig
    ff = puig(3.0, 2.0, 1.0)

    assert ff.pdf(1.0) == pytest.approx(0.03288652, abs=1e-7)
    assert ff.cdf(5.0) + ff.sf(5.0) == pytest.approx(1.0, abs=1e-10)
    assert ff.ppf(0.5) == pytest.approx(ff.median(), rel=1e-6)
    assert ff.support()[0] == 0.0
    assert math.isinf(ff.support()[1])

    rng = np.random.RandomState(7)
    rvs = ff.rvs(size=500, random_state=rng)
    assert rvs.shape == (500,)
    assert np.all(rvs >= 0)


# ============================================================
# 8. Robustez asintótica vs alta precisión (mpmath)
# ============================================================
def test_robustez_asymp_vs_arb():
    for lam in [10.0, 50.0, 100.0, 250.0, 500.0, 1000.0]:
        kv, Tv = 2.4, 0.01
        for frac in [0.7, 0.9, 1.0, 1.1, 1.3]:
            xv = frac * lam
            v_asymp = Puig_pdf(xv, lam, kv, Tv, method="asymp")
            v_arb = Puig_pdf(xv, lam, kv, Tv, method="arb")
            if v_arb > 1e-15:
                assert abs(v_asymp - v_arb) / v_arb < 1e-12
            else:
                assert abs(v_asymp - v_arb) < 1e-15

    for kv in [1.0001, 1.2, 2.0, 3.0, 5.0, 10.0]:
        lam, Tv = 50.0, 0.05
        for frac in [0.8, 1.0, 1.2]:
            xv = frac * lam
            v_asymp = Puig_pdf(xv, lam, kv, Tv, method="asymp")
            v_arb = Puig_pdf(xv, lam, kv, Tv, method="arb")
            if v_arb > 1e-15:
                assert abs(v_asymp - v_arb) / v_arb < 1e-12
            else:
                assert abs(v_asymp - v_arb) < 1e-15


def test_pdf_limites_extremos():
    lam, kv, Tv = 100.0, 2.5, 0.05
    assert Puig_pdf(0.0, lam, kv, Tv, method="asymp") == 0.0
    with pytest.raises(ValueError):
        Puig_pdf(-5.0, lam, kv, Tv, method="asymp")
    assert Puig_pdf(1000.0, lam, kv, Tv, method="asymp") == 0.0
    assert math.isfinite(Puig_pdf(100.0, lam, kv, Tv, method="asymp"))
    assert Puig_pdf(100.0, lam, kv, Tv, method="asymp") > 0.0

    p_a = Puig_pdf(50.0, 50.0, 2.4, 0.05, method="asymp")
    p_b = Puig_pdf(50.0, 50.0, 2.4, 0.05, method="arb")
    assert p_a == pytest.approx(p_b, abs=1e-12)
    p_a2 = Puig_pdf(100.0, 50.0, 2.4, 0.05, method="asymp")
    p_b2 = Puig_pdf(100.0, 50.0, 2.4, 0.05, method="arb")
    assert p_a2 == pytest.approx(p_b2, abs=1e-12)


# ============================================================
# 9. Momentos genéricos de orden n
# ============================================================
def test_momentos():
    d9 = PuigDistribution(2.5, 3.2, 1.4)
    tuple4 = Puig_moments(2.5, 3.2, 1.4)
    vec4 = Puig_moments(2.5, 3.2, 1.4, 4)
    assert len(vec4) == 4
    for i in range(4):
        assert vec4[i] == pytest.approx(tuple4[i])

    vec8 = Puig_moments(2.5, 3.2, 1.4, 8)
    assert len(vec8) == 8
    for i in range(4):
        assert vec8[i] == pytest.approx(vec4[i])

    mu1, mu2, mu3, mu4 = vec4
    assert np.all(vec8[1::2] > 0)

    lamf, kf, Tf = 2.5, 3.2, 1.4
    mu5_rec = ((10.0 - 4.0 + kf) / Tf + lamf ** 2) * mu3 - 3.0 * (1.0 + kf) / Tf ** 2 * mu1
    assert vec8[4] == pytest.approx(mu5_rec, abs=1e-10)
    mu6_rec = ((12.0 - 4.0 + kf) / Tf + lamf ** 2) * mu4 - 4.0 * (2.0 + kf) / Tf ** 2 * mu2
    assert vec8[5] == pytest.approx(mu6_rec, abs=1e-10)

    B_mb, k_mb = 1.0, 2.0
    vec_mb = Puig_moments(0.0, k_mb, 1.0 / B_mb, 8)
    for nn in range(1, 9):
        mu_analitico = (2.0 * B_mb) ** (nn / 2) * gamma((k_mb + nn) / 2) / gamma(k_mb / 2)
        assert vec_mb[nn - 1] == pytest.approx(mu_analitico, abs=1e-10)

    assert Puig_moments(2.5, 3.2, 1.4, 1)[0] == pytest.approx(mu1)
    assert Puig_moments(2.5, 3.2, 1.4, 2)[1] == pytest.approx(mu2)
    assert Puig_moments(2.5, 3.2, 1.4, 3)[2] == pytest.approx(mu3)

    with pytest.raises(ValueError):
        Puig_moments(2.5, 3.2, 1.4, 0)

    st = Puig_stats(2.5, 3.2, 1.4)
    assert st["mean"] == pytest.approx(vec4[0])
    assert st["var"] == pytest.approx(vec4[1] - vec4[0] ** 2)


def test_mb_momentos_analiticos():
    mb = MB(2.0, 1.0)
    vec = mb.moments(8)
    for nn in range(1, 9):
        mu = (2.0) ** (nn / 2) * gamma((2.0 + nn) / 2) / gamma(1.0)
        assert vec[nn - 1] == pytest.approx(mu, abs=1e-10)


# ============================================================
# 10. Marcum-Q asintótico y supervivencia en colas extremas
# ============================================================
def test_marcumq_asymp():
    v = marcumq_asymp(10.0, 15.0, 1.0)
    assert 0.0 < v < 1.0

    lam_t, k_t, T_t = 10.0, 3.2, 1.0
    sT = math.sqrt(T_t)
    a_t = lam_t * sT
    b_extreme = 25.0 * sT
    mq_val = marcumq(a_t, b_extreme, k_t / 2)
    assert mq_val < 1e-30

    x_extreme = [lam_t + 10.0, lam_t + 15.0, lam_t + 20.0]
    S = Puig_surviving(x_extreme, lam_t, k_t, T_t)
    assert np.all(np.isfinite(S))
    assert np.all(S >= 0.0)
    assert np.all(np.diff(S) <= 0.0)


def test_marcumq_frontera_y_bomba():
    lam_t, k_t, T_t = 10.0, 3.2, 1.0
    sT = math.sqrt(T_t)
    x_boundary = [lam_t + 4.0 / sT]
    S_bnd = Puig_surviving(x_boundary, lam_t, k_t, T_t)
    assert math.isfinite(S_bnd[0])
    assert 0.0 < S_bnd[0] < 1.0

    a_t = lam_t * sT
    x_vec = [lam_t + 5.0, lam_t + 8.0, lam_t + 12.0]
    mq_scalar = [marcumq(a_t, xi * sT, k_t / 2) for xi in x_vec]
    mq_vector = marcumq(a_t, np.array(x_vec) * sT, k_t / 2)
    assert np.allclose(mq_scalar, mq_vector, atol=1e-14)

    x_remote = [lam_t + 10.0, lam_t + 12.0, lam_t + 14.0, lam_t + 16.0, lam_t + 18.0]
    S_remote = Puig_surviving(x_remote, lam_t, k_t, T_t)
    assert np.all(np.isfinite(S_remote))
    assert np.all(np.diff(S_remote) <= 0.0)

    mq_small = marcumq(1.0, 3.0, 1.0)  # ab = 3 < 30 -> fallback ncx2
    assert math.isfinite(mq_small)
    assert 0.0 < mq_small < 1.0

    assert marcumq_asymp(10.0, 20.0, 2.0) > marcumq_asymp(10.0, 20.0, 1.0)


# ============================================================
# 11. MB básico
# ============================================================
def test_mb_limite():
    mb = MB(2.0, 1.0)
    # k=2, B=1: chi-2 escala -> pdf en x: x*exp(-x^2/2)
    assert MB_pdf(0.0, 2.0, 1.0) == 0.0
    assert MB_pdf(1.0, 2.0, 1.0) == pytest.approx(math.exp(-0.5), rel=1e-12)
    assert MB_surviving(0.0, 2.0, 1.0) == pytest.approx(1.0)
    assert MB_cumulative(0.0, 2.0, 1.0) == pytest.approx(0.0)