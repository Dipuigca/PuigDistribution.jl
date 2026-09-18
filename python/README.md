# puigdist (Python)

Distribución de **Puig** (χ no central generalizada con dimensión real `k`),
port a Python/NumPy/SciPy del paquete Julia [DistributionsPuig.jl](../../).
Validada contra los *golden files* generados por Julia (`../golden/`).

## Instalación

```bash
pip install .
```

Dependencias: `numpy`, `scipy`, `mpmath`, `pandas`, `matplotlib`.

## Uso rápido

```python
from puigdist import PuigDistribution, Puig_fit
from puigdist.dist import puig   # integración scipy.stats.rv_continuous

d = PuigDistribution(3.0, 2.0, 1.0)         # λ, k, T
d.pdf(2.0)                                  # 0.20216568513008343
d.survival([2.0, 5.0])                      # 0.88672075  0.0306776
d.quantile(0.95)                            # 4.777225027515669
d.moments(8)                                # μ₁ … μ₈
d.rvs(500, rng=42)

ff = puig(3.0, 2.0, 1.0)                    # scipy frozen
ff.pdf([1.0, 2.0]); ff.cdf(5.0); ff.ppf(0.95); ff.stats(moments="mvsk")
```

Métodos: `Puig_pdf(..., method="asymp"|"arb")`, `Puig_cumulative`,
`Puig_surviving`, `Puig_rand`, `Puig_entropy`, `Puig_stats`.

### Ajuste

```python
import pandas as pd, numpy as np
df = pd.DataFrame({ ... variables ... })
res = Puig_fit(df, method="dynamic")     # dynamic | fixed | numerical
res.params.lam, res.params.k, res.params.T
res.indicator     # normas observadas
res.stats         # R², MAE, RMSE, MaxAE, IAE
lam, k, T, s, R2 = res                   # desempaquetado
```

`Puig_fit` acepta `DataFrame` (con `vars`, `time_col`), `ndarray` 2D y
vectores 1D (con `k_fixed`). `effective_dimension(M)` y `empirical_survival`
están disponibles. En el ajuste 1D 3-parámetros, `k` se acota a `[1, 10]`
(diferencia documentada en `../PORTING.md`).

## Tests

```bash
python -m pytest tests -q            # 28 unitarios (port de test/runtests.jl)
python tests/test_golden_core.py     # 4102 comprobaciones contra golden/ Julia
```