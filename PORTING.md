---

editor_options: 
  markdown: 
    wrap: 72
---

# Portes a Python, R y MATLAB

Este repositorio contiene, además del paquete Julia **DistributionsPuig.jl** (`src/`), portes funcionales de la distribución de Puig a Python, R y MATLAB, validados de forma cruzada contra **golden files** generados desde Julia (`golden/`).

Todos los lenguajes implementan el mismo núcleo matemático:

- PDF con Bessel modificada y doble camino numérico `:asymp` (expansión asintótica para el producto `exp × I_ν` con argumentos grandes) y `:arb` (referencia de alta precisión: `mpmath` / Rmpfr / doble con piezo-regular).
- Supervivencia y acumulada por la función Marcum-`Q` generalizada (con emulación de dimensionalidad real de la `χ` no central y expansión asintótica para colas extremas).
- Momentos raw genéricos (recurrencia de tres términos O(n)) y estadísticos (media, varianza, σ, asimetría, kurtosis).
- Ajuste `Puig_fit` (L1/L2) desde datos multivariantes y series 1D.
- Muestreo exacto por la relación `X = sqrt(W/T)`, `W ~ NoncentralChisq(k, λ²T)`.
- Dimensión efectiva, supervivencia empírica y métricas de bondad de ajuste.

## Corrección compartida: momentos con λ grande

La función de Laguerre generalizada `laguerre_real(n, α, x)` es el origen del cálculo de μ₁ y μ₃. La implementación original (Julia) evaluaba `L_n^{(α)}(-y) = coef·e^{-y}·M(α+1+n, α+1, y)` mediante la serie de Kummer de `M(a, b, +y)`. Para `y = λ²T/2` enorme (p. ej. λ=1000, T=0.05 → y=25000), esa serie con argumento positivo no converge dentro de `maxterms` y la cancelación produce μ₁ = 0 (falso «modo oscuro» del perfil). El fallo afectaba por igual a Julia, Python (mpmath) y R (Rmpfr); MATLAB ya era correcto.

Se corrigió en los cuatro lenguajes con **dos ramas numéricamente robustas**:

| Rama | Condición | Fórmula |
|:-----------------------|:-----------------------|:-----------------------|
| Serie (Kummer) | y \< 20 | `M(-n, b, -y) = e^{-y}·M(b+n, b, y)` |
| Asintótica (rama decreciente) | y ≥ 20 | `e^{-y}·M(a,b,y) ≈ Γ(b)/Γ(a)·y^{a-b}·Σₖ (1-a)ₖ(b-a)ₖ/k!·y^{-k}` |

La rama asintótica se validó contra `mpmath.hyper` (corrección exacta de los coeficientes) y evita por completo la cancelación catastrófica. El golden de momentos se regeneró con ese esquema.

## Estado de validación (golden 4102 checks + unitarios)

| Lenguaje | Golden | Unitarios | Comando de validación |
|:-----------------|:-----------------|:-----------------|:-----------------|
| **Julia** (referencia) | regenera `golden/` | 205 | `julia --project=. golden/generate.jl` y `julia --project=. -e "using Pkg; Pkg.test()"` |
| **Python** | 4102 PASS | 28 PASS | `python tests/test_golden_core.py` · `python tests/test_core.py` (desde `python/` con `PYTHONPATH=src`) |
| **R** | 4102 PASS | 108 PASS | `Rscript R/tests/run_golden_tests.R` · `Rscript R/tests/run_unit_tests.R` (desde la raíz) |
| **MATLAB** | 4102 PASS | 106 PASS | `matlab -batch "run('matlab/tests/run_golden_tests.m')"` · `matlab -batch "run('matlab/tests/run_unit_tests.m')"` |

Los golden (véase `golden/generate.jl`) cubren λ ∈ {0, 0.5, 3, 10, 50, 100, 250, 1000}, k ∈ {1.0001, 1.5, 2.4, 3, 10}, T ∈ {0.001, 0.01, 0.05, 0.25, 1}, colas de PDF y supervivencia, cuantiles y el caso límite Maxwell–Boltzmann (`MB`).

## Mapa de archivos

```         
python/src/puigdist/   paquete instalable (pyproject.toml): core, mb, hypergeom,
                       highprec, dist (scipy rv_continuous), fit, plot, __init__
R/                     paquete R: core_types, core_density, highprec, moments,
                       quantiles, fit, plot (ggplot2), wrappers (d/p/q/r),
                       DESCRIPTION, NAMESPACE, tests/
matlab/+puigdist/      paquete MATLAB (28 funciones + clases PuigDistribution/MB):
                       Puig_pdf, Puig_surviving, marcumq(+_asymp), besselix_asymp,
                       Puig_moments, Puig_fit, Puig_plot, effective_dimension, ...
matlab/tests/          run_golden_tests.m, run_unit_tests.m
golden/                golden_pdf/surv/cdf/moments/quantile/mb/mb_moments.csv
```

## Notas específicas por lenguaje

- **Python**: `dist.py` expone `PuigDistribution` como `scipy.stats.rv_continuous` (con `_argcheck`, `_pdf`, `_cdf`, `_sf`, `_ppf`, `_rvs`); `fit.py` y `plot.py` replican `Puig_fit` y `Puig_plot`.
- **R**: el paquete implementa la interfaz idiomática `dpuig/ppuig/qpuig/rpuig` además de las funciones `Puig_*` homólogas. Tests por `Rscript` sin dependencias de terceros (solo base + sugeridas ggplot2).
- **MATLAB**: paquete `+puigdist` de funciones independientes (sin dependencias de Toolboxes); la Marcum-`Q` de orden real usa `ncx2cdf`, y el muestreo `ncx2rnd` (acepta df reales en R2025b). `Puig_fit` usa `fminbnd`/`fminsearch`. `run()` cambia el directorio de trabajo al script, por lo que los tests resuelven rutas con `mfilename('fullpath')`; se ejecutan desde la raíz.
