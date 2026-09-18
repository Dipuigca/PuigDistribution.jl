# Referencia de la API

Esta sección contiene la referencia completa de estructuras, funciones y métodos disponibles en `DistributionsPuig.jl`.

---

## 1. Estructuras de Datos (Tipos)

### `PuigDistribution`
```julia
struct PuigDistribution{S<:Real}
    λ::S  # Norma de las medias de las componentes (λ ≥ 0)
    k::S  # Dimensión efectiva continua (k ≥ 1)
    T::S  # Parámetro de precisión / escala (T > 0, T = 1/σ²)
end
```
Representa una distribución de Puig. Admite promoción automática de tipos y formateo visual con `Base.show`.

---

### `MB`
```julia
struct MB{S<:Real}
    k::S  # Dimensión / grados de libertad reales (k ≥ 1)
    B::S  # Parámetro de escala o dispersión térmica (B > 0, B = 1/T)
end
```
Representa la distribución límite de Maxwell-Boltzmann / $\chi$ generalizada cuando $\lambda \to 0$.

---

### `PuigStats`
```julia
struct PuigStats{S<:Real}
    R²::S     # Coeficiente de determinación sobre la función de supervivencia
    MAE::S    # Error absoluto medio
    RMSE::S   # Raíz del error cuadrático medio
    MaxAE::S  # Error absoluto máximo
    IAE::S    # Error absoluto integrado (regla del trapecio)
end
```
Contenedor tipado de métricas de bondad de ajuste frente a la curva de supervivencia empírica.

---

### `PuigFitResult`
```julia
struct PuigFitResult{S<:Real}
    params::PuigDistribution{S}   # Parámetros ajustados de la distribución
    indicator::Vector{S}          # Vector con los valores del indicador compuesto observado
    stats::PuigStats{S}          # Métricas de bondad de ajuste
end
```
Contenedor del resultado completo del ajuste. Soporta acceso mediante propiedades (`res.params`, `res.indicator`, `res.stats`) y desempaquetado directo de tres elementos:
```julia
dist, ind, st = Puig_fit(df)
```

---

## 2. Núcleo Matemático (`puig_pdf.jl` y `Maxwell_Boltzmann.jl`)

### Densidad, Acumulada y Supervivencia
* `Puig_pdf(x, λ, k, T)` / `Puig_pdf(x, params)`: Función de densidad de probabilidad $f_P(x)$.
* `Puig_cumulative(x, λ, k, T)` / `Puig_cumulative(x, params)`: Función de distribución acumulada $F_P(x) = P(X \le x)$.
* `Puig_surviving(x, λ, k, T)` / `Puig_surviving(x, params)`: Función de supervivencia $S_P(x) = P(X \ge x)$.
* `marcumq(a, b, m)`: Función Marcum-$Q$ generalizada de orden $m$. Para la cola superior lejana ($b - a \ge 4$, $ab \ge 30$) conmuta automáticamente a la expansión asintótica `marcumq_asymp`.
* `marcumq_asymp(a, b, m)`: Expansión asintótica de 2 términos de $Q_m(a, b)$ para cola superior lejana ($b - a \ge 4.0$, $ab \ge 30.0$). Cero cancelaciones catastróficas, $\sim 20\text{ ns}$.

### Momentos y Estadísticos Teóricos
* `Puig_mean(λ, k, T)` / `Puig_mean(params)`: Media teórica $\mu_1 = E[X]$ vía función de Laguerre generalizada.
* `Puig_var(λ, k, T)` / `Puig_var(params)`: Varianza teórica $\mathrm{Var}(X) = E[X^2] - (E[X])^2$.
* `Puig_std(λ, k, T)` / `Puig_std(params)`: Desviación estándar teórica $\sigma = \sqrt{\mathrm{Var}(X)}$.
* `Puig_skewness(λ, k, T)` / `Puig_skewness(params)`: Coeficiente de asimetría (*skewness*) estandarizado $\gamma_1$.
* `Puig_kurtosis(λ, k, T)` / `Puig_kurtosis(params)`: *Kurtosis* teórica $\gamma_2 = E[(X-\mu)^4] / \sigma^4$.
* `Puig_moments(λ, k, T)` / `Puig_moments(params)`: Tupla con los 4 primeros momentos raw $(\mu_1, \mu_2, \mu_3, \mu_4)$.
* `Puig_moments(λ, k, T, n)` / `Puig_moments(params, n)`: Vector con los $n$ primeros momentos raw $[\mu_1, \dots, \mu_n]$ vía recurrencia de tres términos $\mathcal{O}(n)$ exacta.
* `Puig_stats(λ, k, T)` / `Puig_stats(params)`: `NamedTuple` `(mean, var, sig, skewness, kurtosis)`.

### Cuantiles, Muestreo y Otras Propiedades
* `Puig_quantile(p, params; tol=1e-10, maxit=200)`: Función cuantil inversa $Q(p)$ mediante bisección.
* `Puig_ci98(params)`: Intervalo central $(P_1, P_{99})$ (cobertura del 98%).
* `Puig_ci996(params)`: Intervalo central $(P_{0.1}, P_{99.9})$ (cobertura del 99.8%).
* `Puig_logpdf(x, params)`: Log-densidad $\log f_P(x)$ numéricamente estable con `SpecialFunctions.besselix`.
* `Puig_hazard(x, params)`: Función de riesgo instantáneo $h(x) = f(x) / S(x)$.
* `Puig_entropy(params; N=400)`: Entropía diferencial teórica $H(X) = -\int f(x) \log f(x) dx$.
* `rand([rng], params, [n])`: Generación de $n$ números pseudoaleatorios exactos.

### Funciones Límite de Maxwell-Boltzmann ($\lambda = 0$)
* `MB_pdf(x, k, B)` / `MB_pdf(x, params::MB)`
* `MB_surviving(x, k, B)` / `MB_surviving(x, params::MB)`
* `MB_cumulative(x, k, B)` / `MB_cumulative(x, params::MB)`

---

## 3. Módulo de Ajuste y Métricas (`puig_fit.jl`)

### `effective_dimension`
```julia
effective_dimension(M::AbstractMatrix{<:Real}) -> Float64
```
Calcula la dimensión efectiva $k = d^2 / \|R\|_F^2 = d^2 / (d + 2 \sum_{i < j} R_{ij}^2)$ en tiempo $\mathcal{O}(d^2)$ con 0 alocaciones.

---

### `empirical_survival`
```julia
empirical_survival(times::AbstractVector{<:Real}, grid::AbstractVector{<:Real}) -> Vector{Float64}
```
Calcula la curva de supervivencia empírica $S_{\text{emp}}(x) = \frac{1}{N} \sum_{i=1}^N \mathbb{I}(X_i \ge x)$.

---

### `Puig_fit`
```julia
# Ajuste Multivariante (d ≥ 2 variables)
Puig_fit(Data::DataFrame; vars=nothing, time_col=nothing, method=:dynamic) -> PuigFitResult
Puig_fit(Data::AbstractMatrix{<:Real}; times=nothing, method=:dynamic) -> PuigFitResult

# Ajuste Univariante 1D (vector de normas / indicador)
Puig_fit(Data::AbstractVector{<:Real}; method=:dynamic, k_fixed=nothing) -> PuigFitResult
```
Ajusta los datos a la distribución de Puig.

**Métodos para datos multivariantes (`DataFrame` / `Matrix`):**
- `:dynamic` (recomendado): $k$ fijado por la correlación ($k_{\text{corr}}$) y $B$ minimizando $L_1$.
- `:fixed`: $k = d$ fijo (hipótesis de independencia) y $B$ minimizando $L_1$.
- `:numerical`: Optimización conjunta $(k, B)$ con L-BFGS-B minimizando el error cuadrático $L_2$.

**Métodos para datos univariantes 1D (`Vector`):**
- `:dynamic` / `:numerical` (defecto): Optimización conjunta 3D de $(\lambda, k, B)$ minimizando el error cuadrático $L_2$ de la supervivencia.
- `:fixed` o `k_fixed !== nothing`: Fija $k$ en el valor nominal especificado y optimiza conjuntamente $(\lambda, B)$.

---

### Métricas Individuales
* `calc_r2(S_emp, S_mod)`: Coeficiente $R^2$.
* `calc_mae(S_emp, S_mod)`: Error absoluto medio.
* `calc_rmse(S_emp, S_mod)`: Raíz del error cuadrático medio.
* `calc_maxae(S_emp, S_mod)`: Error absoluto máximo.
* `calc_iae(x, S_emp, S_mod)`: Error absoluto integrado mediante la regla del trapecio.
* `calc_all_metrics(x, S_emp, S_mod) -> PuigStats{Float64}`: Empaqueta todas las métricas en un struct `PuigStats`.

---

## 4. Visualización (`puig_plot.jl`)

### `Puig_plot`
```julia
Puig_plot(params::Union{PuigDistribution, MB}, options::String="all") -> Plots.Plot
Puig_plot(params::Union{PuigDistribution, MB}, options::Tuple{Vararg{String}}) -> Plots.Plot
Puig_plot(λ::Real, k::Real, T::Real, options="all") -> Plots.Plot
```
Genera figuras visuales de la distribución con `Plots.jl`. Opciones válidas: `"all"`, `"pdf"`, `"cdf"`, `"surv"`, `"data"`, o tuplas combinadas como `("pdf", "surv")`.
