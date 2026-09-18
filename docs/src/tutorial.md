# Tutorial y Ejemplos de Uso

En esta guía práctica se ilustra el flujo de trabajo completo utilizando `DistributionsPuig.jl`: desde la creación manual de distribuciones hasta el ajuste de datos reales y su visualización.

---

## 1. Carga del Entorno

```julia
using Pkg
Pkg.activate(".")

include("puig_pdf.jl")
include("puig_fit.jl")
include("puig_plot.jl")

using DataFrames
using Random
using Plots
```

---

## 2. Definición y Evaluación Básica

Para definir una distribución con parámetros $\lambda = 3.0$, $k = 2.4$ y $T = 0.8$:

```julia
dist = PuigDistribution(3.0, 2.4, 0.8)

# Rejilla de evaluación
x = 0.0:0.05:6.0

# 1. Densidad de probabilidad (PDF)
f_vals = Puig_pdf(x, dist)

# 2. Distribución acumulada (CDF)
F_vals = Puig_cumulative(x, dist)

# 3. Función de supervivencia (CCDF)
S_vals = Puig_surviving(x, dist)

# 4. Tasa de fallo o riesgo instantáneo (Hazard)
h_vals = Puig_hazard(x, dist)
```

### Estadísticos Teóricos y Cuantiles

```julia
# Momentos y estadísticos descriptivos teóricos
st = Puig_stats(dist)
println("Media teórica: ", round(st.mean, digits=4))
println("Varianza:      ", round(st.var, digits=4))
println("Asimetría:     ", round(st.skewness, digits=4))
println("Kurtosis:      ", round(st.kurtosis, digits=4))

# Cuantiles específicos
q_inter = Puig_quantile([0.25, 0.50, 0.75], dist)
p1, p99 = Puig_ci98(dist)

# Generación de 10,000 muestras aleatorias
muestras = rand(dist, 10000)
```

---

## 3. Ajuste de Datos Multivariantes (`Puig_fit`)

Supongamos que disponemos de un `DataFrame` con el rendimiento de atletas en tres disciplinas correlacionadas (por ejemplo, Sentadilla, Press de Banca y Peso Muerto):

```julia
Random.seed!(42)
N = 1000

# Simulación de un factor latente común + ruido individual
factor = randn(N)
squat = 180.0 .+ 25.0 .* (0.8 .* factor .+ 0.6 .* randn(N))
bench = 120.0 .+ 18.0 .* (0.75 .* factor .+ 0.66 .* randn(N))
deadlift = 220.0 .+ 30.0 .* (0.85 .* factor .+ 0.53 .* randn(N))

df_powerlifting = DataFrame(
    Sentadilla = squat,
    Banca      = bench,
    PesoMuerto = deadlift
)
```

### Ejecución del Ajuste

```julia
# Ajuste usando el método dinámico recomendado
resultado = Puig_fit(df_powerlifting; method=:dynamic)

# Desempaquetado directo
modelo, indicador_obs, metricas = resultado

println("--- Parámetros Estimados ---")
println("λ estimado: ", round(modelo.λ, digits=3))
println("k efectivo: ", round(modelo.k, digits=3))
println("T estimado: ", round(modelo.T, digits=6))

println("\n--- Métricas de Bondad de Ajuste ---")
println("R²:    ", round(metricas.R², digits=5))
println("MAE:   ", round(metricas.MAE, digits=5))
println("RMSE:  ", round(metricas.RMSE, digits=5))
println("IAE:   ", round(metricas.IAE, digits=5))
```

---

## 4. Comparación de Estrategias de Ajuste

Podemos comparar el rendimiento de los tres métodos disponibles:

```julia
# 1. Dinámico: k estimado por correlación (Frobenius), B por L1
res_dyn = Puig_fit(df_powerlifting; method=:dynamic)

# 2. Fijo: k = d = 3 (hipótesis de independencia clásica), B por L1
res_fix = Puig_fit(df_powerlifting; method=:fixed)

# 3. Numérico: optimización conjunta de (k, B) mediante L-BFGS-B (L2)
res_num = Puig_fit(df_powerlifting; method=:numerical)

println("R² Dinámico (k=$(round(res_dyn.params.k, digits=2))):  $(round(res_dyn.stats.R², digits=5))")
println("R² Fijo     (k=3.00): $(round(res_fix.stats.R², digits=5))")
println("R² Numérico (k=$(round(res_num.params.k, digits=2))):  $(round(res_num.stats.R², digits=5))")
```

El modelo dinámico proporciona un $R^2$ significativamente superior al modelo con $k=3$ fijo, confirmando el impacto de capturar la correlación entre disciplinas.

---

## 5. Visualización de Resultados (`puig_plot.jl`)

### Panel Completo 2x2
```julia
# Genera panel con PDF, CDF, Supervivencia y Estadísticos
fig_completa = Puig_plot(modelo, "all")
savefig(fig_completa, "ajuste_puig_panel.png")
```

### Gráficos Individuales o Combinados
```julia
# Solo función de supervivencia
fig_surv = Puig_plot(modelo, "surv")

# Combinación personalizada: Densidad y Supervivencia
fig_combo = Puig_plot(modelo, ("pdf", "surv"))
```

