# DistributionsPuig.jl

[![Julia](https://img.shields.io/badge/Julia-1.10%2B-blue.svg)](https://julialang.org) [![Licence](https://img.shields.io/badge/Licence-MIT-green.svg)](LICENSE) [![Tests](https://img.shields.io/badge/Tests-205%2F205%20passing-brightgreen.svg)](test/runtests.jl) [![Documentation](https://img.shields.io/badge/docs-Documenter.jl-blue.svg)](docs/)

Implementación en Julia de la **Distribución de Puig** (distribución $\chi$ no central generalizada con parámetro de dimensión real continuo $k \in \mathbb{R}^+$, $k \ge 1$), diseñada para el modelado y clasificación de rendimiento multivariante, análisis de variables compuestas y normas de vectores correlacionados.

------------------------------------------------------------------------

## 📖 Fundamentación Matemática

La distribución de Puig generaliza la distribución $\chi$ no central permitiendo que el número de grados de libertad o dimensión efectiva $k$ tome valores continuos reales ($k \ge 1$), capturando la correlación y dependencia estadística entre variables componentes.

### 1. Función de Densidad de Probabilidad (PDF)

Para una variable aleatoria real no negativa $X \ge 0$:

$$f_P(x; \lambda, k, T) = T \frac{x^{k/2}}{\lambda^{k/2 - 1}} \exp\!\left(-\frac{T}{2}(x^2 + \lambda^2)\right) I_{k/2 - 1}(x \lambda T)$$

donde: - $\lambda \ge 0$: Norma euclídea del vector de medias de las componentes ($\lambda = \|\boldsymbol{\mu}\| = \sqrt{\sum \mu_j^2}$). - $k \ge 1$: Dimensión efectiva continua ($k \in [1, d]$). - $T > 0$: Parámetro de precisión / escala inversa ($T = 1/\sigma^2$). - $I_\nu(z)$: Función de Bessel modificada de primera especie de orden $\nu = k/2 - 1$.

### 2. Relación con la distribución $\chi^2$ no central y Muestreo Exacto

Si $\mathbf{Y} \sim \mathcal{N}_k(\boldsymbol{\mu}, \sigma^2 \mathbf{I})$ y $X = \|\mathbf{Y}\|$, entonces:

$$X = \sqrt{\frac{W}{T}}, \quad \text{donde } W \sim \text{NoncentralChisq}(k, \lambda^2 T)$$

Esta relación fundamenta el **muestreo pseudoaleatorio exacto** (`rand`) sin aproximaciones numéricas.

### 3. Supervivencia y Acumulada nativas (Función Marcum-$Q$)

La función de supervivencia $S(x) = P(X \ge x)$ y la acumulada $F(x) = P(X \le x)$ se calculan de forma nativa mediante la función Marcum-$Q$ generalizada de orden $m = k/2$:

$$S(x) = Q_{k/2}\left(\lambda \sqrt{T}, x \sqrt{T}\right) = \operatorname{ccdf}\left(\text{NoncentralChisq}(k, \lambda^2 T), x^2 T\right)$$ $$F(x) = 1 - S(x) = \operatorname{cdf}\left(\text{NoncentralChisq}(k, \lambda^2 T), x^2 T\right)$$

Para la **cola superior lejana** ($x \gg \lambda$, es decir $b - a \ge 4$ y $ab \ge 30$), `marcumq` conmuta automáticamente a la expansión asintótica `marcumq_asymp`, que elimina cancelaciones catastróficas y falsos suelos de ruido numérico en probabilidades extremas ($\sim 10^{-50}$).

### 4. Dimensión Efectiva Optimizada ($\mathcal{O}(d^2)$, 0 alocaciones)

A partir de la matriz de correlación $R$ de $d$ disciplinas, la dimensión efectiva (razón de participación) se evalúa analíticamente mediante la norma de Frobenius:

$$k = \frac{\left(\sum_{i=1}^d \lambda_i\right)^2}{\sum_{i=1}^d \lambda_i^2} = \frac{d^2}{\|R\|_F^2} = \frac{d^2}{d + 2 \sum_{i < j} R_{ij}^2}$$

Esta identidad evita calcular autovalores vía LAPACK ($\mathcal{O}(d^3)$), reduciendo la complejidad a $\mathcal{O}(d^2)$ con **cero alocaciones de memoria**.

### 5. Estabilidad Numérica Dual (`:asymp` / `:arb`)

Para valores grandes de $\lambda$ ($\lambda > 30$), el producto directo $\exp \times I_\nu$ produce $0.0 \times \infty = \text{NaN}$ en `Float64`. La función `Puig_pdf` ofrece dos métodos seleccionables mediante el keyword `method`:

``` julia
Puig_pdf(x, λ, k, T; method=:asymp)   # Por defecto: expansión asintótica 5 términos (DLMF 10.40.1)
Puig_pdf(x, λ, k, T; method=:arb)     # ArbNumerics 128-bit (referencia de alta precisión)
```

- **`:asymp`** *(por defecto)*: Para $z = x\lambda T \ge 200$ usa expansión asintótica; para $z < 200$ usa `SpecialFunctions.besselix`. Precisión relativa $< 10^{-12}$, ultrarrápido en `Float64` puro.
- **`:arb`**: `ArbNumerics` con precisión de 128 bits ($\sim 38$ dígitos). Referencia de validación.

> **Nota**: `Distributions.pdf(dist, x)` utiliza siempre el método `:asymp`.

------------------------------------------------------------------------

## 🚀 Instalación y Carga

``` julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using DistributionsPuig
using Distributions
```

------------------------------------------------------------------------

## 💡 Ejemplos de Uso

### 1. Interfaz Nativa con `Distributions.jl`

`PuigDistribution` es un subtipo formal de `ContinuousUnivariateDistribution`, por lo que es 100% compatible con las funciones del ecosistema `Distributions.jl`:

``` julia
using DistributionsPuig
using Distributions

# Crear una distribución: λ = 3.5, k = 2.4, T = 0.25
d = PuigDistribution(3.5, 2.4, 0.25)

# Evaluación
pdf(d, 2.0)            # Densidad de probabilidad (usa :asymp internamente)
logpdf(d, 2.0)         # Log-densidad
cdf(d, 2.0)            # Función acumulada
ccdf(d, 2.0)           # Supervivencia P(X ≥ x)
quantile(d, 0.95)      # Cuantil 95%

# API custom con elección de método
Puig_pdf(2.0, d)                    # asymp (default)
Puig_pdf(2.0, d; method=:arb)      # ArbNumerics 128-bit

# Momentos y estadísticos teóricos exactos
mean(d)                # Media (vía Laguerre generalizado)
var(d)                 # Varianza
std(d)                 # Desviación típica
skewness(d)            # Asimetría
kurtosis(d)            # Kurtosis

# Momentos raw de orden genérico (hasta n=8, n=20, etc.)
moments_4 = Puig_moments(d)      # Tupla (μ₁, μ₂, μ₃, μ₄) — retrocompatible
moments_8 = Puig_moments(d, 8)   # Vector [μ₁, …, μ₈] vía recurrencia O(n) exacta

# Muestreo exacto
muestras = rand(d, 1000)
```

------------------------------------------------------------------------

### 2. Ajuste de Modelos (`Puig_fit`)

La función `Puig_fit` realiza la estimación automática de parámetros sobre tres formatos de datos:

#### A. Desde `DataFrame` (Multivariante)

``` julia
using DataFrames

df = DataFrame(
    Sentadilla = [180.0, 195.0, 160.0, 210.0, 175.0, 190.0, 205.0, 185.0, 170.0, 200.0],
    Banca      = [120.0, 135.0, 110.0, 145.0, 115.0, 130.0, 140.0, 125.0, 118.0, 138.0],
    PesoMuerto = [220.0, 240.0, 200.0, 260.0, 215.0, 235.0, 255.0, 230.0, 210.0, 250.0]
)

# Ajuste dinámico (recomendado)
resultado = Puig_fit(df; method=:dynamic)

# Desempaquetado directo idiomático
distribucion, indicador_obs, metricas = resultado

println("Parámetros: λ=$(distribucion.λ), k=$(distribucion.k), T=$(distribucion.T)")
println("R² de supervivencia: ", metricas.R²)
println("Error absoluto medio (MAE): ", metricas.MAE)
```

#### B. Desde `Matrix`

``` julia
mat = Matrix{Float64}(df)
resultado = Puig_fit(mat; method=:dynamic)
```

#### C. Desde `Vector` 1D (Optimización conjunta 3D)

Cuando solo se dispone de la serie univariante de normas observadas:

``` julia
normas = [285.3, 310.2, 270.8, 335.1, 290.4, 305.7, 328.0, 300.5, 275.9, 320.1]

# Optimización conjunta 3D de (λ, k, B)
res_1d = Puig_fit(normas)

# O fijando la dimensión nominal conocida (ej. k = 3 disciplinas)
res_k3 = Puig_fit(normas; k_fixed=3.0)
```

#### Métodos de ajuste disponibles

| Método (`method`) | Descripción | Optimización |
|:---|:---|:---|
| `:dynamic` *(defecto)* | $k$ fijado analíticamente por correlación ($k_{\text{corr}}$). | $B$ minimizando distancia $L_1$. |
| `:fixed` | $k = d$ fijado a la dimensión geométrica (independencia). | $B$ minimizando distancia $L_1$. |
| `:numerical` | Optimización conjunta de $(k, B)$ con restricciones de caja. | Algoritmo L-BFGS-B minimizando $L_2$. |

------------------------------------------------------------------------

### 3. Visualización Integrada (`Puig_plot`)

El módulo `src/puig_plot.jl` proporciona gráficos desacoplados basados en `Plots.jl`:

``` julia
using Plots

# Panel 2x2 completo (PDF, CDF, Supervivencia y Ficha de Estadísticos)
Puig_plot(distribucion, "all")

# Gráficos específicos
Puig_plot(distribucion, "surv")

# Combinaciones personalizadas
Puig_plot(distribucion, ("pdf", "surv"))
```

------------------------------------------------------------------------

## 🧪 Pruebas Automatizadas

El paquete incluye una suite completa de 205 pruebas unitarias con cobertura del núcleo matemático, algoritmos de ajuste, casos límite, interfaz de `Distributions.jl`, robustez asintótica vs ArbNumerics, momentos genéricos de orden $n$ y la expansión asintótica de Marcum-$Q$ en colas extremas:

``` bash
julia --project=. -e "using Pkg; Pkg.test()"
```

------------------------------------------------------------------------

## 📚 Estructura del Proyecto

```         
00-libreria puigdist/
├── Project.toml              # Metadatos del paquete y dependencias
├── Manifest.toml             # Árbol de dependencias resuelto
├── README.md                 # Esta guía de uso y documentación rápida
├── contexto.md               # Memoria técnica completa y roadmap del proyecto
├── src/
│   ├── DistributionsPuig.jl  # Módulo principal
│   ├── Maxwell_Boltzmann.jl  # Caso límite λ → 0 (modelo MB)
│   ├── puig_pdf.jl           # Núcleo analítico: PDF con opción method (:asymp/:arb), estabilidad numérica
│   ├── types.jl              # Conformidad con Distributions.jl
│   ├── puig_fit.jl           # Algoritmos de ajuste y dimensión efectiva
│   └── puig_plot.jl          # Visualización con Plots.jl
├── test/
│   └── runtests.jl           # 205 tests unitarios automatizados
└── docs/                     # Documentación con Documenter.jl
    ├── make.jl               # Compilador de la web estática
    └── src/                  # Fuentes Markdown (index, theory, tutorial, api)
```

------------------------------------------------------------------------

## 📑 Referencia Académica

Si utilizas esta librería o la metodología en tu investigación, por favor cita el artículo original publicado en *Symmetry* (MDPI):

> **Puig Castro, D., Coronado Ferrer, A., Castro Palacio, J. C., Fernández de Córdoba, P., Ortigosa, N., & Sánchez Pérez, E. A. (2025).** *Non-Centered Chi Distributions as Models for Fair Assessment in Sports Performance*. **Symmetry**, 17(7), 1039. <https://doi.org/10.3390/sym17071039>

### BibTeX

``` bibtex
@article{puigcastro2025noncentered,
  author  = {Puig Castro, Diego and Coronado Ferrer, Ana and Castro Palacio, Juan Carlos and Fern{\'a}ndez de C{\'o}rdoba, Pedro and Ortigosa, Nuria and S{\'a}nchez P{\'e}rez, Enrique A.},
  title   = {Non-Centered Chi Distributions as Models for Fair Assessment in Sports Performance},
  journal = {Symmetry},
  year    = {2025},
  volume  = {17},
  number  = {7},
  pages   = {1039},
  doi     = {10.3390/sym17071039},
  url     = {https://www.mdpi.com/2073-8994/17/7/1039}
}
```

Para la generalización continua a dimensión real $k \in \mathbb{R}^+$ y análisis de dependencia dimensional:

> **Puig, D.** *Generalized Non-Central Chi Distribution with Real-Valued Dimensionality for Powerlifting Performance Classification*. Instituto Universitario de Matemática Pura y Aplicada (IUAMPA), Universitat Politècnica de València.

------------------------------------------------------------------------

## 📄 Licencia

Este proyecto está bajo la Licencia MIT.
