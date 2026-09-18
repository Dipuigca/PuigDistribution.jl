# Contexto y Mapa General del Proyecto — DistributionsPuig.jl

**Fecha de actualización:** 11 de septiembre de 2026\
**Autor del proyecto:** Diego Puig (Instituto Universitario de Matemática Pura y Aplicada — IUAMPA, Universitat Politècnica de València)\
**Estado actual:** Paquete modularizado y limpio (`DistributionsPuig.jl`), 205/205 tests pasando, interfaz completa con `Distributions.jl`. `Puig_pdf` admite keyword `method` para elegir entre convergencia asintótica (`:asymp`, por defecto) y ArbNumerics (`:arb`). La interfaz `Distributions.jl` usa exclusivamente `:asymp`. Expansión asintótica `marcutq_asymp` para colas superiores lejanas. Momentos genéricos de orden $n$ vía recurrencia de tres términos $\mathcal{O}(n)$ exacta.

------------------------------------------------------------------------

## 1. Resumen Ejecutivo y Propósito

El proyecto implementa en Julia la **Distribución de Puig** (`DistributionsPuig.jl`), una generalización continua de la distribución $\chi$ no central con dimensión real continua $k \in \mathbb{R}^+$ ($k \ge 1$).

### 1.1. Problema que resuelve

En el modelado multivariante de rendimiento deportivo y análisis de variables compuestas, es habitual combinar $d$ variables mediante su norma euclídea: $$X = \|\mathbf{Y}\| = \sqrt{\sum_{j=1}^d Y_j^2}$$

Los enfoques tradicionales asumen que las variables componentes son independientes, modelando la norma con una distribución $\chi$ no central estándar con grados de libertad enteros $k = d$. Sin embargo, en disciplinas físicas, biomédicas y psicométricas, las componentes exhiben correlaciones moderadas o altas (por ejemplo, en *powerlifting*: Sentadilla, Press de Banca y Peso Muerto).

Asumir independencia ($k = d$) sobreestima la dimensionalidad efectiva del espacio y subestima las probabilidades de cola. La **Distribución de Puig** resuelve esto estimando una **dimensión efectiva continua** $k \in [1, d]$ a partir de la matriz de correlación empírica, lo que reduce el error absoluto integrado (IAE) en hasta un **42%** frente al modelo tradicional.

------------------------------------------------------------------------

## 2. Mapa y Estructura Actual del Repositorio

El repositorio se encuentra estructurado bajo el estándar idiomático de paquetes de Julia:

```         
00-libreria puigdist/
├── Project.toml              # Definición del paquete (UUID, versión 0.1.0, deps y compat)
├── Manifest.toml             # Árbol de dependencias resuelto (lockfile del entorno Julia)
├── README.md                 # Guía de inicio rápido, fórmulas y documentación oficial
├── CITATION.bib              # Archivo estándar de cita bibliográfica (BibTeX)
├── contexto.md               # [Este documento] Mapa y memoria técnica integral del proyecto
│
├── src/                      # Código fuente del paquete DistributionsPuig
│   ├── DistributionsPuig.jl  # Módulo principal (imports, includes y API pública)
│   ├── Maxwell_Boltzmann.jl  # Caso límite λ → 0 (modelo MB / χ generalizada)
│   ├── puig_pdf.jl           # Núcleo analítico: densidades, supervivencia, momentos y estabilidad
│   ├── types.jl              # Integración formal con Distributions.jl (ContinuousUnivariateDistribution)
│   ├── puig_fit.jl           # Pipeline de ajuste (DataFrame, Matrix, Vector), k Frobenius y PuigStats
│   └── puig_plot.jl          # Visualización desacoplada con Plots.jl
│
├── test/                     # Pruebas automatizadas del paquete
│   └── runtests.jl           # Suite de 100 tests unitarios verificados con Pkg.test()
│
└── docs/                     # Documentación generada con Documenter.jl
    ├── Project.toml          # Dependencias del entorno de documentación
    ├── Manifest.toml         # Lockfile de dependencias de docs
    ├── make.jl               # Script de compilación de la documentación HTML
    ├── src/                  # Documentos fuente en Markdown
    │   ├── index.md          # Portada y visión general
    │   ├── theory.md         # Fundamentación matemática y deducciones analíticas
    │   ├── tutorial.md       # Guía práctica y ejemplos de uso
    │   └── api.md            # Referencia exhaustiva de funciones y tipos
    └── build/                # Sitio estático HTML generado (listo para GitHub Pages)
```

### 2.1. Inventario detallado de archivos actuales

| Archivo / Directorio | Líneas / Tamaño | Rol en el proyecto | Estado |
|:---|:---|:---|:---|
| `Project.toml` | 31 líneas / 859 B | Metadatos del paquete, dependencias directas y rangos `[compat]` | ✅ Vigente |
| `Manifest.toml` | 57.7 KB | Árbol completo de dependencias resueltas en Julia | ✅ Vigente |
| `README.md` | 236 líneas / 9.8 KB | Guía de bienvenida, fórmulas matemáticas, quickstart y badges | ✅ Actualizado |
| `CITATION.bib` | 13 líneas / 498 B | Entrada BibTeX estándar del artículo publicado en *Symmetry* | ✅ Vigente |
| `contexto.md` | Documento maestro | Registro exhaustivo de arquitectura, matemáticas y roadmap | ✅ Actualizado |
| `src/DistributionsPuig.jl` | 46 líneas / 1.7 KB | Punto de entrada del módulo, `using`, `include` y `export` | ✅ Operativo |
| `src/Maxwell_Boltzmann.jl` | 81 líneas / 3.2 KB | Struct `MB` y métodos analíticos para el límite $\lambda \to 0$ | ✅ Operativo |
| `src/puig_pdf.jl` | 557 líneas / 23 KB | `PuigDistribution`, PDF con opción `method` (`:asymp`/`:arb`), Marcum-Q, momentos | ✅ Operativo |
| `src/types.jl` | 72 líneas / 3.6 KB | Subtipo `ContinuousUnivariateDistribution` (`Distributions.jl`) | ✅ Operativo |
| `src/puig_fit.jl` | 405 líneas / 15.2 KB | `Puig_fit`, `effective_dimension` ($\mathcal{O}(d^2)$), `PuigStats`, `PuigFitResult` | ✅ Operativo |
| `src/puig_plot.jl` | 96 líneas / 4.3 KB | Funciones de visualización desacopladas basadas en `Plots.jl` | ✅ Operativo |
| `test/runtests.jl` | 340 líneas / 12.0 KB | 160 pruebas unitarias para la API custom, Distributions.jl y robustez asintótica | ✅ 156/156 pasando |
| `docs/make.jl` | 23 líneas / 624 B | Script de generación de documentación HTML | ✅ Operativo |
| `docs/src/*.md` | 4 archivos / $\sim 19$ KB | Páginas fuente de la documentación (index, theory, tutorial, api) | ✅ Vigentes |
| `docs/build/` | Directorio HTML | Compilación web de la documentación con buscador y temas CSS | ✅ Compilado |

------------------------------------------------------------------------

## 3. Fundamentación Matemática del Núcleo

### 3.1. Función de Densidad de Probabilidad (PDF)

Para una variable aleatoria real no negativa $X \ge 0$: $$f_P(x; \lambda, k, T) = T \frac{x^{k/2}}{\lambda^{k/2 - 1}} \exp\!\left(-\frac{T}{2}(x^2 + \lambda^2)\right) I_{k/2 - 1}(x \lambda T)$$

Donde los tres parámetros son: - $\lambda \ge 0$: Parámetro de no centralidad, correspondiente a la norma euclídea del vector de medias de las disciplinas componentes: $\lambda = \|\boldsymbol{\mu}\| = \sqrt{\sum_{j=1}^d \mu_j^2}$. - $k \ge 1$: Dimensión efectiva continua ($k \in [1, d]$). - $T > 0$: Parámetro de precisión / escala inversa ($T = 1/\sigma^2$). - $I_\nu(z)$: Función de Bessel modificada de primera especie de orden $\nu = k/2 - 1$.

### 3.2. Relación con la distribución $\chi^2$ no central y Muestreo Exacto

Si $\mathbf{Y} \sim \mathcal{N}_k(\boldsymbol{\mu}, \sigma^2 \mathbf{I})$ y $X = \|\mathbf{Y}\|$, se cumple la equivalencia exacta: $$X = \sqrt{\frac{W}{T}}, \quad \text{donde } W \sim \text{NoncentralChisq}(k, \lambda^2 T)$$

Esta formulación permite que `Base.rand` genere muestras pseudoaleatorias exactas sin recurrir a aproximaciones de Monte Carlo:

``` julia
W = rand(rng, NoncentralChisq(k, λ^2 * T))
return sqrt(W / T)
```

### 3.3. Supervivencia y Acumulada nativas (Función Marcum-$Q$)

La función de supervivencia $S(x) = P(X \ge x)$ y la función acumulada $F(x) = P(X \le x)$ se calculan de manera nativa sin dependencias externas en C, R o Fortran mediante la función Marcum-$Q$ generalizada de orden $m = k/2$: $$S(x) = Q_{k/2}\left(\lambda \sqrt{T}, x \sqrt{T}\right) = \mathrm{ccdf}\left(\text{NoncentralChisq}(k, \lambda^2 T), x^2 T\right)$$ $$F(x) = 1 - S(x) = \mathrm{cdf}\left(\text{NoncentralChisq}(k, \lambda^2 T), x^2 T\right)$$

### 3.4. Dimensión efectiva optimizada en $\mathcal{O}(d^2)$ con 0 alocaciones

A partir de la matriz de correlación muestral $R \in \mathbb{R}^{d \times d}$, la dimensión efectiva (razón de participación) se define como: $$k = \frac{\left(\sum_{i=1}^d \lambda_i\right)^2}{\sum_{i=1}^d \lambda_i^2}$$

Aprovechando que la traza de $R$ es constante ($\mathrm{tr}(R) = \sum \lambda_i = d$) y que la suma de los cuadrados de los autovalores es idénticamente la norma de Frobenius al cuadrado $\|R\|_F^2 = d + 2\sum_{i < j} R_{ij}^2$, se obtiene la forma cerrada: $$k = \frac{d^2}{d + 2 \sum_{i < j} R_{ij}^2}$$

**Ventaja computacional:** Evita la descomposición espectral completa (`eigvals` vía LAPACK, que es $\mathcal{O}(d^3)$ y aloca memoria), reduciendo el cómputo a un simple bucle triangular sobre la matriz de correlación en $\mathcal{O}(d^2)$ con 0 alocaciones en el heap.

### 3.5. Momentos Teóricos Analíticos

- **Segundo momento no centrado (forma cerrada):** $$\mu_2 = E[X^2] = \frac{k}{T} + \lambda^2$$
- **Cuarto momento no centrado (forma cerrada):** $$\mu_4 = E[X^4] = \mu_2^2 + \frac{2k}{T^2} + \frac{4\lambda^2}{T}$$
- **Primer momento (Media) y Tercer momento (vía funciones de Laguerre generalizadas):** $$\mu_1 = E[X] = \sqrt{\frac{\pi}{2T}} L_{1/2}^{(k/2 - 1)}\left(-\frac{\lambda^2 T}{2}\right)$$ $$\mu_3 = E[X^3] = \frac{3}{T} \sqrt{\frac{\pi}{2T}} L_{3/2}^{(k/2 - 1)}\left(-\frac{\lambda^2 T}{2}\right)$$

La evaluación de las funciones de Laguerre de orden fraccionario $L_n^{(\alpha)}(z)$ se realiza mediante la función hipergeométrica confluente de Kummer $M(a, b, z)$ implementada con precisión arbitraria en `src/puig_pdf.jl` (`kummer_M_stable`).

A partir de $\mu_1, \dots, \mu_4$ se deducen de forma exacta: $$\mathrm{Var}(X) = \mu_2 - \mu_1^2, \quad \sigma = \sqrt{**\operatorname**{Var}(X)}$$ $$\text{Skewness } (\gamma_1) = \frac{\mu_3 - 3\mu_1\mu_2 + 2\mu_1^3}{\sigma^3}$$ $$\text{Kurtosis } (\gamma_2) = \frac{\mu_4 - 4\mu_1\mu_3 + 6\mu_1^2\mu_2 - 3\mu_1^4}{\sigma^4}$$

### 3.6. Caso Límite $\lambda \to 0$ (Maxwell-Boltzmann / $\chi$ Generalizada)

Cuando $\lambda = 0$, la densidad colapsa continuamente al modelo de Maxwell-Boltzmann con parámetro de dispersión térmica $B = 1/T$: $$f(x; k, B) = \frac{2^{1 - k/2} B^{-k/2}}{\Gamma(k/2)} x^{k-1} \exp\!\left(-\frac{x^2}{2B}\right)$$ $$S(x) = \frac{\Gamma(k/2, x^2 / (2B))}{\Gamma(k/2)}$$

### 3.7. Justificación de las restricciones $k \ge 1$, $\lambda \ge 0$, $T > 0$

- $\lambda \ge 0$: Representa la norma de un vector de medias ($\|\boldsymbol{\mu}\| \ge 0$). Para $\lambda < 0$, el argumento $x\lambda T$ se volvería negativo y $I_\nu(z)$ adquiriría valores complejos para $\nu$ no entero.
- $k \ge 1$: Número mínimo de dimensiones del espacio. Para $k < 1$, el factor $x^{k/2} I_{k/2-1}(x\lambda T) \sim x^{k-1}$ diverge en el origen impidiendo la normalización de la probabilidad.
- $T > 0$: Inversa de varianza ($T = 1/\sigma^2$). $T \le 0$ produciría una integral de densidad divergente.

### 3.8. Estabilidad Numérica y Método Dual (`:asymp` / `:arb`)

Cuando $\lambda > 30$ o $\lambda > 100$ (típico en puntuaciones totales de powerlifting donde $\lambda \approx 100-300$), se produce una cancelación catastrófica extrema en `Float64`: - $\exp(-T/2 \cdot (x^2 + \lambda^2))$ genera subdesbordamiento (*underflow*) a `0.0`. - $I_\nu(x\lambda T)$ genera desbordamiento (*overflow*) a `Inf`. - Su producto directo produce `0.0 * Inf = NaN`.

**Solución dual con elección de método:**

La función `Puig_pdf` acepta un keyword `method::Symbol` que permite al usuario elegir entre dos estrategias de cálculo:

``` julia
Puig_pdf(x, λ, k, T; method=:asymp)   # Por defecto: expansión asintótica (recomendado)
Puig_pdf(x, λ, k, T; method=:arb)     # ArbNumerics 128-bit (referencia de alta precisión)
```

1.  **`:asymp`** *(por defecto)*: Expansión asintótica de 5 términos (DLMF 10.40.1) para $z = x\lambda T \ge 200$, y `SpecialFunctions.besselix` para $z < 200$. Ambas absorben el decaimiento exponencial en $\exp(-T/2 \cdot (x-\lambda)^2)$, evitando desbordamientos. Precisión relativa $< 10^{-12}$. Método ultrarrápido en `Float64` puro.

2.  **`:arb`**: Utiliza `ArbNumerics` fijado a **128 bits** ($\sim 38$ dígitos significativos, exponentes de hasta $\pm 10^{10^9}$). Ambos factores se calculan con rango suficiente y su producto es analíticamente exacto. Útil como referencia de validación.

**Nota**: La interfaz `Distributions.jl` (`pdf(dist, x)`) utiliza exclusivamente el método `:asymp` para garantizar compatibilidad y rendimiento óptimo.

**`Puig_logpdf`**: Opera en `Float64` de forma ultrarrápida empleando la función de Bessel exponencialmente escalada \`SpecialFunctions.besselix(\nu, z) = I\_\nu(z)e\^{-z}\$, absorbiendo el término cuadrático: $\exp(-T/2 \cdot (x-\lambda)^2)$.

------------------------------------------------------------------------

## 4. Arquitectura de Código y API del Paquete

### 4.1. Estructuras de Datos Principales

``` julia
# Distribución de Puig (subtipo de Distributions.jl)
struct PuigDistribution{S<:Real} <: ContinuousUnivariateDistribution
    λ::S  # Norma de medias (λ ≥ 0)
    k::S  # Dimensión continua (k ≥ 1)
    T::S  # Precisión / Escala inversa (T > 0)
end

# Caso límite λ = 0
struct MB{S<:Real}
    k::S  # Dimensión real (k ≥ 1)
    B::S  # Parámetro térmico (B = 1/T > 0)
end

# Métricas de ajuste
struct PuigStats{S<:Real}
    R²::S; MAE::S; RMSE::S; MaxAE::S; IAE::S
end

# Contenedor de resultado con desempaquetado directo
struct PuigFitResult{S<:Real}
    params::PuigDistribution{S}
    indicator::Vector{S}
    stats::PuigStats{S}
end
```

### 4.2. Desempaquetado Idiomático

`PuigFitResult` implementa `Base.iterate`, permitiendo:

``` julia
dist, ind, stats = Puig_fit(df)
println("k estimado: ", dist.k)
println("R² ajuste:  ", stats.R²)
```

### 4.3. Interfaz Dual de Métodos

| Funcionalidad | API Custom Vectorizada | Interfaz Distributions.jl |
|:---|:---|:---|
| Densidad (PDF) | `Puig_pdf(x, dist; method=:asymp)` | `pdf(dist, x)` ← siempre `:asymp` |
| Log-densidad | `Puig_logpdf(x, dist)` | `logpdf(dist, x)` |
| Acumulada (CDF) | `Puig_cumulative(x, dist)` | `cdf(dist, x)` |
| Supervivencia (CCDF) | `Puig_surviving(x, dist)` | `ccdf(dist, x)` |
| Cuantil | `Puig_quantile(p, dist)` | `quantile(dist, p)` |
| Muestreo | `rand(dist, n)` | `rand(rng, dist, n)` |
| Media | `Puig_mean(dist)` | `mean(dist)` |
| Varianza | `Puig_var(dist)` | `var(dist)` |
| Desviación típica | `Puig_std(dist)` | `std(dist)` |
| Asimetría (Skewness) | `Puig_skewness(dist)` | `skewness(dist)` |
| Kurtosis | `Puig_kurtosis(dist)` | `kurtosis(dist)` |
| Parámetros | `(dist.λ, dist.k, dist.T)` | `params(dist)` |
| Soporte | `[0.0, Inf)` | `insupport(dist, x)`, `minimum`, `maximum` |

### 4.4. Pipeline de Ajuste (`Puig_fit`)

La función `Puig_fit` admite tres tipos de entrada de datos: 1. **`DataFrame`**: Detecta automáticamente las columnas numéricas (o las pasadas en `vars`), computa las normas por fila (o toma `time_col`), calcula $\lambda$ analíticamente como $\|\bar{\mathbf{X}}\|$ y ejecuta el ajuste. 2. **`Matrix`**: Filas = observaciones, columnas = variables componentes. 3. **`Vector`**: Serie 1D de normas observadas. Realiza optimización continua 3D $(\log \lambda, \log(k-1), \log B)$ mediante Nelder-Mead de `Optim.jl`. Si se especifica `k_fixed`, fija la dimensión nominal y optimiza conjuntamente $(\lambda, B)$.

**Métodos de ajuste disponibles:** - `:dynamic` *(recomendado)*: $k$ fijado analíticamente por la correlación ($k_{\text{corr}}$) y $B$ optimizado minimizando la distancia $L_1$ sobre la curva de supervivencia. - `:fixed`: $k = d$ fijo (supuesto de independencia clásica) y $B$ minimizando $L_1$. - `:numerical`: Optimización conjunta $(k, B)$ con restricciones de caja mediante L-BFGS-B minimizando el error cuadrático $L_2$.

### 4.5. Visualización (`Puig_plot`)

Módulo desacoplado en `src/puig_plot.jl` basado en `Plots.jl`: - `Puig_plot(dist, "all")`: Panel 2x2 completo con PDF, CDF, Supervivencia y Ficha de Estadísticos/Cuantiles. - `Puig_plot(dist, "pdf")`, `Puig_plot(dist, "cdf")`, `Puig_plot(dist, "surv")`, `Puig_plot(dist, "data")`. - `Puig_plot(dist, ("pdf", "surv"))`: Composición matricial de las vistas seleccionadas.

------------------------------------------------------------------------

## 5. Estado de Verificación y Tests (100/100 Pasando)

La suite de pruebas en `test/runtests.jl` cubre exhaustivamente tanto el núcleo matemático como la interfaz del ecosistema:

```         
Suite de Tests: DistributionsPuig.jl (156 tests pasando en ~12s)
├── 0. Módulo y tipos (PuigDistribution, MB)
├── 1. effective_dimension (1D, ortogonal, colineal, fórmula Frobenius)
├── 2. empirical_survival (condiciones de frontera y monotonía)
├── 3. Métricas y struct PuigStats (ajuste perfecto y errores controlados)
├── 4. Puig_fit sobre DataFrame y Matrix (:dynamic, :fixed, :numerical, vars, time_col, 1D con/sin k_fixed)
├── 5. Desempaquetado e Indexación de PuigFitResult (dist, ind, st = Puig_fit(...))
├── 6. Validación de Errores y Casos Límite (N < 10, métodos desconocidos, k_fixed < 1, DataFrames sin columnas numéricas)
├── 7. Interfaz Distributions.jl (pdf, logpdf, cdf, ccdf, quantile, mean, var, std, skewness, kurtosis, rand, insupport, params)
└── 8. Robustez asintótica vs ArbNumerics (λ extremos, k extremos, colas, coherencia :asymp/:arb)
```

Comprobación ejecutada con éxito mediante:

``` bash
julia --project=. -e "using Pkg; Pkg.test()"
```

------------------------------------------------------------------------

## 6. Contexto Científico y Metodológico (Origen del Proyecto)

El paquete tiene su base en la línea de investigación desarrollada por **Diego Puig** y colaboradores en el **IUAMPA (Universitat Politècnica de València)**.

### 6.1. Artículo original publicado (*Symmetry*, 2025)

El fundamento metodológico inicial del uso de la distribución $\chi$ no central en evaluación deportiva se encuentra publicado en la revista *Symmetry* (MDPI):

> **Puig Castro, D., Coronado Ferrer, A., Castro Palacio, J. C., Fernández de Córdoba, P., Ortigosa, N., & Sánchez Pérez, E. A. (2025).**\
> *Non-Centered Chi Distributions as Models for Fair Assessment in Sports Performance*.\
> **Symmetry**, 17(7), 1039.\
> DOI: [10.3390/sym17071039](https://doi.org/10.3390/sym17071039)\
> URL: <https://www.mdpi.com/2073-8994/17/7/1039>

En este trabajo fundacional se demuestra la idoneidad de la distribución $\chi$ no central para modelar variables compuestas asimétricas en rendimiento deportivo, estableciendo percentiles equitativos de clasificación frente a la dispersión estocástica.

### 6.2. Generalización a dimensión real continua ($k \in \mathbb{R}^+$)

El paquete `DistributionsPuig.jl` expande dicho trabajo resolviendo la limitación de la independencia entre disciplinas: - **Estudio extendido**: *"Generalized Non-Central Chi Distribution with Real-Valued Dimensionality for Powerlifting Performance Classification"*. - **Dataset**: 11.382 levantadores competitivos de la *International Powerlifting Federation* (IPF, temporadas 2023–2024) en 16 categorías de peso corporal (8 masculinas, 8 femeninas). - **Disciplinas evaluadas**: Sentadilla (*Squat*), Press de Banca (*Bench Press*) y Peso Muerto (*Deadlift*). - **Indicadores analizados**: - **DPR** (*Discipline/Per-Weight Ratio*): Levantamientos normalizados por el peso corporal de cada atleta. - **DP** (*Discipline/Mean Ratio*): Levantamientos normalizados por la media de su categoría de peso. - **DPNM** (*Discipline Performance Norm in Kilograms*): Norma euclídea cruda en kilogramos absolutos. - **Conclusiones clave**: - En la norma cruda en kilogramos (DPNM), la dimensión efectiva $k$ es muy cercana a $3.0$ (el espacio de fuerza absoluta es casi isotrópico). - En los indicadores normalizados (DPR y DP), las componentes presentan correlaciones significativas ($\rho \approx 0.7-0.9$), lo que reduce la dimensión efectiva a valores $k \approx 1.8 - 2.4$. - Asumir $k = 3$ fijo comete errores graves de sobreestimación en las colas. La distribución de Puig con dimensión dinámica reduce el error absoluto integrado (IAE) en hasta un **42%**.

------------------------------------------------------------------------

## 7. Diagnóstico Técnico y Tareas Pendientes

### 7.1. Ajuste y Compilación de `docs/make.jl` ✅ (Completado)

El script `docs/make.jl` ha sido actualizado para cargar el paquete directamente desde el entorno raíz:

``` julia
pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter
using DistributionsPuig
```

La documentación HTML fue compilada con éxito mediante `Documenter.jl` en `docs/build/` con código de salida 0.

### 7.2. Generación del nuevo `README.md` Oficial ✅ (Completado)

Se ha creado el archivo `README.md` oficial del paquete en la raíz, incluyendo: - Badges de versión de Julia (1.10+), licencia MIT, tests (100/100 passing) y documentación. - Resumen matemático riguroso (PDF, relación con $\chi^2$ no central, Marcum-$Q$, Frobenius $\mathcal{O}(d^2)$ y estabilidad dual). - Guía de inicio rápido para `Distributions.jl` y la API custom vectorizada. - Tutorial de ajuste multivariante con `DataFrames` y desempaquetado directo (`dist, ind, st = Puig_fit(df)`). - Ejemplos de visualización con `Puig_plot`. - Cita académica del trabajo de Diego Puig (IUAMPA - UPV).

### 7.3. Integración Continua (CI) 📅 (Pendiente)

Añadir el flujo de trabajo `.github/workflows/CI.yml` para validar automáticamente los tests en GitHub Actions sobre Julia 1.10, 1.11 y 1.12 en Windows, Linux y macOS.

### 7.4. Publicación en el Julia General Registry 📅 (Pendiente)

Una vez configurado CI y revisada la documentación, registrar `DistributionsPuig` en el registro general de Julia mediante el bot `@JuliaRegistrator`.

------------------------------------------------------------------------

## 8. Guía Rápida para Desarrolladores

### Activar el entorno y ejecutar los tests

``` julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

### Uso básico del paquete

``` julia
using DistributionsPuig
using Distributions

# 1. Definir una distribución
d = PuigDistribution(3.5, 2.4, 0.25)

# 2. Funciones básicas
x = 2.0
f = pdf(d, x)
S = ccdf(d, x)
q = quantile(d, 0.95)
muestras = rand(d, 1000)

# 3. Ajuste de un DataFrame
using DataFrames
df = DataFrame(Squat=rand(100), Bench=rand(100), Deadlift=rand(100))
modelo, ind, stats = Puig_fit(df; method=:dynamic)
```
