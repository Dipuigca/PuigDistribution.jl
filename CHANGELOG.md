# Changelog

Todas las modificaciones relevantes del mono-repo se registran aquí.
Formato basado en [Keep a Changelog](https://keepachangelog.com/es/1.1.0/).
Versionado siguiendo [SemVer](https://semver.org/lang/es/).

## [No publicado]

### Añadido

- Documentación completa de funciones en los portes R y Python:
  - **R**: roxygen `@param`/`@return`/`@export` en todas las funciones
    exportadas y generación de los `man/*.Rd` con roxygen2 (44 temas de ayuda).
  - **Python**: docstrings en los métodos públicos de `PuigDistribution`, `MB`,
    los hooks scipy de `puig_gen` y las funciones de ajuste (`calc_*`,
    `Puig_fit`) y sus helpers internos.
- `CHANGELOG.md` (este fichero).

### Cambiado

- **R**: el paquete se reestructura a la disposición estándar de CRAN: el
  código fuente pasa de la raíz `R/*.R` a `R/R/*.R`. Los tests
  (`R/tests/*.R`) ahora hacen `source` desde `R/R/`.
- **R**: corregido `importFrom` en `NAMESPACE` (se importaban desde `stats`
  funciones base como `gamma`, `lgamma`, `besselI`, `mean`, lo que impedía
  `pkgload::load_all`).
- **R**: añadido wrapper público `kummer_M_stable` (estaba exportado en
  `NAMESPACE` pero no definido).

## [0.1.0] — 2026-09-17

### Añadido

- Distribución de Puig (χ no central generalizada con dimensión k continua):
  PDF, log-PDF, supervivencia, acumulada, cuantiles, momentos, estadísticos,
  función de riesgo, entropía y muestreo.
- Implementación de referencia en Julia (`src/`) con alta precisión
  (ArbNumerics) para la función de Kummer y las funciones de Laguerre
  generalizadas de orden fraccionario.
- Caso límite Maxwell-Boltzmann (λ = 0) con densidad de chi generalizada
  de dimensión real k.
- Ajuste de datos a la distribución de Puig (`Puig_fit`) con métodos
  `:dynamic`, `:fixed` y `:numerical`, dimensión efectiva y métricas de
  bondad de ajuste (R², MAE, RMSE, MaxAE, IAE).
- Generador de golden files (`golden/generate.jl`) y validación cruzada
  entre lenguajes.
- Documentación técnica del bache del kernel (`laguerre_real`) en los momentos
  con λ grande.