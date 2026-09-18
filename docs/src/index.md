# DistributionsPuig.jl

Bienvenido a la documentación oficial de **DistributionsPuig.jl**, una librería en Julia para el análisis, evaluación y ajuste de la **Distribución de Puig** (distribución $\chi$ no central generalizada con dimensión real $k \in \mathbb{R}^+$).

---

## 🎯 Motivación y Contexto

En el modelado estadístico de rendimiento deportivo y análisis multivariante, frecuentemente se combinan múltiples variables continuas en un indicador compuesto mediante su norma euclídea $\|\mathbf{X}\|$. 

Los enfoques tradicionales asumen que las distintas disciplinas o componentes son estadísticamente independientes, modelando la norma resultante con una distribución $\chi$ no central estándar con $d$ grados de libertad enteros ($k = d$).

Sin embargo, en datos reales (por ejemplo, sentadilla, press de banca y peso muerto en *powerlifting*), las componentes presentan correlaciones significativas. La **Distribución de Puig** resuelve este problema estimando una **dimensión efectiva continua** $k \in [1, d]$, capturando fielmente la redundancia de información y reduciendo drásticamente el error de ajuste empírico.

---

## ⚡ Características Principales

1. **Núcleo Matemático Exacto**:
   - Densidad de probabilidad ($f_P$), función de distribución acumulada ($F_P$) y función de supervivencia ($S_P$).
   - Implementación nativa de la función Marcum-$Q$ sin dependencias de entornos externos (como R o C).
   - Momentos teóricos exactos ($\mu_1, \dots, \mu_4$), varianza, asimetría (*skewness*) y *kurtosis* mediante funciones de Laguerre generalizadas con transformación de Kummer.

2. **Pipeline de Ajuste Robusto (`Puig_fit`)**:
   - Estimación ultra-rápida de la dimensión efectiva $k$ en $\mathcal{O}(d^2)$ mediante la norma de Frobenius con 0 alocaciones.
   - Tres estrategias de ajuste: `:dynamic` (recomendado/paper), `:fixed` (modelo clásico) y `:numerical` (optimización conjunta L-BFGS-B).
   - Métricas de bondad de ajuste estructuradas ($R^2$, MAE, RMSE, MaxAE, IAE).
   - Soporte polimórfico para `DataFrame`, `Matrix` y `Vector`.

3. **Visualización Integrada (`puig_plot.jl`)**:
   - Generación de paneles comparativos 2x2 y gráficos individuales desacoplados del núcleo numérico.

4. **Caso Límite $\lambda \to 0$**:
   - Tratamiento analítico continuo de la distribución límite de Maxwell-Boltzmann / $\chi$ generalizada.

---

## 📚 Estructura de la Documentación

- [Fundamentación Teórica](theory.md): Formulación matemática rigurosa, propiedades analíticas, cálculo de momentos y dimensión efectiva.
- [Tutorial y Ejemplos](tutorial.md): Guía paso a paso desde la definición de parámetros hasta el ajuste multivariante con DataFrames.
- [Referencia de la API](api.md): Documentación detallada de todas las funciones y tipos exportados.

