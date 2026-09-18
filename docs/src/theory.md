# Fundamentación Teórica

La **Distribución de Puig** generaliza la distribución $\chi$ no central continua al relajar la restricción de grados de libertad enteros, permitiendo un parámetro de dimensionalidad real $k \in \mathbb{R}^+$ ($k \ge 1$).

---

## 1. Definición y Función de Densidad (PDF)

Sea una variable aleatoria continua no negativa $X \in [0, \infty)$. Su función de densidad de probabilidad está parametrizada por tres magnitudes:
- $\lambda \ge 0$: parámetro de no centralidad (norma de las medias de las componentes).
- $k \ge 1$: dimensión efectiva continua.
- $T > 0$: parámetro de escala o precisión ($T = 1/\sigma^2$).

La PDF viene dada por:

```math
f_P(x; \lambda, k, T) = T \frac{x^{k/2}}{\lambda^{k/2 - 1}} \exp\!\left(-\frac{T}{2}(x^2 + \lambda^2)\right) I_{k/2 - 1}(x \lambda T), \quad x \ge 0
```

donde $I_\nu(z)$ denota la función de Bessel modificada de primera especie de orden $\nu = k/2 - 1$.

---

## 2. Identidad Fundamental con la Distribución $\chi^2$ No Central

Si consideramos un vector gaussiano $d$-dimensional $\mathbf{Y} \sim \mathcal{N}(\boldsymbol{\mu}, \sigma^2 \mathbf{I})$ y su norma euclídea $X = \|\mathbf{Y}\|$, se cumple que:

```math
\frac{X^2}{\sigma^2} \sim \text{NoncentralChisq}\left(k, \frac{\|\boldsymbol{\mu}\|^2}{\sigma^2}\right)
```

Haciendo la correspondencia $T = 1/\sigma^2$ y $\lambda = \|\boldsymbol{\mu}\|$, la variable aleatoria de Puig es exactamente:

```math
X = \sqrt{\frac{W}{T}}, \quad \text{donde } W \sim \text{NoncentralChisq}(k, \lambda^2 T)
```

Esta equivalencia permite realizar muestreo directo exacto y fundamenta el cálculo numéricamente estable de la función de distribución acumulada y de supervivencia.

---

## 3. Función de Supervivencia y Función Marcum-$Q$

La función de supervivencia $S(x) = P(X \ge x)$ viene dada por la función Marcum-$Q$ generalizada de orden $m = k/2$:

```math
S(x) = Q_{k/2}\left(\lambda \sqrt{T}, x \sqrt{T}\right)
```

En `DistributionsPuig.jl`, la evaluación de $Q_m(a, b)$ se calcula de forma nativa sin llamadas a librerías externas mediante el complemento de la acumulada de la $\chi^2$ no central:

```math
S(x) = \operatorname{ccdf}\left(\text{NoncentralChisq}(k, \lambda^2 T), x^2 T\right)
```

La función de distribución acumulada (CDF) es complementaria:

```math
F(x) = 1 - S(x) = \operatorname{cdf}\left(\text{NoncentralChisq}(k, \lambda^2 T), x^2 T\right)
```

---

## 4. Dimensión Efectiva de Correlación

Cuando se combinan $d$ variables correlacionadas en un indicador compuesto $X = \|\mathbf{Y}\|$, el espacio latente posee un número efectivo de grados de libertad menor que la dimensión geométrica $d$.

Inspirado en la razón de participación empleada en neurociencia y física estadística, la **dimensión efectiva** se define sobre los autovalores $\lambda_1, \dots, \lambda_d$ de la matriz de correlación $R$:

```math
k = \frac{\left(\sum_{i=1}^d \lambda_i\right)^2}{\sum_{i=1}^d \lambda_i^2}
```

Dado que la traza de la matriz de correlación es constante ($\operatorname{tr}(R) = \sum_{i=1}^d \lambda_i = d$) y que la suma de los cuadrados de los autovalores coincide idénticamente con el cuadrado de la norma de Frobenius $\|R\|_F^2$:

```math
\sum_{i=1}^d \lambda_i^2 = \|R\|_F^2 = d + 2 \sum_{i < j} R_{ij}^2
```

Se obtiene la expresión analítica simplificada:

```math
k = \frac{d^2}{d + 2 \sum_{i < j} R_{ij}^2}
```

Esta identidad permite evaluar la dimensión efectiva en tiempo $\mathcal{O}(d^2)$ con **cero alocaciones de memoria**, evitando el cómputo espectral mediante LAPACK ($\mathcal{O}(d^3)$).

---

## 5. Momentos Teóricos y Funciones de Laguerre Generalizadas

Los momentos teóricos no centrados $\mu_n = E[X^n]$ se obtienen analíticamente:

### Momentos Pares (Forma Cerrada Polinómica)
- **Segundo momento**:
  ```math
  \mu_2 = E[X^2] = \frac{k}{T} + \lambda^2
  ```
- **Cuarto momento**:
  ```math
  \mu_4 = E[X^4] = \mu_2^2 + \frac{2k}{T^2} + \frac{4\lambda^2}{T}
  ```

### Momentos Impares (Funciones de Laguerre Reales)
Para órdenes no enteros, los polinomios de Laguerre generalizados se extienden a funciones de Laguerre $L_n^{(\alpha)}(z)$ evaluadas vía la función hipergeométrica confluente de Kummer $M(a, b, z)$:

```math
L_n^{(\alpha)}(x) = \frac{\Gamma(n + \alpha + 1)}{\Gamma(n + 1)\Gamma(\alpha + 1)} \exp(x) M(\alpha + 1 + n, \alpha + 1, -x)
```

- **Media (Primer Momento $\mu_1$)**:
  ```math
  \mu_1 = E[X] = \sqrt{\frac{\pi}{2T}} L_{1/2}^{(k/2 - 1)}\left(-\frac{\lambda^2 T}{2}\right)
  ```
- **Tercer Momento ($\mu_3$)**:
  ```math
  \mu_3 = E[X^3] = \frac{3}{T} \sqrt{\frac{\pi}{2T}} L_{3/2}^{(k/2 - 1)}\left(-\frac{\lambda^2 T}{2}\right)
  ```

A partir de los momentos no centrados $\mu_1, \dots, \mu_4$, se calculan la varianza $\sigma^2 = \mu_2 - \mu_1^2$, la desviación estándar $\sigma$, el coeficiente de asimetría (*skewness*) $\gamma_1$ y la *kurtosis* $\gamma_2$.

---

## 6. Caso Límite $\lambda \to 0$ (Maxwell-Boltzmann / $\chi$ Generalizada)

Cuando $\lambda = 0$, la distribución de Puig colapsa continuamente a la distribución de Maxwell-Boltzmann generalizada con parámetro $B = 1/T$:

```math
f(x; k, B) = \frac{2^{1 - k/2} B^{-k/2}}{\Gamma(k/2)} x^{k-1} \exp\!\left(-\frac{x^2}{2B}\right), \quad x \ge 0
```

La supervivencia en este caso se expresa mediante la función gamma incompleta superior regularizada:

```math
S(x) = \frac{\Gamma(k/2, x^2 / (2B))}{\Gamma(k/2)}
```

---

## 7. Expansión Asintótica de la Función Marcum-$Q$ para Colas Superiores

### 7.1. Diagnóstico del Problema

Para argumentos moderados o grandes ($a \ge 10$ o $b \ge 15$), la evaluación de $Q_m(a, b)$ mediante `pnchisq` (la rutina subyacente de `Distributions.jl`) puede emitir advertencias `"full precision may not have been achieved in 'pnchisq'"` y producir ruido numérico espurio del orden de $10^{-14}$ en la cola superior lejana cuando la probabilidad real es $\sim 10^{-50}$.

### 7.2. Expansión Asintótica de la Cola Superior ($b > a$, $b - a \ge 4$, $ab \ge 30$)

Aprovechando la expansión asintótica de la función de Bessel modificada $I_{\nu}(z)$ (DLMF 10.40.1), se obtiene una fórmula directa con cancelaciones catastróficas nulas:

```math
Q_m(a, b) \sim \left(\frac{b}{a}\right)^{m - 1/2} \frac{\exp\!\left(-\frac{(b-a)^2}{2}\right)}{\sqrt{2\pi}(b-a)} \left[1 - \frac{b + a + \frac{4(m-1/2)^2 - 1}{4(b-a)}}{2ab(b-a)} + \mathcal{O}\!\left(\frac{1}{(ab(b-a))^2}\right)\right]
```

**Propiedades:**
- Factoriza directamente $\exp(-(b-a)^2/2)$, idéntico al diseño de `Puig_pdf`.
- Para $(b-a)^2/2 > 740$, devuelve limpiamente $0.0$ sin ruido de precisión.
- Complejidad computacional: $\sim 20\text{ ns}$ frente a $> 5\,\mu\text{s}$ de `pnchisq`.

### 7.3. Conmutación Automática

La función `marcumq(a, b, m)` conmuta automáticamente al régimen asintótico cuando se cumplen las condiciones $b - a \ge 4.0$ y $ab \ge 30.0$. Para la zona de transición ($|b - a| < 4$ o $ab < 30$), se conserva la evaluación mediante `NoncentralChisq`, donde `pnchisq` es estable porque las probabilidades son del orden $0.05 \dots 0.95$.

---

## 8. Momentos de Orden $n \in \mathbb{N}$ — Recurrencia de Tres Términos

### 8.1. Fórmula Analítica General

Dado $X \sim \text{Puig}(\lambda, k, T)$, el $n$-ésimo momento raw $\mu_n = \mathbb{E}[X^n]$ es:

```math
\mu_n = \left(\frac{2}{T}\right)^{n/2} \Gamma\!\left(\frac{n}{2} + 1\right) L_{n/2}^{(k/2 - 1)}\!\left(-\frac{\lambda^2 T}{2}\right)
```

donde $L_p^{(\alpha)}(z)$ es la función generalizada de Laguerre.

### 8.2. Relación de Recurrencia de Tres Términos ($\mathcal{O}(n)$ exacta)

A partir de la relación de contigüidad de las funciones de Laguerre:

```math
\mu_n = \left(\frac{2n - 4 + k}{T} + \lambda^2\right) \mu_{n-2} - \frac{(n-2)(n + k - 4)}{T^2} \mu_{n-4}, \quad n \ge 4
```

**Semillas:**
- $\mu_0 = 1.0$
- $\mu_1 = \sqrt{\frac{\pi}{2T}} \, L_{1/2}^{(k/2-1)}\!\left(-\frac{\lambda^2 T}{2}\right)$
- $\mu_2 = \frac{k}{T} + \lambda^2$
- $\mu_3 = \frac{3}{T}\sqrt{\frac{\pi}{2T}} \, L_{3/2}^{(k/2-1)}\!\left(-\frac{\lambda^2 T}{2}\right)$

**Propiedades:**
- Solo se necesitan **2 llamadas** a `laguerre_real` ($\mu_1$ y $\mu_3$); el resto se obtiene con operaciones aritméticas elementales.
- Error relativo $< 2 \times 10^{-16}$ (precisión máquina estricta).
- Para $\lambda = 0$ (Maxwell-Boltzmann), la fórmula cerrada independiente es: $\mu_n = (2B)^{n/2} \, \Gamma\!\left(\frac{k+n}{2}\right) / \Gamma\!\left(\frac{k}{2}\right)$.

