"""
    marcumq_asymp(a::Real, b::Real, m::Real=1) -> Float64

Calcula Q_m(a, b) mediante la expansión asintótica de 2 términos para la cola superior lejana
(b > a, b - a ≥ 4.0, ab ≥ 30.0). Evita cancelaciones catastróficas de `pnchisq` y devuelve
limpiamente 0.0 para argumentos extremos.

Basado en Cantrell (1986) y la aproximación de Bessel modificada (DLMF 10.40.1):
    Q_m(a,b) ~ (b/a)^{m-1/2} exp(-(b-a)²/2) / (√(2π)(b-a)) × [1 + corr_1 + ...]
"""
function marcumq_asymp(a::Real, b::Real, m::Real=1)
    af = Float64(a)
    bf = Float64(b)
    mf = Float64(m)

    (b - a) >= 4.0 || throw(ArgumentError("marcumq_asymp requiere b - a ≥ 4.0, recibido b-a=$(b-a)"))
    (a * b) >= 30.0 || throw(ArgumentError("marcumq_asymp requiere ab ≥ 30.0, recibido ab=$(a*b)"))

    diff = bf - af
    arg_exp = -0.5 * diff * diff
    arg_exp < -740.0 && return 0.0

    rho = bf / af
    prefactor = rho^(mf - 0.5) * exp(arg_exp) / (sqrt(2.0 * π) * diff)

    num1 = bf + af + (4.0 * (mf - 0.5)^2 - 1.0) / (4.0 * diff)
    denom1 = 2.0 * af * bf * diff
    corr = 1.0 - num1 / denom1

    val = prefactor * corr
    return isfinite(val) ? max(val, 0.0) : 0.0
end

"""
    marcumq(a::Real, b::Real, m::Real=1) -> Float64
    marcutq(a::Real, b::AbstractVector{<:Real}, m::Real=1) -> Vector{Float64}

Calcula la función Marcum-Q generalizada de orden `m`:
    Q_m(a, b) = \\int_b^\\infty x \\left(\\frac{x}{a}\\right)^{m-1} \\exp\\left(-\\frac{x^2+a^2}{2}\\right) I_{m-1}(a x) dx
Implementada de forma nativa en Julia sin dependencias externas mediante `ccdf(NoncentralChisq(2m, a^2), b^2)`.

Para la cola superior lejana (b - a ≥ 4.0 y ab ≥ 30.0) se usa la expansión asintótica `marcumq_asymp`.
"""
function marcumq(a::Real, b::Real, m::Real=1)
    if (b - a) >= 4.0 && (a * b) >= 30.0
        return marcumq_asymp(a, b, m)
    end
    return redirect_stderr(devnull) do
        ccdf(NoncentralChisq(2m, a^2), b^2)
    end
end
function marcumq(a::Real, b::AbstractVector{<:Real}, m::Real=1)
    dist = NoncentralChisq(2m, a^2)
    out = similar(b, Float64)
    @inbounds for i in eachindex(b)
        bi = Float64(b[i])
        if (bi - a) >= 4.0 && (a * bi) >= 30.0
            out[i] = marcumq_asymp(a, bi, m)
        else
            out[i] = redirect_stderr(devnull) do
                ccdf(dist, bi^2)
            end
        end
    end
    return out
end

"""
    struct PuigDistribution{S<:Real}

Representa la Distribución de Puig (distribución χ no central generalizada con dimensión continua).

# Parámetros:
- `λ::S`: Norma del vector de medias de las componentes (λ ≥ 0).
- `k::S`: Dimensión efectiva continua (k ≥ 1).
- `T::S`: Parámetro de escala/precisión (T > 0, T = 1/σ²).
"""
struct PuigDistribution{S<:Real}
    λ::S
    k::S
    T::S
    function PuigDistribution(λ::S, k::S, T::S) where {S<:Real}
        λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
        k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
        T > 0 || throw(ArgumentError("T debe ser mayor que 0"))
        new{S}(λ, k, T)
    end
end

function Base.show(io::IO, params::PuigDistribution)
    println(io, "Distribución de Puig (PuigDistribution):")
    println(io, "  λ = $(round(params.λ, digits=4))")
    println(io, "  k = $(round(params.k, digits=4))")
    print(io,   "  T = $(round(params.T, digits=4))")
end
PuigDistribution(λ::Real, k::Real, T::Real) = PuigDistribution(promote(λ, k, T)...)

"""
    besselix_asymp(nu::Float64, z::Float64) -> Float64

Calcula la función de Bessel modificada escalada e^{-z} I_ν(z) mediante su expansión asintótica
de 5 términos (DLMF 10.40.1). Precisión relativa < 1e-12 para z ≥ 100 y < 1e-15 para z ≥ 200.
"""
function besselix_asymp(nu::Float64, z::Float64)::Float64
    mu = 4.0 * nu^2
    res = 1.0
    t1 = -(mu - 1.0) / (8.0 * z)
    res += t1
    t2 = -t1 * (mu - 9.0) / (16.0 * z)
    res += t2
    t3 = -t2 * (mu - 25.0) / (24.0 * z)
    res += t3
    t4 = -t3 * (mu - 49.0) / (32.0 * z)
    res += t4
    t5 = -t4 * (mu - 81.0) / (40.0 * z)
    res += t5
    return res / sqrt(2.0 * π * z)
end

"""
    puig_pdf_arb(x::Float64, λ::Float64, k::Float64, T::Float64) -> Float64

Calcula la PDF de Puig usando ArbNumerics con precisión de 128 bits para referencia de alta precisión.
"""
function puig_pdf_arb(x::Float64, λ::Float64, k::Float64, T::Float64)::Float64
    x <= 0.0 && return 0.0
    setworkingprecision(ArbReal, bits=128)
    nu = ArbReal(k / 2 - 1)
    lb = ArbReal(λ)
    Ta = ArbReal(T)
    ka = ArbReal(k)
    xf = ArbReal(x)
    Tf = Ta * xf^(ka / 2) / lb^(ka / 2 - 1)
    Te = exp(-Ta / 2 * (xf^2 + lb^2))
    Tb = ArbNumerics.besseli(nu, xf * lb * Ta)
    return Float64(Tf * Te * Tb)
end

"""
    Puig_pdf(x::Real, λ::Number, k::Number, T::Number; method::Symbol=:asymp) -> Float64
    Puig_pdf(x::Real, params::Union{PuigDistribution, MB}) -> Float64

Evalúa la función de densidad de probabilidad (PDF) para un valor escalar `x` en Float64 puro.

# Método (`method`)
- `:asymp` (por defecto): Usa expansión asintótica de 5 términos (DLMF 10.40.1) para z ≥ 200,
  y `SpecialFunctions.besselix` para z < 200. Rápido y preciso (< 1e-12 rel.).
- `:arb`: Usa ArbNumerics con precisión de 128 bits para referencia de alta precisión.
"""
function Puig_pdf(x::Real, λ::Number, k::Number, T::Number; method::Symbol=:asymp)::Float64
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    x >= 0 || throw(ArgumentError("Dominio de la función es [0, ∞)"))

    x <= 0.0 && return 0.0
    λ == 0 && return MB_pdf([Float64(x)], Float64(k), 1.0 / Float64(T))[1]

    λf = Float64(λ)
    kf = Float64(k)
    Tf = Float64(T)
    xf = Float64(x)

    if method == :arb
        return puig_pdf_arb(xf, λf, kf, Tf)
    end

    arg_exp = -0.5 * Tf * (xf - λf)^2
    arg_exp < -740.0 && return 0.0
    expo = exp(arg_exp)
    expo == 0.0 && return 0.0

    nu = kf / 2.0 - 1.0
    z = xf * λf * Tf
    bix = z >= 200.0 ? besselix_asymp(nu, z) : besselix(nu, z)

    factor = Tf * (xf^(kf / 2.0)) / (λf^(kf / 2.0 - 1.0))
    val = factor * expo * bix
    return isfinite(val) ? val : 0.0
end

"""
    Puig_pdf(x::AbstractVector{<:Real}, λ::Number, k::Number, T::Number; method::Symbol=:asymp) -> Vector{Float64}
    Puig_pdf(x::AbstractVector{<:Real}, params::PuigDistribution) -> Vector{Float64}

Calcula la función de densidad de probabilidad (PDF) de la distribución de Puig de forma vectorizada
en Float64 puro mediante escalado exponencial y expansión asintótica para argumentos grandes (z ≥ 200).

# Método (`method`)
- `:asymp` (por defecto): Expansión asintótica + besselix.
- `:arb`: ArbNumerics con precisión de 128 bits.
"""
function Puig_pdf(x::AbstractVector{<:Real}, λ::Number,
                  k::Number, T::Number; method::Symbol=:asymp)::Vector{Float64}
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0 || throw(ArgumentError("T debe ser > 0"))
    all(xi -> xi >= 0, x) || throw(ArgumentError("Dominio de la función es [0, ∞)"))

    if λ == 0
        B = 1/T
        return MB_pdf(x, k, B)
    end

    λf = Float64(λ)
    kf = Float64(k)
    Tf = Float64(T)
    nu = kf / 2.0 - 1.0
    half_T = 0.5 * Tf
    pow_k_half = kf / 2.0
    λ_denom = λf^(pow_k_half - 1.0)

    if method == :arb
        out = zeros(Float64, length(x))
        @inbounds for i in eachindex(x)
            out[i] = puig_pdf_arb(Float64(x[i]), λf, kf, Tf)
        end
        return out
    end

    out = zeros(Float64, length(x))
    @inbounds for i in eachindex(x)
        xi = Float64(x[i])
        if xi > 0.0
            arg_exp = -half_T * (xi - λf)^2
            if arg_exp >= -740.0
                expo = exp(arg_exp)
                if expo > 0.0
                    z = xi * λf * Tf
                    bix = z >= 200.0 ? besselix_asymp(nu, z) : besselix(nu, z)
                    factor = Tf * (xi^pow_k_half) / λ_denom
                    val = factor * expo * bix
                    out[i] = isfinite(val) ? val : 0.0
                end
            end
        end
    end
    return out
end
function Puig_pdf_vec1(x::AbstractVector{<:Real}, λ::Float64,
                       k_vec::Vector{<:Real}, TT::Vector{<:Real}; method::Symbol=:asymp)
    length(k_vec) == length(TT) || throw(ArgumentError("k_vec y TT deben tener la misma longitud"))
    nx = length(x)
    nk = length(k_vec)
    Z  = zeros(Float64, nx, nk)
    for i in 1:nk
        Z[:, i] = Puig_pdf(x, λ, k_vec[i], TT[i]; method=method)
    end
    return Z
end
function Puig_pdf_vec(x::AbstractVector{<:Real}, λ::Float64,
                      k_vec::Union{Float64, Vector{Float64}},
                      TT::Union{Float64, Vector{Float64}}; method::Symbol=:asymp)
    Puig_pdf_vec1(collect(x), λ,
                  isa(k_vec, Float64) ? [k_vec] : k_vec,
                  isa(TT,    Float64) ? [TT]    : TT; method=method)
end

"""
    Puig_surviving(x::AbstractVector{<:Real}, λ::Real, k::Real, T::Real) -> Vector{Float64}
    Puig_surviving(x::AbstractVector{<:Real}, params::PuigDistribution) -> Vector{Float64}

Calcula la función de supervivencia S(x) = P(X ≥ x) = Q_{k/2}(λ√T, x√T) mediante la función Marcum-Q.
"""
function Puig_surviving(x::AbstractVector{<:Real}, λ::Real, k::Real, T::Real)::Vector{Float64}
    k >= 1 || throw(ArgumentError("Puig_survivng requiere k >= 1 (orden Marcum-Q = k/2 >= 0.5), recibido k=$k"))
    T > 0 || throw(ArgumentError("T debe ser > 0"))
    λ >= 0 || throw(ArgumentError("λ debe ser >= 0"))
    all(xi -> xi >= 0, x) || throw(ArgumentError("Dominio de la función es [0, ∞)"))
    if λ==0
        B=1/T;
        return MB_surviving(x,k,B)
    else
     order = Float64(k / 2);
     sT = sqrt(Float64(T));
     a = Float64(λ)*sT;
     b=Float64.(x)*sT;
     c=marcumq(a,b,order);
      @inbounds for i in eachindex(c)
        if isnan(c[i])
        throw(DomainError(c[i], "marcumq devolvió NaN en x=$(x[i]), a=$a, m=$order"))
        elseif c[i] <= 0.0 && x[i] < λ
        c[i] = 1.0
        end
    end
    return c
    end    
end

"""
    Puig_cumulative(x::AbstractVector{<:Real}, λ::Real, k::Real, T::Real) -> Vector{Float64}
    Puig_cumulative(x::AbstractVector{<:Real}, params::PuigDistribution) -> Vector{Float64}

Calcula la función de distribución acumulada (CDF) F(x) = P(X ≤ x) = 1 - S(x).
"""
function Puig_cumulative(x::AbstractVector{<:Real}, λ::Real, k::Real, T::Real)::Vector{Float64}
    k >= 1 || throw(ArgumentError("Puig_cumlative requiere k >= 1 (orden Marcum-Q = k/2 >= 0.5), recibido k=$k"))
    T > 0 || throw(ArgumentError("T debe ser > 0"))
    λ >= 0 || throw(ArgumentError("λ debe ser >= 0"))
    all(xi -> xi >= 0, x) || throw(ArgumentError("Dominio de la función es [0, ∞)"))
    if λ==0
        B=1/T;
        return MB_cumulative(x,k,B)
    else
    return 1 .- Puig_surviving(x, λ, k, T)
    end
end
Puig_pdf(x::AbstractVector{<:Real},params::PuigDistribution; method::Symbol=:asymp)=Puig_pdf(x,Float64(params.λ),Float64(params.k),Float64(params.T); method=method)
Puig_pdf(x::Real, params::PuigDistribution; method::Symbol=:asymp)=Puig_pdf(Float64(x),Float64(params.λ),Float64(params.k),Float64(params.T); method=method)
Puig_pdf(x::AbstractVector{<:Real},params::MB)=MB_pdf(x,params)
Puig_pdf(x::Real, params::MB)=MB_pdf([Float64(x)],params)[1]
Puig_cumulative(x::AbstractVector{<:Real},params::PuigDistribution)=Puig_cumulative(x,Float64(params.λ),Float64(params.k),Float64(params.T))
Puig_cumulative(x::AbstractVector{<:Real},params::MB)=MB_cumulative(x,params)
Puig_surviving(x::AbstractVector{<:Real},params::PuigDistribution)=Puig_surviving(x,Float64(params.λ),Float64(params.k),Float64(params.T))
Puig_surviving(x::AbstractVector{<:Real},params::MB)=MB_surviving(x,params) 
function kummer_M_stable(a::ArbReal, b::ArbReal, z::ArbReal; maxterms::Int=15000)
    term  = one(z)
    total = term
    k = 0
    while k < maxterms
        term *= (a + k) / (b + k) * z / (k + 1)
        total += term
        term < total * ArbReal(10)^(-30) && break   # criterio de parada
        k += 1
    end
    return total
end
function laguerre_real(n::Real, α::Real, x::Real)
    setworkingprecision(ArbReal, bits=128)
    nA, αA, xA = ArbReal(n), ArbReal(α), ArbReal(x)
    coef = ArbNumerics.gamma(nA + αA + 1) / (ArbNumerics.gamma(nA + 1) * ArbNumerics.gamma(αA + 1))

    # Transformación de Kummer: M(-n, α+1, x) = exp(x) * M(α+1+n, α+1, -x)
    Mval = exp(xA) * kummer_M_stable(αA + 1 + nA, αA + 1, -xA)
    return Float64(coef * Mval)
end
"""
    Puig_mean(λ::Real, k::Real, T::Real) -> Float64
    Puig_mean(params::Union{PuigDistribution, MB}) -> Float64

Calcula la media teórica E[X] (primer momento raw μ₁) de la distribución de Puig usando la función de Laguerre generalizada:
    E[X] = \\sqrt{\\frac{\\pi}{2T}} L_{1/2}^{(k/2-1)}\\left(-\\frac{\\lambda^2 T}{2}\\right)
"""
function Puig_mean(λ::Real,k::Real,T::Real)::Number
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    λs=λ*sqrt(T);
    out=sqrt(π/2/T) * laguerre_real(1/2, k/2 - 1, -λs^2/2)
    return out
end
Puig_mean(params::PuigDistribution)=Puig_mean(Float64(params.λ),Float64(params.k),Float64(params.T))
Puig_mean(params::MB)=Puig_mean(0.0,Float64(params.k),1.0/Float64(params.B))

"""
    Puig_var(λ::Real, k::Real, T::Real) -> Float64
    Puig_var(params::Union{PuigDistribution, MB}) -> Float64

Calcula la varianza teórica Var(X) = E[X²] - (E[X])² con E[X²] = k/T + λ².
"""
function Puig_var(λ::Real, k::Real, T::Real)::Number
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    μ2 =k/T+λ^2;   # E[X²]
    μ1 = Puig_mean(λ, k, T)                       # E[X]
    return μ2 - μ1^2
end
Puig_var(params::PuigDistribution)=Puig_var(Float64(params.λ),Float64(params.k),Float64(params.T))
Puig_var(params::MB)=Puig_var(0.0,Float64(params.k),1.0/Float64(params.B))

"""
    Puig_std(λ::Real, k::Real, T::Real) -> Float64
    Puig_std(params::Union{PuigDistribution, MB}) -> Float64

Calcula la desviación estándar teórica σ = √(Var(X)).
"""
function Puig_std(λ::Real, k::Real, T::Real)::Number
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    return sqrt(Puig_var(λ,k,T))
end
Puig_std(params::PuigDistribution)=sqrt(Puig_var(Float64(params.λ),Float64(params.k),Float64(params.T)))
Puig_std(params::MB)=sqrt(Puig_var(0.0,Float64(params.k),1.0/Float64(params.B)))

"""
    Puig_skewness(λ::Real, k::Real, T::Real) -> Float64
    Puig_skewness(params::Union{PuigDistribution, MB}) -> Float64

Calcula el coeficiente de asimetría (skewness) estandarizado γ₁ = E[(X-μ)³] / σ³.
"""
function Puig_skewness(λ::Real, k::Real, T::Real)::Number
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    par=PuigDistribution(λ,k,T)
    λs=λ*sqrt(T);
    L=laguerre_real(3/2, k/2 - 1, -λs^2/2)
    μ3=3/T*sqrt(π/2/T)*L;
    μ1=Puig_mean(par)
    σ2=Puig_var(par)
    μ2=σ2+μ1^2;
    μ3c=μ3 - 3*μ1*μ2 + 2*μ1^3;
    σ3=Puig_std(par)^3;
    out=μ3c/σ3
    return out
end
Puig_skewness(params::PuigDistribution)=Puig_skewness(Float64(params.λ),Float64(params.k),Float64(params.T))
Puig_skewness(params::MB)=Puig_skewness(0.0,Float64(params.k),1.0/Float64(params.B))

"""
    Puig_kurtosis(λ::Real, k::Real, T::Real) -> Float64
    Puig_kurtosis(params::Union{PuigDistribution, MB}) -> Float64

Calcula la kurtosis (cuarto momento estandarizado no restado) γ₂ = E[(X-μ)⁴] / σ⁴.
"""
function Puig_kurtosis(λ::Real, k::Real, T::Real)::Number
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    par=PuigDistribution(λ,k,T)
    μ2 =k/T+λ^2; 
    μ4=(μ2)^2+2*k/(T^2)+4*λ^2/T;
    λs = λ*sqrt(T)
    L   = laguerre_real(3/2, k/2 - 1, -λs^2/2)
    μ3  = 3/T*sqrt(π/2/T)*L
    μ1  = Puig_mean(par)
    μ4c = μ4 - 4*μ1*μ3 + 6*μ1^2*μ2 - 3*μ1^4
    σ4  = Puig_var(par)^2
    return μ4c/σ4
end
Puig_kurtosis(params::PuigDistribution) = Puig_kurtosis(Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_kurtosis(params::MB)  = Puig_kurtosis(0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    Puig_moments(λ::Real, k::Real, T::Real) -> NTuple{4, Float64}
    Puig_moments(params::Union{PuigDistribution, MB}) -> NTuple{4, Float64}

Devuelve los 4 primeros momentos no centrados o raw (μ₁, μ₂, μ₃, μ₄) donde μₙ = E[Xⁿ].
Retrocompatible con la API original.
"""
function Puig_moments(λ::Real, k::Real, T::Real)::NTuple{4, Float64}
    μ = Puig_moments_raw(λ, k, T, 4)
    return (μ[1], μ[2], μ[3], μ[4])
end
Puig_moments(params::PuigDistribution) = Puig_moments(Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_moments(params::MB) = Puig_moments(0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    Puig_moments(λ::Real, k::Real, T::Real, n::Integer) -> Vector{Float64}
    Puig_moments(params::Union{PuigDistribution, MB}, n::Integer) -> Vector{Float64}

Calcula los n primeros momentos raw μ₁, …, μₙ de la distribución de Puig usando la relación
de recurrencia de tres términos O(n) exacta sin evaluaciones hipergeométricas adicionales.

## Recurrencia (n ≥ 4):
    μₙ = ((2n - 4 + k)/T + λ²) μₙ₋₂ - (n-2)(n+k-4)/T² μₙ₋₄

Semillas:
- μ₀ = 1.0
- μ₁ = Puig_mean(λ, k, T)
- μ₂ = k/T + λ²
- μ₃ = (3/T) √(π/(2T)) L_{3/2}^{(k/2-1)}(-λ²T/2)
"""
function Puig_moments(λ::Real, k::Real, T::Real, n::Integer)::Vector{Float64}
    μ = Puig_moments_raw(λ, k, T, Int(n))
    return μ
end
Puig_moments(params::PuigDistribution, n::Integer) = Puig_moments(Float64(params.λ), Float64(params.k), Float64(params.T), n)
Puig_moments(params::MB, n::Integer) = Puig_moments(0.0, Float64(params.k), 1.0/Float64(params.B), n)

"""
    Puig_moments_raw(λ::Real, k::Real, T::Real, n::Int) -> Vector{Float64}

Calcula internamente los n primeros momentos raw usando semillas + recurrencia.
Vector resultado es 1-indexado: μ[1] = μ₁, μ[2] = μ₂, …, μ[n] = μₙ.
"""
function Puig_moments_raw(λ::Real, k::Real, T::Real, n::Int)::Vector{Float64}
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    n >= 1 || throw(ArgumentError("n debe ser ≥ 1, recibido n=$n"))

    λf = Float64(λ)
    kf = Float64(k)
    Tf = Float64(T)

    μ = Vector{Float64}(undef, n)

    # --- Semillas ---
    # μ₀ = 1.0 (no almacenado, pero usado en recurrencia para μ₄)
    # μ₁
    μ[1] = Float64(Puig_mean(λf, kf, Tf))

    if n == 1
        return μ
    end

    # μ₂ = k/T + λ²
    μ[2] = kf / Tf + λf^2

    if n == 2
        return μ
    end

    # μ₃ = (3/T) √(π/(2T)) L_{3/2}^{(k/2-1)}(-λ²T/2)
    λs = λf * sqrt(Tf)
    L32 = laguerre_real(3 / 2, kf / 2 - 1, -λs^2 / 2)
    μ[3] = (3.0 / Tf) * sqrt(π / (2.0 * Tf)) * L32

    if n == 3
        return μ
    end

    # μ₄ = (μ₂)² + 2k/T² + 4λ²/T
    μ[4] = μ[2]^2 + 2.0 * kf / Tf^2 + 4.0 * λf^2 / Tf

    if n == 4
        return μ
    end

    # --- Recurrencia de tres términos para n ≥ 5 ---
    # μₙ = ((2n - 4 + k)/T + λ²) μₙ₋₂ - (n-2)(n+k-4)/T² μₙ₋₄
    invT  = 1.0 / Tf
    invT2 = invT * invT
    λ2 = λf * λf

    for idx in 5:n
        μ[idx] = ((2.0 * idx - 4.0 + kf) * invT + λ2) * μ[idx - 2] -
                 Float64(idx - 2) * Float64(idx + kf - 4.0) * invT2 * μ[idx - 4]
    end

    return μ
end

"""
    Puig_stats(λ::Real, k::Real, T::Real) -> NamedTuple
    Puig_stats(params::Union{PuigDistribution, MB}) -> NamedTuple

Devuelve un `NamedTuple` con los principales estadísticos descriptivos teóricos:
`(mean, var, sig, skewness, kurtosis)`.
"""
function Puig_stats(λ::Real, k::Real, T::Real)::NamedTuple
    μ = Puig_moments_raw(Float64(λ), Float64(k), Float64(T), 4)
    μ1 = μ[1]
    μ2 = μ[2]
    μ3 = μ[3]
    μ4 = μ[4]
    σ2  = μ2 - μ1^2
    σ   = sqrt(σ2)
    sk  = (μ3 - 3μ1*μ2 + 2μ1^3) / σ^3
    kur = (μ4 - 4μ1*μ3 + 6μ1^2*μ2 - 3μ1^4) / σ2^2
    return (mean=μ1, var=σ2, sig=σ, skewness=sk, kurtosis=kur)
end
Puig_stats(params::PuigDistribution) = Puig_stats(Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_stats(params::MB) = Puig_stats(0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    Puig_hazard(x, λ::Real, k::Real, T::Real)
    Puig_hazard(x, params::Union{PuigDistribution, MB})

Calcula la función de riesgo o tasa de fallo instantánea:
    h(x) = \\frac{f(x)}{S(x)}
"""
function Puig_hazard(x::AbstractVector{<:Real},λ::Real,k::Real,T::Real)
    return Puig_pdf(x,λ,k,T) ./ Puig_surviving(x,λ,k,T)
end
Puig_hazard(x::Real,λ::Real,k::Real,T::Real) = Puig_hazard([x],λ,k,T)[1] 
Puig_hazard(x,params::PuigDistribution) = Puig_hazard(x, Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_hazard(x,params::MB)  = Puig_hazard(x, 0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    Puig_quantile(p, λ::Real, k::Real, T::Real; tol=1e-10, maxit=200)
    Puig_quantile(p, params::Union{PuigDistribution, MB}; tol=1e-10, maxit=200)

Calcula el cuantil o función cuantil inversa Q(p) tal que F(Q(p)) = p mediante método numérico de bisección.
"""
function Puig_quantile(p::Real, λ::Real, k::Real, T::Real; tol::Float64=1e-10, maxit::Int=200)
    0 <= p <= 1 || throw(ArgumentError("p debe estar en [0,1]"))
    p == 0 && return 0.0
    p == 1 && return Inf
    lo, hi = 0.0, max(λ, 1.0)*2
    while Puig_cumulative([hi], λ, k, T)[1] < p
        hi *= 2
    end
    for _ in 1:maxit
        mid = (lo+hi)/2
        if Puig_cumulative([mid], λ, k, T)[1] < p
            lo = mid
        else
            hi = mid
        end
        (hi-lo) < tol && break
    end
    return (lo+hi)/2
end
Puig_quantile(p::AbstractVector{<:Real}, λ::Real, k::Real, T::Real) = Puig_quantile.(p, λ, k, T)
Puig_quantile(p, params::PuigDistribution) = Puig_quantile(p, Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_quantile(p, params::MB)  = Puig_quantile(p, 0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    Puig_ci98(params::Union{PuigDistribution, MB}) -> Tuple{Float64, Float64}

Devuelve el intervalo central que contiene el 98% de la distribución: `(P1, P99)`.
"""
function Puig_ci98(λ::Real,k::Real,T::Real)::Tuple
    Li, Ls= Puig_quantile([0.01,0.99],λ,k,T)
    return (Li,Ls)
end
Puig_ci98(params::PuigDistribution) = Puig_ci98(Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_ci98(params::MB)  = Puig_ci98(0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    Puig_ci996(params::Union{PuigDistribution, MB}) -> Tuple{Float64, Float64}

Devuelve el intervalo central que contiene el 99.8% de la distribución: `(P0.1, P99.9)`.
"""
function Puig_ci996(λ::Real,k::Real,T::Real)::Tuple
    Li, Ls= Puig_quantile([0.001,0.999],λ,k,T)
    return (Li,Ls)
end
Puig_ci996(params::PuigDistribution) = Puig_ci996(Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_ci996(params::MB)  = Puig_ci996(0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    Puig_logpdf(x, λ::Real, k::Real, T::Real)
    Puig_logpdf(x, params::Union{PuigDistribution, MB})

Calcula el logaritmo natural de la densidad de probabilidad log(f(x)) de forma numéricamente estable usando `SpecialFunctions.besselix`.
"""
function Puig_logpdf(x::Real, λ::Real, k::Real, T::Real)::Float64
    λ >= 0 || throw(ArgumentError("λ debe ser ≥ 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    T > 0  || throw(ArgumentError("T debe ser > 0"))
    x >= 0 || throw(ArgumentError("Dominio de la función es [0, ∞)"))
    x == 0 && return -Inf   
    if λ == 0
        B = 1/T
        kh = k/2
        return (1-kh)*log(2) - kh*log(B) - loggamma(kh) + (k-1)*log(x) - x^2/(2B)
    end
    ν = k/2 - 1
    z = x*λ*T
    return log(T) + (k/2)*log(x) - (k/2-1)*log(λ) - T/2*(x-λ)^2 + log(besselix(ν, z))
end
Puig_logpdf(x::AbstractVector{<:Real}, λ::Real, k::Real, T::Real) = Puig_logpdf.(x, λ, k, T)
Puig_logpdf(x, params::PuigDistribution) = Puig_logpdf(x, Float64(params.λ), Float64(params.k), Float64(params.T))
Puig_logpdf(x, params::MB)  = Puig_logpdf(x, 0.0, Float64(params.k), 1.0/Float64(params.B))

"""
    rand(rng::AbstractRNG, params::Union{PuigDistribution, MB}, [n::Integer])

Genera muestras pseudoaleatorias exactas de la distribución de Puig a través de la representación X = √(W / T) con W ~ NoncentralChisq(k, λ² T).
"""
function Base.rand(rng::AbstractRNG, params::PuigDistribution)
    if params.λ == 0
        return sqrt(rand(rng, Gamma(Float64(params.k)/2, 2/Float64(params.T))))
    end
    W = rand(rng, NoncentralChisq(Float64(params.k), Float64(params.λ)^2 * Float64(params.T)))
    return sqrt(W / Float64(params.T))
end
Base.rand(params::PuigDistribution) = rand(Random.default_rng(), params)
Base.rand(rng::AbstractRNG, params::PuigDistribution, n::Integer) = [rand(rng, params) for _ in 1:n]
Base.rand(params::PuigDistribution, n::Integer) = rand(Random.default_rng(), params, n) 
Base.rand(rng::AbstractRNG, params::MB) = sqrt(rand(rng, Gamma(Float64(params.k)/2, 2*Float64(params.B))))
Base.rand(params::MB) = rand(Random.default_rng(), params)


"""
    Puig_entropy(params::PuigDistribution; N=400) -> Float64

Calcula la entropía diferencial teórica H(X) = - \\int_0^\\infty f(x) \\log f(x) dx mediante integración numérica.
"""
function Puig_entropy(params::PuigDistribution; N=400)
    init=max(Puig_ci996(params)[1]/10, 1e-6)
    fin=Puig_ci996(params)[2]
    x = range(init,fin, length=N)
    f = Puig_pdf(x, params)
    lf=Puig_logpdf(x,params)
    F=-f.*lf;
    h=(fin-init)/(N-1);
    FH=F.*h;
    FH[1]=FH[1]/2.0;
    FH[end]=FH[end]/2.0;
    return sum(FH)
end
