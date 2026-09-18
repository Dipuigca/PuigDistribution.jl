"""
    struct MB{S<:Real}

Representa la distribución límite de Maxwell-Boltzmann (o distribución χ generalizada con dimensión real `k`) cuando λ → 0.

# Parámetros:
- `k::S`: Dimensión / grados de libertad reales (k ≥ 1).
- `B::S`: Parámetro de escala o dispersión térmica (B > 0, donde B = 1/T).
"""
struct MB{S<:Real}
    k::S
    B::S
    function MB(k::S, B::S) where {S<:Real}
        k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
        B > 0 || throw(ArgumentError("B debe ser mayor que 0"))
        new{S}(k, B)
    end
end
MB(k::Real, B::Real) = MB(promote(k, B)...)

"""
    MB_pdf(x::AbstractVector{<:Real}, k::Real, B::Real) -> Vector{Float64}
    MB_pdf(x::AbstractVector{<:Real}, params::MB) -> Vector{Float64}

Calcula la función de densidad de probabilidad (PDF) para el caso límite λ → 0:
    f(x) = \\frac{2^{1-k/2} B^{-k/2}}{\\Gamma(k/2)} x^{k-1} \\exp\\left(-\\frac{x^2}{2B}\\right)
"""
function MB_pdf(x::AbstractVector{<:Real}, k::Real, B::Real)::Vector{Float64}
    B > 0 || throw(ArgumentError("B debe ser > 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    all(xi -> xi >= 0, x) || throw(ArgumentError("Dominio de la función es [0, ∞)"))
    k=Float64(k); B=Float64(B);
    kh=k/2;
    E=exp.(-x.^2/2/B);
    G=gamma(kh);
    C=2^(1-kh)*B^(-kh)/G;
    xt=x.^(k-1);
    return xt.*E.*C
end

function Base.show(io::IO, params::MB)
    println(io, "MB (λ=0):")
    println(io, "  k = $(round(params.k, digits=4))")
    print(io,   "  B = $(round(params.B, digits=4)) (T = $(round(1.0/params.B, digits=4)))")
end

"""
    MB_surviving(x::AbstractVector{<:Real}, k::Real, B::Real) -> Vector{Float64}
    MB_surviving(x::AbstractVector{<:Real}, params::MB) -> Vector{Float64}

Calcula la función de supervivencia S(x) = P(X ≥ x) para el caso límite λ → 0 utilizando la función gamma incompleta superior.
"""
function MB_surviving(x::AbstractVector{<:Real}, k::Real, B::Real)::Vector{Float64}
    B > 0 || throw(ArgumentError("B debe ser > 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    all(xi -> xi >= 0, x) || throw(ArgumentError("Dominio de la función es [0, ∞)"))
    k=Float64(k); B=Float64(B);
    kh=k/2;
    P=x.^2/2/B;
    G=1/gamma(kh);
    Up=gamma.(kh,P)
    return G.*Up
end

"""
    MB_cumulative(x::AbstractVector{<:Real}, k::Real, B::Real) -> Vector{Float64}
    MB_cumulative(x::AbstractVector{<:Real}, params::MB) -> Vector{Float64}

Calcula la función de distribución acumulada F(x) = 1 - S(x) para el caso límite λ → 0.
"""
function MB_cumulative(x::AbstractVector{<:Real}, k::Real, B::Real)::Vector{Float64}
    B > 0 || throw(ArgumentError("B debe ser > 0"))
    k >= 1 || throw(ArgumentError("k debe ser ≥ 1"))
    all(xi -> xi >= 0, x) || throw(ArgumentError("Dominio de la función es [0, ∞)"))    
    out=1.0.- MB_surviving(x,k,B)
    return out
end
MB_pdf(x::AbstractVector{<:Real},params::MB)=MB_pdf(x,Float64(params.k),Float64(params.B))
MB_cumulative(x::AbstractVector{<:Real},params::MB)=MB_cumulative(x,Float64(params.k),Float64(params.B))
MB_surviving(x::AbstractVector{<:Real},params::MB)=MB_surviving(x,Float64(params.k),Float64(params.B))
