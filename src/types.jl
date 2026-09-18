# ─── PuigDistribution como subtipo de ContinuousUnivariateDistribution ────────
# Permite usar toda la ecosistema Distributions.jl (pdf, cdf, fit, etc.)

function Distributions.params(d::PuigDistribution)
    return (d.λ, d.k, d.T)
end

Distributions.minimum(d::PuigDistribution) = 0.0
Distributions.maximum(d::PuigDistribution) = Inf
Distributions.insupport(d::PuigDistribution, x::Real) = x >= 0

Distributions.pdf(d::PuigDistribution, x::Real) = Puig_pdf(Float64(x), d; method=:asymp)
Distributions.logpdf(d::PuigDistribution, x::Real) = Puig_logpdf(x, d)

function Distributions.pdf(d::PuigDistribution, x::AbstractVector{<:Real})
    return Puig_pdf(x, d; method=:asymp)
end

function Distributions.logpdf(d::PuigDistribution, x::AbstractVector{<:Real})
    return Puig_logpdf(collect(Float64, x), d)
end

# ─── Acumulada y supervivencia ──────────────────────────────────────────────

Distributions.cdf(d::PuigDistribution, x::Real) = Puig_cumulative([Float64(x)], d)[1]
Distributions.ccdf(d::PuigDistribution, x::Real) = Puig_surviving([Float64(x)], d)[1]

function Distributions.cdf(d::PuigDistribution, x::AbstractVector{<:Real})
    return Puig_cumulative(collect(Float64, x), d)
end

function Distributions.ccdf(d::PuigDistribution, x::AbstractVector{<:Real})
    return Puig_surviving(collect(Float64, x), d)
end

# ─── Cuantiles ──────────────────────────────────────────────────────────────

Distributions.quantile(d::PuigDistribution, p::Real) = Puig_quantile(p, d)
Distributions.quantile(d::PuigDistribution, p::AbstractVector{<:Real}) = Puig_quantile(p, d)

# ─── Momentos y estadísticos ────────────────────────────────────────────────

Distributions.mean(d::PuigDistribution) = Puig_mean(d)
Distributions.var(d::PuigDistribution) = Puig_var(d)
Distributions.std(d::PuigDistribution) = Puig_std(d)
Distributions.skewness(d::PuigDistribution) = Puig_skewness(d)
Distributions.kurtosis(d::PuigDistribution) = Puig_kurtosis(d)

# ─── Muestreo ───────────────────────────────────────────────────────────────

function Distributions._rand!(rng::AbstractRNG, d::PuigDistribution, A::AbstractVector{Float64})
    for i in eachindex(A)
        A[i] = rand(rng, d)
    end
    return A
end

# ─── Soporte para el módulo Puig ────────────────────────────────────────────
# Las funciones Puig_* ya están definidas en puig_pdf.jl y operan sobre
# PuigDistribution. Las funciones Distributions.* despachan a ellas.

# Convenience: PuigDistribution(λ, k, T) ya está definido en puig_pdf.jl.
# Aquí no necesitamos constructor adicional porque el struct ya tiene el inner
# constructor que valida y promueve tipos.

# Soporte para advertencias silenciosas en marcumq (NoncentralChisq)
function Base.show(io::IO, ::MIME"text/plain", d::PuigDistribution)
    print(io, "PuigDistribution(λ=$(d.λ), k=$(d.k), T=$(d.T))")
end

# ─── Ajuste según interfaz Distributions.jl (fit) ──────────────────────────
function Distributions.fit(::Type{PuigDistribution}, data::Union{DataFrame, AbstractMatrix{<:Real}, AbstractVector{<:Real}}; kwargs...)
    return Puig_fit(data; kwargs...).params
end
