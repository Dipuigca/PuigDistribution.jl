# --- helpers internos: un x-grid razonable y cada subplot por separado ---

function _puig_xgrid(params::PuigDistribution; npoints::Int=400, pad_frac::Float64=0.05)
    x_min, x_max = Puig_ci996(params)
    margen = pad_frac * (x_max - x_min)
    return collect(range(max(0.0, x_min - margen), x_max, length=npoints))
end
function _plot_pdf(params::PuigDistribution)
    x = _puig_xgrid(params)
    plot(x, Puig_pdf(x, params), label="pdf", xlabel="x", ylabel="f(x)",
         title="Densidad (pdf)", legend=false)
end
function _plot_cdf(params::PuigDistribution)
    x = _puig_xgrid(params)
    plot(x, Puig_cumulative(x, params), label="cdf", xlabel="x", ylabel="F(x)",
         title="Acumulada (cdf)", legend=false)
end
function _plot_surv(params::PuigDistribution)
    x = _puig_xgrid(params)
    plot(x, Puig_surviving(x, params), label="S(x)", xlabel="x", ylabel="S(x)",
         title="Supervivencia", legend=false)
end
function _plot_data(params::PuigDistribution)
    st = Puig_stats(params)
    Limi, Limsu = Puig_ci98(params)
    q1, q2, q3 = Puig_quantile([0.25, 0.5, 0.75], params)

    texto = """
    λ = $(round(params.λ, digits=4))   k = $(round(params.k, digits=4))   T = $(round(params.T, digits=4))

    media     = $(round(st.mean, digits=4))
    varianza  = $(round(st.var, digits=4))
    σ        = $(round(st.sig, digits=4))
    skewness  = $(round(st.skewness, digits=4))
    kurtosis  = $(round(st.kurtosis, digits=4))

    P1  = $(round(Limi, digits=4))    P99 = $(round(Limsu, digits=4))

    Q1 (25%) = $(round(q1, digits=4))
    Q2 (50%) = $(round(q2, digits=4))
    Q3 (75%) = $(round(q3, digits=4))
    """

    p = plot(framestyle=:none, legend=false, title="Parámetros y estadísticos")
    annotate!(p, 0.0, 0.5, text(texto, :left, 9))
    return p
end
const _PUIG_PLOT_FUNS = Dict(
    "pdf"  => _plot_pdf,
    "cdf"  => _plot_cdf,
    "surv" => _plot_surv,
    "data" => _plot_data,
)
# --- API pública ---

"""
    Puig_plot(params::Union{PuigDistribution, MB}, options::String="all") -> Plots.Plot
    Puig_plot(params::Union{PuigDistribution, MB}, options::Tuple{Vararg{String}}) -> Plots.Plot
    Puig_plot(λ::Real, k::Real, T::Real, options="all") -> Plots.Plot

Genera representaciones gráficas de la distribución de Puig o Maxwell-Boltzmann usando `Plots.jl`.

# Opciones (`options`):
- `"all"` (por defecto): Panel comparativo 2x2 conteniendo Densidad (PDF), Acumulada (CDF), Supervivencia (CCDF) y Cuadro de Estadísticos y Cuantiles.
- `"pdf"`: Gráfico individual de la función de densidad f(x).
- `"cdf"`: Gráfico individual de la función de distribución acumulada F(x).
- `"surv"`: Gráfico individual de la función de supervivencia S(x).
- `"data"`: Panel textual con parámetros, momentos (media, varianza, asimetría, kurtosis) y cuantiles (P1, Q1, Q2, Q3, P99).
- Tupla de opciones (ej. `("pdf", "surv")`): Genera un layout con la combinación exacta solicitada.
"""
function Puig_plot(params::PuigDistribution, options::String="all")
    if options == "all"
        return plot(_plot_pdf(params), _plot_cdf(params), _plot_surv(params), _plot_data(params),
                    layout=(2,2), size=(900,700))
    end
    haskey(_PUIG_PLOT_FUNS, options) || throw(ArgumentError(
        "options debe ser \"all\", \"pdf\", \"cdf\", \"surv\" o \"data\"; recibido \"$options\""))
    return _PUIG_PLOT_FUNS[options](params)
end
function Puig_plot(params::PuigDistribution, options::Tuple{Vararg{String}})
    "all" in options && throw(ArgumentError("\"all\" no está permitido dentro de una tupla de opciones"))
    all(o -> haskey(_PUIG_PLOT_FUNS, o), options) || throw(ArgumentError(
        "Todas las opciones de la tupla deben ser \"pdf\", \"cdf\", \"surv\" o \"data\"; recibido $options"))

    plots = [_PUIG_PLOT_FUNS[o](params) for o in options]
    n = length(plots)
    ncols = n <= 2 ? n : ceil(Int, sqrt(n))
    nrows = ceil(Int, n/ncols)
    return plot(plots..., layout=(nrows, ncols), size=(450*ncols, 350*nrows))
end
Puig_plot(params::MB, options="all") = Puig_plot(PuigDistribution(0.0, Float64(params.k), 1.0/Float64(params.B)), options)
function Puig_plot(λ::Real,k::Real,T::Real, options::Union{String,Tuple{Vararg{String}}}="all")
    χ=PuigDistribution(λ,k,T);
    return Puig_plot(χ,options)
end
