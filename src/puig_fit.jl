## Estructura para las Estadísticas / Métricas de Ajuste

"""
    struct PuigStats{S<:Real}

Almacena las métricas de bondad de ajuste del modelo frente a la supervivencia empírica:
- `R²`: Coeficiente de determinación sobre la supervivencia.
- `MAE`: Error absoluto medio.
- `RMSE`: Raíz del error cuadrático medio.
- `MaxAE`: Error absoluto máximo.
- `IAE`: Error absoluto integrado (regla del trapecio).
"""
struct PuigStats{S<:Real}
    R²::S
    MAE::S
    RMSE::S
    MaxAE::S
    IAE::S

    function PuigStats(R²::S, MAE::S, RMSE::S, MaxAE::S, IAE::S) where {S<:Real}
        new{S}(R², MAE, RMSE, MaxAE, IAE)
    end
end

PuigStats(R²::Real, MAE::Real, RMSE::Real, MaxAE::Real, IAE::Real) = 
    PuigStats(promote(Float64(R²), Float64(MAE), Float64(RMSE), Float64(MaxAE), Float64(IAE))...)

function Base.show(io::IO, s::PuigStats)
    println(io, "Ajuste PuigStats:")
    println(io, "  R²    = $(round(s.R², digits=5))")
    println(io, "  MAE   = $(round(s.MAE, digits=5))")
    println(io, "  RMSE  = $(round(s.RMSE, digits=5))")
    println(io, "  MaxAE = $(round(s.MaxAE, digits=5))")
    print(io,   "  IAE   = $(round(s.IAE, digits=5))")
end


## Estructura del Resultado del Ajuste

"""
    struct PuigFitResult{S<:Real}

Contenedor del resultado del ajuste:
- `params::PuigDistribution{S}`: Caracterización de la distribución de Puig (λ, k, T).
- `indicator::Vector{S}`: Vector con los valores del indicador compuesto (normas/tiempos).
- `stats::PuigStats{S}`: Métricas de bondad de ajuste (R², MAE, RMSE, MaxAE, IAE).

Soporta acceso por propiedades (`res.params`, `res.indicator`, `res.stats`) 
y desempaquetado directo: `dist, ind, st = Puig_fit(df)`.
"""
struct PuigFitResult{S<:Real}
    params::PuigDistribution{S}
    indicator::Vector{S}
    stats::PuigStats{S}
end

function Base.show(io::IO, r::PuigFitResult)
    println(io, "PuigFitResult:")
    println(io, "  Distribución (PuigDistribution): λ = $(round(r.params.λ, digits=4)), k = $(round(r.params.k, digits=4)), T = $(round(r.params.T, digits=5))")
    println(io, "  Indicador: $(length(r.indicator)) observaciones (rango: [$(round(minimum(r.indicator), digits=3)), $(round(maximum(r.indicator), digits=3))])")
    println(io, "  Métricas (PuigStats):")
    println(io, "    R²    = $(round(r.stats.R², digits=5))")
    println(io, "    MAE   = $(round(r.stats.MAE, digits=5))")
    println(io, "    RMSE  = $(round(r.stats.RMSE, digits=5))")
    println(io, "    MaxAE = $(round(r.stats.MaxAE, digits=5))")
    print(io,   "    IAE   = $(round(r.stats.IAE, digits=5))")
end

# Soporte para desempaquetado: params, indicator, fit_stats = Puig_fit(...)
Base.iterate(r::PuigFitResult, state=1) = state == 1 ? (r.params, 2) : (state == 2 ? (r.indicator, 3) : (state == 3 ? (r.stats, 4) : nothing))
Base.length(::PuigFitResult) = 3
Base.firstindex(::PuigFitResult) = 1
Base.lastindex(::PuigFitResult) = 3
Base.getindex(r::PuigFitResult, i::Int) = (r.params, r.indicator, r.stats)[i]


## Métricas de bondad de ajuste (evaluadas sobre la supervivencia)

calc_r2(S_emp::AbstractVector{<:Real}, S_mod::AbstractVector{<:Real}) = 
    1.0 - sum((S_emp .- S_mod).^2) / sum((S_emp .- mean(S_emp)).^2)

calc_mae(S_emp::AbstractVector{<:Real}, S_mod::AbstractVector{<:Real}) = 
    mean(abs.(S_emp .- S_mod))

calc_rmse(S_emp::AbstractVector{<:Real}, S_mod::AbstractVector{<:Real}) = 
    sqrt(mean((S_emp .- S_mod).^2))

calc_maxae(S_emp::AbstractVector{<:Real}, S_mod::AbstractVector{<:Real}) = 
    maximum(abs.(S_emp .- S_mod))

function calc_iae(x::AbstractVector{<:Real}, S_emp::AbstractVector{<:Real}, S_mod::AbstractVector{<:Real})
    diffs = abs.(S_emp .- S_mod)
    dx = diff(x)
    return sum(dx .* (diffs[1:end-1] .+ diffs[2:end]) ./ 2.0)
end

function calc_all_metrics(x::AbstractVector{<:Real}, S_emp::AbstractVector{<:Real}, S_mod::AbstractVector{<:Real})::PuigStats{Float64}
    return PuigStats(
        calc_r2(S_emp, S_mod),
        calc_mae(S_emp, S_mod),
        calc_rmse(S_emp, S_mod),
        calc_maxae(S_emp, S_mod),
        calc_iae(x, S_emp, S_mod)
    )
end


## Dimensión efectiva a partir de la matriz de correlación

"""
    effective_dimension(M::AbstractMatrix{<:Real}) -> Float64

Calcula la dimensión efectiva (razón de participación) a partir de la matriz de correlación R:
    k = (Σ λᵢ)² / Σ (λᵢ²) = d² / ‖R‖_F² = d² / (d + 2 ∑_{i < j} R_{ij}²)
donde k ∈ [1, d] (d = número de columnas/variables).

Aprovecha la identidad de Frobenius y la traza constante de la matriz de correlación 
para evitar el cálculo de autovalores, reduciendo el coste a O(d²) con 0 alocaciones.
"""
function effective_dimension(M::AbstractMatrix{<:Real})::Float64
    _, d = size(M)
    d <= 1 && return 1.0
    R = cor(M)
    
    s_off = 0.0
    @inbounds for j in 2:d
        for i in 1:(j-1)
            s_off += R[i, j]^2
        end
    end
    
    denom = d + 2.0 * s_off
    return denom == 0.0 ? 1.0 : (Float64(d)^2 / denom)
end


## Función de supervivencia empírica

"""
    empirical_survival(times::AbstractVector{<:Real}, grid::AbstractVector{<:Real}) -> Vector{Float64}

Calcula la supervivencia empírica S_emp(x) = P(X >= x) evaluada sobre `grid`.
"""
function empirical_survival(times::AbstractVector{<:Real}, grid::AbstractVector{<:Real})::Vector{Float64}
    n = length(times)
    n == 0 && return zeros(Float64, length(grid))
    return [count(>=(g), times) / n for g in grid]
end


# Ajuste principal desde Matrix

"""
    Puig_fit(Data::AbstractMatrix{<:Real}; 
             times::Union{Nothing, AbstractVector{<:Real}} = nothing, 
             method::Symbol = :dynamic) -> PuigFitResult

Ajusta los datos multivariantes (filas = observaciones, columnas = variables) a la distribución de Puig.

# Métodos (`method`):
- `:dynamic` (recomendado): k fijado por la dimensión efectiva de la correlación (k_corr) y B estimado minimizando la distancia L1.
- `:numerical`: k y B optimizados conjuntamente mediante L-BFGS-B minimizando el error cuadrático L2.
- `:fixed`: k fijado a la dimensión geométrica d (columnas) y B estimado minimizando L1.
"""
function Puig_fit(Data::AbstractMatrix{<:Real}; 
                  times::Union{Nothing, AbstractVector{<:Real}} = nothing, 
                  method::Symbol = :dynamic)::PuigFitResult
    
    n_obs, d = size(Data)
    n_obs >= 10 || throw(ArgumentError("Se requieren al menos 10 observaciones para el ajuste (recibidas: $n_obs)"))

    # 1. Vector de tiempos / norma euclídea por fila si no se especifica
    t_raw = times === nothing ? sqrt.(sum(Data .^ 2, dims=2)[:, 1]) : collect(Float64, times)
    valid_idx = findall(t -> isfinite(t) && t >= 0, t_raw)
    t = t_raw[valid_idx]
    
    # 2. λ = norma del vector de medias de las componentes
    col_means = [mean(Data[valid_idx, j]) for j in 1:d]
    λ = sqrt(sum(col_means .^ 2))

    # 3. Dimensión efectiva k_corr
    k_corr = effective_dimension(Data[valid_idx, :])

    # Rejilla para optimización y supervivencia empírica
    min_t, max_t = extrema(t)
    x_opt = collect(range(min_t, max_t, length=150))
    S_emp_opt = empirical_survival(t, x_opt)

    sdt = std(t)
    B0 = (isfinite(sdt) && sdt > 0) ? sdt : 1.0

    k_est = k_corr
    B_est = B0

    if method == :dynamic || method == :fixed
        k_target = (method == :fixed) ? Float64(d) : k_corr

        # Minimización L1 sobre B ∈ [0.01, 5·B0]
        err_l1 = function(B)
            B <= 0 && return 1e12
            T_val = 1.0 / (B^2)
            try
                S_mod = Puig_surviving(x_opt, λ, k_target, T_val)
                return sum(abs.(S_emp_opt .- S_mod))
            catch
                return 1e12
            end
        end

        res = Optim.optimize(err_l1, 0.01, 5.0 * B0, Brent())
        B_est = Optim.minimizer(res)
        k_est = k_target

    elseif method == :numerical
        # Minimización L2 conjunta (k, B) mediante L-BFGS-B
        err_l2 = function(par)
            k_val, B_val = par[1], par[2]
            (k_val < 1.0 || B_val <= 0.0) && return 1e12
            T_val = 1.0 / (B_val^2)
            try
                S_mod = Puig_surviving(x_opt, λ, k_val, T_val)
                any(!isfinite, S_mod) && return 1e12
                return sum((S_emp_opt .- S_mod) .^ 2)
            catch
                return 1e12
            end
        end

        lower = [1.0, 0.001]
        upper = [max(3.0, Float64(d)), 10.0 * B0]
        init_par = [clamp(k_corr, 1.0, upper[1]), B0]

        res = try
            Optim.optimize(err_l2, lower, upper, init_par, Fminbox(LBFGS()))
        catch
            nothing
        end

        if res !== nothing && Optim.converged(res)
            k_est, B_est = Optim.minimizer(res)
        else
            # Fallback a dinámico
            k_est = clamp(k_corr, 1.0, Float64(d))
            B_est = B0
        end
    else
        throw(ArgumentError("Método desconocido :$method. Opciones válidas: :dynamic, :numerical, :fixed"))
    end

    # 4. Cálculo de T final y métricas sobre rejilla fina
    T_est = 1.0 / (B_est^2)
    x_eval = collect(range(min_t, max_t, length=500))
    S_emp_eval = empirical_survival(t, x_eval)
    S_mod_eval = Puig_surviving(x_eval, λ, k_est, T_est)
    fit_metrics = calc_all_metrics(x_eval, S_emp_eval, S_mod_eval)

    params = PuigDistribution(Float64(λ), Float64(k_est), Float64(T_est))
    return PuigFitResult(params, t, fit_metrics)
end


## Ajuste desde DataFrame

"""
    Puig_fit(df::DataFrame; 
             vars::Union{Nothing, Vector{Symbol}, Vector{String}} = nothing, 
             time_col::Union{Nothing, Symbol, String} = nothing, 
             method::Symbol = :dynamic) -> PuigFitResult

Ajusta las columnas de un `DataFrame` a la distribución de Puig.

# Argumentos:
- `vars`: Columnas que componen las dimensiones/disciplinas (por defecto todas las numéricas excepto `time_col`).
- `time_col`: Columna opcional que ya contiene la norma/tiempo calculada. Si no se pasa, se calcula como `√(∑ vars²)`.
- `method`: `:dynamic` (defecto), `:numerical`, o `:fixed`.
"""
function Puig_fit(df::DataFrame; 
                  vars::Union{Nothing, Vector{Symbol}, Vector{String}} = nothing, 
                  time_col::Union{Nothing, Symbol, String} = nothing, 
                  method::Symbol = :dynamic)::PuigFitResult
    
    # 1. Identificar columnas de variables
    var_names = if vars !== nothing
        Symbol.(vars)
    else
        # Si no se pasan, tomar todas las columnas numéricas excluyendo time_col
        exclude = time_col === nothing ? Symbol[] : [Symbol(time_col)]
        [col for col in propertynames(df) if eltype(df[!, col]) <: Real && !(col in exclude)]
    end

    length(var_names) >= 1 || throw(ArgumentError("Se requiere al menos una columna numérica para el ajuste"))

    # 2. Convertir a matriz
    mat = Matrix{Float64}(df[:, var_names])

    # 3. Vector de tiempos si existe
    times = time_col !== nothing ? Vector{Float64}(df[!, Symbol(time_col)]) : nothing

    return Puig_fit(mat; times=times, method=method)
end


# Ajuste 1D / Vectorial (Univariante)

"""
    Puig_fit(Data::AbstractVector{<:Real}; 
             method::Symbol = :dynamic, 
             k_fixed::Union{Nothing, Real} = nothing) -> PuigFitResult

Ajusta una muestra 1D de observaciones (normas observadas o indicador compuesto) a la distribución de Puig estimando conjuntamente los parámetros (λ, k, T).

# Argumentos:
- `Data`: Vector con las observaciones unidimensionales.
- `method`:
  - `:dynamic` (por defecto) o `:numerical`: Optimización conjunta de los 3 parámetros (λ, k, B) minimizando el error cuadrático L2 frente a la supervivencia empírica.
  - `:fixed`: Fija k en `k_fixed` (o k=1.0 si no se especifica) y optimiza conjuntamente (λ, B).
- `k_fixed`: Dimensión conocida a priori (ej. k=3.0 para 3 disciplinas). Si se proporciona, fuerza el comportamiento con k fijo.
"""
function Puig_fit(Data::AbstractVector{<:Real}; 
                  method::Symbol = :dynamic, 
                  k_fixed::Union{Nothing, Real} = nothing)::PuigFitResult
    
    valid_idx = findall(x -> isfinite(x) && x >= 0, Data)
    t = Float64.(Data[valid_idx])
    n_obs = length(t)
    n_obs >= 10 || throw(ArgumentError("Se requieren al menos 10 observaciones para el ajuste (recibidas: $n_obs)"))

    # 1. Rejilla y supervivencia empírica
    min_t, max_t = extrema(t)
    x_opt = collect(range(min_t, max_t, length=60))
    S_emp_opt = empirical_survival(t, x_opt)

    # 2. Inicialización robusta por momentos
    m1 = mean(t)
    m2 = mean(t .^ 2)
    s = std(t)
    B0 = (isfinite(s) && s > 0) ? s : 1.0
    λ0 = sqrt(max(0.0, m1^2 - s^2))
    k0 = clamp((m2 - λ0^2) / B0^2, 1.0, 10.0)

    is_fixed_k = (method == :fixed) || (k_fixed !== nothing)
    k_target = k_fixed !== nothing ? Float64(k_fixed) : 1.0
    k_target >= 1.0 || throw(ArgumentError("k_fixed debe ser >= 1.0 (recibido: $k_target)"))

    k_est, λ_est, B_est = 1.0, λ0, B0

    if is_fixed_k
        k_est = k_target
        # Optimización 2D sobre (u=log λ, w=log B)
        err_2d = function(p)
            λ_val = exp(p[1])
            B_val = exp(p[2])
            T_val = 1.0 / (B_val^2)
            try
                S_mod = Puig_surviving(x_opt, λ_val, k_target, T_val)
                any(!isfinite, S_mod) && return 1e12
                return sum((S_emp_opt .- S_mod) .^ 2)
            catch
                return 1e12
            end
        end

        init_p = [log(max(1e-3, λ0)), log(max(1e-3, B0))]
        res = Optim.optimize(err_2d, init_p, NelderMead(), Optim.Options(iterations=500))
        p_opt = Optim.minimizer(res)
        λ_est = exp(p_opt[1])
        B_est = exp(p_opt[2])

    elseif method == :dynamic || method == :numerical
        # Optimización 3D conjunta sobre (u=log λ, v=log(k-1), w=log B)
        err_3d = function(p)
            λ_val = exp(p[1])
            k_val = 1.0 + exp(p[2])
            B_val = exp(p[3])
            T_val = 1.0 / (B_val^2)
            try
                S_mod = Puig_surviving(x_opt, λ_val, k_val, T_val)
                any(!isfinite, S_mod) && return 1e12
                return sum((S_emp_opt .- S_mod) .^ 2)
            catch
                return 1e12
            end
        end

        init_p = [log(max(1e-3, λ0)), log(max(1e-3, k0 - 1.0)), log(max(1e-3, B0))]
        res = Optim.optimize(err_3d, init_p, NelderMead(), Optim.Options(iterations=600))
        p_opt = Optim.minimizer(res)
        λ_est = exp(p_opt[1])
        k_est = 1.0 + exp(p_opt[2])
        B_est = exp(p_opt[3])
    else
        throw(ArgumentError("Método desconocido :$method. Opciones válidas: :dynamic, :numerical, :fixed"))
    end

    # 3. Evaluación sobre rejilla fina y empaquetado
    T_est = 1.0 / (B_est^2)
    x_eval = collect(range(min_t, max_t, length=500))
    S_emp_eval = empirical_survival(t, x_eval)
    S_mod_eval = Puig_surviving(x_eval, λ_est, k_est, T_est)
    fit_metrics = calc_all_metrics(x_eval, S_emp_eval, S_mod_eval)

    params = PuigDistribution(Float64(λ_est), Float64(k_est), Float64(T_est))
    return PuigFitResult(params, t, fit_metrics)
end
