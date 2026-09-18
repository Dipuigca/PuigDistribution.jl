# ============================================================================
# Genera los "golden files" de referencia para la validación cruzada de los
# portes a Python, R y MATLAB.
#
# Uso:  julia --project=.. golden/generate.jl
# (Desde la raíz del repo:  julia --project=. golden/generate.jl )
#
# Emite en ./golden:
#   golden_pdf.csv          λ,k,T,x,pdf_asymp,pdf_arb
#   golden_surv.csv        λ,k,T,x,su
#   golden_cdf.csv         λ,k,T,x,cu
#   golden_moments.csv     λ,k,T,n,mu
#   golden_quantile.csv    λ,k,T,p,q
#   golden_mb.csv          k,B,x,mb_pdf,mb_surv,mb_cdf
#   golden_mb_moments.csv  k,B,n,mb_mu
# ============================================================================

using DistributionsPuig

const OUT = joinpath(@__DIR__)

# Parámetros de referencia  (el fichero PATH incluye rangos extremos y colas)
const LAM = [0.0, 0.5, 3.0, 10.0, 50.0, 100.0, 250.0, 1000.0]
const K   = [1.0001, 1.5, 2.4, 3.0, 10.0]
const T   = [0.001, 0.01, 0.05, 0.25, 1.0]
const PDF_FRACS = [0.7, 0.9, 1.0, 1.1, 1.3]     # x = frac·λ  (λ>0)
const SURV_FRACS = [0.2, 0.8, 1.0, 1.5, 3.0]    # x = frac·λ  (λ>0)
const X0 = [0.05, 0.5, 1.0, 2.0, 5.0]           # x cuando λ = 0

const QP = [0.01, 0.25, 0.5, 0.75, 0.95, 0.99]

# Combos de momentos Puig (λ>0)
const MOM_COMBOS = [
    (0.5, 1.5, 0.05),
    (3.0, 1.0001, 0.05),
    (3.0, 2.4, 0.25),
    (10.0, 3.0, 0.05),
    (50.0, 2.4, 0.01),
    (100.0, 2.5, 0.05),
    (250.0, 2.4, 0.01),
    (1000.0, 3.0, 0.05),
]

# Combos MB
const MB_COMBOS = [(1.5, 1.0), (2.4, 10.0), (3.0, 5.0), (5.0, 2.0), (10.0, 0.5)]

fmt(v::Real) = repr(Float64(v))

function write_csv(name::String, header::Tuple, rows)
    open(joinpath(OUT, name), "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join(r, ","))
        end
    end
    println("✓ ", name, "  ", length(rows), " filas")
    return nothing
end

# ----- PDF ----------------------------------------------------------------
function gen_pdf()
    rows = Tuple{String,String,String,String,String,String}[]
    for λ in LAM, kv in K, Tv in T
        if λ == 0
            xs = X0
        else
            xs = PDF_FRACS .* λ
        end
        for x in xs
            a = Puig_pdf(x, λ, kv, Tv; method=:asymp)
            b = λ == 0 ? a : Puig_pdf(x, λ, kv, Tv; method=:arb)
            push!(rows, (fmt(λ), fmt(kv), fmt(Tv), fmt(x), fmt(a), fmt(b)))
        end
    end
    write_csv("golden_pdf.csv", ("λ", "k", "T", "x", "pdf_asymp", "pdf_arb"), rows)
end

# ----- Supervivencia / Acumulada -------------------------------------------
function gen_surv()
    rows = Tuple{String,String,String,String,String}[]
    for λ in LAM, kv in K, Tv in T
        xs = λ == 0 ? X0 : SURV_FRACS .* λ
        for x in xs
            su = Puig_surviving([x], λ, kv, Tv)[1]
            cu = Puig_cumulative([x], λ, kv, Tv)[1]
            push!(rows, (fmt(λ), fmt(kv), fmt(Tv), fmt(x), fmt(su)))
        end
    end
    write_csv("golden_surv.csv", ("λ", "k", "T", "x", "su"), rows)

    rows2 = Tuple{String,String,String,String,String}[]
    for λ in LAM, kv in K, Tv in T
        xs = λ == 0 ? X0 : SURV_FRACS .* λ
        for x in xs
            cu = Puig_cumulative([x], λ, kv, Tv)[1]
            push!(rows2, (fmt(λ), fmt(kv), fmt(Tv), fmt(x), fmt(cu)))
        end
    end
    write_csv("golden_cdf.csv", ("λ", "k", "T", "x", "cu"), rows2)
end

# ----- Momentos Puig --------------------------------------------------------
function gen_moments()
    rows = Tuple{String,String,String,String,String}[]
    for (λ, kv, Tv) in MOM_COMBOS, n in 1:8
        μ = Puig_moments(λ, kv, Tv, n)
        push!(rows, (fmt(λ), fmt(kv), fmt(Tv), string(n), fmt(μ[n])))
    end
    write_csv("golden_moments.csv", ("λ", "k", "T", "n", "mu"), rows)
end

# ----- Cuantiles ------------------------------------------------------------
function gen_quantile()
    rows = Tuple{String,String,String,String,String}[]
    for (λ, kv, Tv) in MOM_COMBOS, p in QP
        q = Puig_quantile(p, λ, kv, Tv)
        push!(rows, (fmt(λ), fmt(kv), fmt(Tv), fmt(p), fmt(q)))
    end
    write_csv("golden_quantile.csv", ("λ", "k", "T", "p", "q"), rows)
end

# ----- Maxwell-Boltzmann ----------------------------------------------------
function gen_mb()
    rows = Tuple{String,String,String,String,String,String}[]
    for (kv, Bv) in MB_COMBOS, x in X0
        p = MB_pdf([x], kv, Bv)[1]
        s = MB_surviving([x], kv, Bv)[1]
        c = MB_cumulative([x], kv, Bv)[1]
        push!(rows, (fmt(kv), fmt(Bv), fmt(x), fmt(p), fmt(s), fmt(c)))
    end
    write_csv("golden_mb.csv", ("k", "B", "x", "mb_pdf", "mb_surv", "mb_cdf"), rows)

    rows2 = Tuple{String,String,String,String}[]
    for (kv, Bv) in MB_COMBOS, n in 1:8
        μ = Puig_moments(0.0, kv, 1.0 / Bv, n)
        push!(rows2, (fmt(kv), fmt(Bv), string(n), fmt(μ[n])))
    end
    write_csv("golden_mb_moments.csv", ("k", "B", "n", "mb_mu"), rows2)
end

println("Generando golden files de DistributionsPuig.jl …")
gen_pdf()
gen_surv()
gen_moments()
gen_quantile()
gen_mb()
println("Listo. Golden files en: ", OUT)