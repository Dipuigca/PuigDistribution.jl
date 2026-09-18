pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter
using DistributionsPuig

makedocs(
    sitename = "DistributionsPuig.jl",
    modules  = [DistributionsPuig],
    authors  = "Diego Puig",
    remotes  = nothing,
    warnonly = true,
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true,
    ),
    pages = [
        "Inicio" => "index.md",
        "Fundamentación Teórica" => "theory.md",
        "Tutorial y Ejemplos" => "tutorial.md",
        "Referencia de la API" => "api.md"
    ]
)
