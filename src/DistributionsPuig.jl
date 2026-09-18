module DistributionsPuig

using ArbNumerics
using Distributions
using SpecialFunctions
using Random
using Statistics
using LinearAlgebra
using DataFrames
using Optim
using Plots

include("Maxwell_Boltzmann.jl")
include("puig_pdf.jl")
include("types.jl")
include("puig_fit.jl")
include("puig_plot.jl")

# ─── API vectorizada (custom) ──────────────────────────────────────────────
export PuigDistribution, MB
export Puig_pdf, Puig_pdf_vec, Puig_pdf_vec1
export Puig_surviving, Puig_cumulative, Puig_logpdf
export Puig_mean, Puig_var, Puig_std, Puig_skewness, Puig_kurtosis
export Puig_moments, Puig_stats, Puig_hazard
export Puig_quantile, Puig_ci98, Puig_ci996, Puig_entropy
export Puig_plot
export Puig_fit, PuigFitResult, PuigStats
export MB_pdf, MB_surviving, MB_cumulative

# ─── Funciones auxiliares internas (exportadas para testing) ────────────────
export effective_dimension, empirical_survival
export calc_r2, calc_mae, calc_rmse, calc_maxae, calc_iae, calc_all_metrics
export marcumq, marcumq_asymp, besselix_asymp

# ─── API Distributions.jl (no re-exportar: el usuario carga Distributions) ─
# Las siguientes funciones están definidas en src/types.jl y se dispatchan
# automáticamente sobre PuigDistribution cuando el usuario hace `using Distributions`:
#   Distributions.pdf, Distributions.cdf, Distributions.ccdf,
#   Distributions.logpdf, Distributions.quantile, Distributions.rand,
#   Distributions.mean, Distributions.var, Distributions.std,
#   Distributions.skewness, Distributions.kurtosis,
#   Distributions.params, Distributions.minimum, Distributions.maximum,
#   Distributions.insupport

end
