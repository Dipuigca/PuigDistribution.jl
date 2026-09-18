using Test
using DistributionsPuig
using Distributions

@testset "Suite de Tests: DistributionsPuig.jl" begin

    # ============================================================
    # 0. Tests de carga del módulo y tipos
    # ============================================================
    @testset "Módulo y tipos" begin
        d = PuigDistribution(3.0, 2.0, 1.0)
        @test d isa PuigDistribution
        @test d.λ == 3.0
        @test d.k == 2.0
        @test d.T == 1.0

        mb = MB(2.0, 1.0)
        @test mb isa MB
        @test mb.k == 2.0
        @test mb.B == 1.0
    end

    # ============================================================
    # 1. Tests para effective_dimension
    # ============================================================
    @testset "effective_dimension" begin
        using LinearAlgebra, Random

        # Caso 1D
        M1 = randn(50, 1)
        @test effective_dimension(M1) ≈ 1.0

        # Caso Variables Incorreladas (R ≈ I_3 -> k ≈ 3)
        Random.seed!(123)
        M_uncorr = randn(10000, 3)
        @test effective_dimension(M_uncorr) ≈ 3.0 atol=0.1

        # Caso Variables Perfectamente Correladas (k = 1)
        x_base = randn(100)
        M_corr = hcat(x_base, x_base, x_base)
        @test effective_dimension(M_corr) ≈ 1.0 atol=1e-5

        # Verificación exacta con fórmula conocida
        d = 3
        s_off = 0.8^2 + 0.6^2 + 0.7^2
        k_exacto = Float64(d)^2 / (d + 2.0 * s_off)
        @test k_exacto ≈ 1.5050167224 atol=1e-6
    end

    # ============================================================
    # 2. Tests para empirical_survival
    # ============================================================
    @testset "empirical_survival" begin
        t = [1.0, 2.0, 3.0, 4.0, 5.0]
        grid = [0.0, 1.0, 3.0, 5.0, 6.0]
        s = empirical_survival(t, grid)
        
        @test s[1] == 1.0  # Todos >= 0.0
        @test s[2] == 1.0  # Todos >= 1.0
        @test s[3] == 0.6  # [3, 4, 5] -> 3/5 = 0.6
        @test s[4] == 0.2  # [5] -> 1/5 = 0.2
        @test s[5] == 0.0  # Ninguno >= 6.0

        # Monotonía no creciente
        @test issorted(s, rev=true)
    end

    # ============================================================
    # 3. Tests para métricas y struct PuigStats
    # ============================================================
    @testset "Métricas y struct PuigStats" begin
        using Statistics

        # Ajuste perfecto
        S_true = [1.0, 0.8, 0.5, 0.2, 0.0]
        x_grid = [1.0, 2.0, 3.0, 4.0, 5.0]
        st_perf = calc_all_metrics(x_grid, S_true, S_true)
        
        @test st_perf isa PuigStats{Float64}
        @test st_perf.R² ≈ 1.0
        @test st_perf.MAE ≈ 0.0
        @test st_perf.RMSE ≈ 0.0
        @test st_perf.MaxAE ≈ 0.0
        @test st_perf.IAE ≈ 0.0

        # Discrepancia controlada
        S_approx = S_true .+ 0.05
        st_err = calc_all_metrics(x_grid, S_true, S_approx)
        @test st_err.MAE ≈ 0.05
        @test st_err.MaxAE ≈ 0.05
        @test st_err.RMSE ≈ 0.05
        @test st_err.R² < 1.0
    end

    # ============================================================
    # 4. Tests de ajuste con DataFrame y Matrices
    # ============================================================
    @testset "Puig_fit sobre DataFrame y Matrix" begin
        using DataFrames, Random, Statistics

        Random.seed!(42)
        N = 150
        X1 = 150.0 .+ randn(N) .* 20.0
        X2 = 100.0 .+ 0.8 .* X1 .+ randn(N) .* 10.0
        X3 = 200.0 .+ 0.5 .* X1 .+ randn(N) .* 15.0
        
        df = DataFrame(Squat=X1, Bench=X2, Deadlift=X3)

        # 4.1 Método dinámico (:dynamic)
        res_dyn = Puig_fit(df; method=:dynamic)
        @test res_dyn isa PuigFitResult
        @test res_dyn.params isa PuigDistribution
        @test res_dyn.indicator isa Vector{Float64}
        @test res_dyn.stats isa PuigStats
        @test length(res_dyn.indicator) == N
        @test res_dyn.params.λ > 0
        @test 1.0 <= res_dyn.params.k <= 3.0
        @test res_dyn.params.T > 0
        @test res_dyn.stats.R² > 0.95

        # 4.2 Método fijo (:fixed con k=3)
        res_fix = Puig_fit(df; method=:fixed)
        @test res_fix.params.k == 3.0
        @test res_fix.stats.R² > 0.90

        # 4.3 Método numérico (:numerical)
        res_num = Puig_fit(df; method=:numerical)
        @test 1.0 <= res_num.params.k <= 3.0
        @test res_num.stats.R² > 0.95

        # 4.4 Subconjunto de variables con `vars`
        res_sub = Puig_fit(df; vars=[:Squat, :Bench])
        @test res_sub.params.k <= 2.0
        @test length(res_sub.indicator) == N

        # 4.5 Con columna de tiempo / indicador precalculado
        df_with_time = copy(df)
        df_with_time.Total = sqrt.(df.Squat.^2 .+ df.Bench.^2 .+ df.Deadlift.^2)
        res_time = Puig_fit(df_with_time; vars=[:Squat, :Bench, :Deadlift], time_col=:Total)
        @test res_time.indicator == df_with_time.Total

        # 4.6 Ajuste directo desde Matrix
        mat = Matrix{Float64}(df)
        res_mat = Puig_fit(mat; method=:dynamic)
        @test res_mat.params.λ ≈ res_dyn.params.λ
        @test res_mat.params.k ≈ res_dyn.params.k
        @test res_mat.params.T ≈ res_dyn.params.T

        # 4.7 Ajuste desde Vector 1D (3 parámetros conjuntos)
        vec_ind = df_with_time.Total
        res_vec = Puig_fit(vec_ind)
        @test res_vec isa PuigFitResult
        @test length(res_vec.indicator) == N
        @test res_vec.params.λ > 0
        @test res_vec.params.k >= 1.0
        @test res_vec.params.T > 0
        @test res_vec.stats.R² > 0.95

        # 4.8 Ajuste 1D sobre muestra sintética con k_fixed
        Random.seed!(77)
        dist_sint = PuigDistribution(40.0, 2.5, 0.05)
        muestra_1d = rand(dist_sint, 500)
        
        res_1d_opt = Puig_fit(muestra_1d)
        @test res_1d_opt.params.λ ≈ 40.0 atol=5.0
        @test res_1d_opt.params.k >= 1.0
        @test res_1d_opt.stats.R² > 0.98

        res_1d_kfix = Puig_fit(muestra_1d; k_fixed=2.5)
        @test res_1d_kfix.params.k == 2.5
        @test res_1d_kfix.stats.R² > 0.98
    end

    # ============================================================
    # 5. Tests de Desempaquetado e Indexación
    # ============================================================
    @testset "Desempaquetado de PuigFitResult" begin
        using DataFrames, Random

        Random.seed!(99)
        df_test = DataFrame(A = 10.0 .+ randn(50), B = 20.0 .+ randn(50))
        res = Puig_fit(df_test)

        # Desempaquetado tipo tupla
        p, ind, st = res
        @test p === res.params
        @test ind === res.indicator
        @test st === res.stats

        # Indexación por posición
        @test res[1] === res.params
        @test res[2] === res.indicator
        @test res[3] === res.stats
        @test length(res) == 3
    end

    # ============================================================
    # 6. Tests de Validación y Manejo de Errores
    # ============================================================
    @testset "Manejo de Errores y Casos Límite" begin
        using DataFrames

        # Muestra demasiado pequeña (N < 10)
        df_small = DataFrame(A=rand(5), B=rand(5))
        @test_throws ArgumentError Puig_fit(df_small)
        @test_throws ArgumentError Puig_fit(rand(5))

        # Método de ajuste desconocido
        df_ok = DataFrame(A=rand(20), B=rand(20))
        @test_throws ArgumentError Puig_fit(df_ok; method=:metodo_inexistente)
        @test_throws ArgumentError Puig_fit(rand(20); method=:metodo_inexistente)

        # k_fixed inválido (< 1.0)
        @test_throws ArgumentError Puig_fit(rand(20); k_fixed=0.5)

        # DataFrame sin columnas numéricas
        df_non_numeric = DataFrame(Nombre=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"])
        @test_throws ArgumentError Puig_fit(df_non_numeric)
    end

    # ============================================================
    # 7. Tests de la interfaz Distributions.jl
    # ============================================================
    @testset "Interfaz Distributions.jl" begin
        d = PuigDistribution(3.0, 2.0, 1.0)

        # pdf / logpdf
        @test Distributions.pdf(d, 1.0) isa Float64
        @test Distributions.pdf(d, 0.0) == 0.0
        @test Distributions.logpdf(d, 1.0) isa Float64
        @test Distributions.logpdf(d, 0.0) == -Inf

        # cdf / ccdf
        @test Distributions.cdf(d, 0.0) == 0.0
        @test Distributions.ccdf(d, 0.0) == 1.0
        @test Distributions.cdf(d, 100.0) ≈ 1.0 atol=1e-6
        @test Distributions.ccdf(d, 100.0) ≈ 0.0 atol=1e-6
        @test Distributions.cdf(d, 5.0) + Distributions.ccdf(d, 5.0) ≈ 1.0 atol=1e-10

        # quantile
        @test Distributions.quantile(d, 0.5) isa Float64
        @test Distributions.quantile(d, 0.0) == 0.0
        @test Distributions.quantile(d, 1.0) == Inf

        # mean / var / std
        m = Distributions.mean(d)
        v = Distributions.var(d)
        s = Distributions.std(d)
        @test m isa Number
        @test v isa Number
        @test s isa Number
        @test v ≈ s^2 atol=1e-6
        @test m > 0
        @test v > 0

        # skewness / kurtosis
        @test Distributions.skewness(d) isa Number
        @test Distributions.kurtosis(d) isa Number

        # rand
        rng = Random.MersenneTwister(42)
        samples = Distributions.rand(rng, d, 100)
        @test length(samples) == 100
        @test all(x -> x >= 0, samples)
        @test mean(samples) ≈ m atol=1.0

        # insupport / minimum / maximum
        @test Distributions.insupport(d, 0.0) == true
        @test Distributions.insupport(d, -1.0) == false
        @test Distributions.minimum(d) == 0.0
        @test Distributions.maximum(d) == Inf

        # params
        λ, k, T = Distributions.params(d)
        @test λ == 3.0
        @test k == 2.0
        @test T == 1.0
    end

    # ============================================================
    # 8. Robustez asintótica vs ArbNumerics (method=:asymp / :arb)
    # ============================================================
    @testset "Robustez asintótica vs ArbNumerics" begin
        # La función puig_pdf_arb es interna del módulo y sirve como referencia
        # de alta precisión (ArbNumerics 128-bit) para validar :asymp

        # 8.1 Casos extremos de λ
        for λv in [10.0, 50.0, 100.0, 250.0, 500.0, 1000.0]
            kv = 2.4
            Tv = 0.01
            for frac in [0.7, 0.9, 1.0, 1.1, 1.3]
                xv = frac * λv
                v_asymp = Puig_pdf(xv, λv, kv, Tv; method=:asymp)
                v_arb   = Puig_pdf(xv, λv, kv, Tv; method=:arb)
                if v_arb > 1e-15
                    err = abs(v_asymp - v_arb) / v_arb
                    @test err < 1e-12
                else
                    @test abs(v_asymp - v_arb) < 1e-15
                end
            end
        end

        # 8.2 Casos extremos de k
        for kv in [1.0001, 1.2, 2.0, 3.0, 5.0, 10.0]
            λv = 50.0
            Tv = 0.05
            for frac in [0.8, 1.0, 1.2]
                xv = frac * λv
                v_asymp = Puig_pdf(xv, λv, kv, Tv; method=:asymp)
                v_arb   = Puig_pdf(xv, λv, kv, Tv; method=:arb)
                if v_arb > 1e-15
                    err = abs(v_asymp - v_arb) / v_arb
                    @test err < 1e-12
                end
            end
        end

        # 8.3 Colas extremas (lejos de la media)
        λv = 100.0; kv = 2.5; Tv = 0.05
        @test Puig_pdf(0.0, λv, kv, Tv; method=:asymp) == 0.0
        @test_throws ArgumentError Puig_pdf(-5.0, λv, kv, Tv; method=:asymp)
        @test Puig_pdf(1000.0, λv, kv, Tv; method=:asymp) == 0.0
        @test isfinite(Puig_pdf(100.0, λv, kv, Tv; method=:asymp))
        @test Puig_pdf(100.0, λv, kv, Tv; method=:asymp) > 0.0

        # 8.4 Coherencia: asymp y arb dan resultados idénticos en zona intermedia
        d_test = PuigDistribution(50.0, 2.4, 0.05)
        @test Puig_pdf(50.0, d_test; method=:asymp) ≈ Puig_pdf(50.0, d_test; method=:arb) atol=1e-12
        @test Puig_pdf(100.0, d_test; method=:asymp) ≈ Puig_pdf(100.0, d_test; method=:arb) atol=1e-12

        # 8.5 Distributions.pdf usa :asymp (verificación de que no cambia)
        @test Distributions.pdf(d_test, 50.0) ≈ Puig_pdf(50.0, d_test; method=:asymp)
    end

    # ============================================================
    # 9. Momentos Genéricos de Orden n
    # ============================================================
    @testset "Momentos Genéricos de Orden n" begin
        using SpecialFunctions

        # 9.1 Retrocompatibilidad: Puig_moments(params) == Puig_moments(params, 4)
        d9 = PuigDistribution(2.5, 3.2, 1.4)
        tuple4 = Puig_moments(d9)
        vec4   = Puig_moments(d9, 4)
        @test length(vec4) == 4
        @test vec4[1] ≈ tuple4[1]
        @test vec4[2] ≈ tuple4[2]
        @test vec4[3] ≈ tuple4[3]
        @test vec4[4] ≈ tuple4[4]

        # 9.2 Puig_moments con n posicional
        vec8 = Puig_moments(d9, 8)
        @test length(vec8) == 8
        @test vec8[1] ≈ vec4[1]
        @test vec8[2] ≈ vec4[2]
        @test vec8[3] ≈ vec4[3]
        @test vec8[4] ≈ vec4[4]

        # 9.3 Puig_moments con n positional (retorna Vector{Float64})
        vec8_dup = Puig_moments(d9, 8)
        @test vec8_dup == vec8

        # 9.4 Los momentos pares son positivos y crecientes para λ > 0
        @test all(m -> m > 0, vec8[2:2:end])

        # 9.5 Verificación: μ₅ debe ser coherente con μ₁, μ₂, μ₃, μ₄
        # Usamos la relación de recurrencia manualmente para verificar μ₅:
        # μ₅ = ((2*5 - 4 + k)/T + λ²) * μ₃ - (5-2)*(5+k-4)/T² * μ₁
        λf, kf, Tf = 2.5, 3.2, 1.4
        μ1, μ2, μ3, μ4 = vec4
        μ5_rec = ((10.0 - 4.0 + kf) / Tf + λf^2) * μ3 - 3.0 * (1.0 + kf) / Tf^2 * μ1
        @test vec8[5] ≈ μ5_rec atol=1e-10

        # 9.6 Verificación: μ₆ recurrencia manual
        μ6_rec = ((12.0 - 4.0 + kf) / Tf + λf^2) * μ4 - 4.0 * (2.0 + kf) / Tf^2 * μ2
        @test vec8[6] ≈ μ6_rec atol=1e-10

        # 9.7 Límite Maxwell-Boltzmann (λ = 0): μₙ = (2B)^{n/2} Γ((k+n)/2) / Γ(k/2)
        B_mb = 1.0
        k_mb = 2.0
        mb_dist = MB(k_mb, B_mb)
        vec_mb = Puig_moments(mb_dist, 8)
        for nn in 1:8
            μ_analitico = (2.0 * B_mb)^(nn / 2) * gamma((k_mb + nn) / 2) / gamma(k_mb / 2)
            @test vec_mb[nn] ≈ μ_analitico atol=1e-10
        end

        # 9.8 Puig_moments funciona con n=1, n=2, n=3
        vec1 = Puig_moments(d9, 1)
        @test length(vec1) == 1
        @test vec1[1] ≈ vec4[1]

        vec2 = Puig_moments(d9, 2)
        @test length(vec2) == 2
        @test vec2[2] ≈ vec4[2]

        vec3 = Puig_moments(d9, 3)
        @test length(vec3) == 3
        @test vec3[3] ≈ vec4[3]

        # 9.9 Puig_moments para MB con n positional
        vec_mb_dup = Puig_moments(mb_dist, 6)
        @test length(vec_mb_dup) == 6
        @test vec_mb_dup[1] ≈ vec_mb[1]

        # 9.10 Validación de errores
        @test_throws ArgumentError Puig_moments(d9, 0)

        # 9.11 Puig_stats usa momentos nativamente (coherencia)
        st = Puig_stats(d9)
        @test st.mean ≈ vec4[1]
        @test st.var  ≈ vec4[2] - vec4[1]^2
    end

    # ============================================================
    # 10. Marcum-Q Asintótico y Supervivencia en Colas Extremas
    # ============================================================
    @testset "Marcum-Q Asintótico y Supervivencia en Colas Extremas" begin
        # 10.1 marcumq_asymp produce resultados razonables
        @test marcumq_asymp(10.0, 15.0, 1.0) > 0.0
        @test marcumq_asymp(10.0, 15.0, 1.0) < 1.0

        # 10.2 En colas extremas, marcutq no produce falsos positivos ~1e-14
        λ_t, k_t, T_t = 10.0, 3.2, 1.0
        sT = sqrt(T_t)
        a_t = λ_t * sT
        b_extreme = 25.0 * sT  # b >> a
        mq_val = marcumq(a_t, b_extreme, k_t / 2)
        @test mq_val < 1e-30   # Debe ser extremadamente pequeño

        # 10.3 Supervivencia en cola superior lejana no produce NaN
        x_extreme = [λ_t + 10.0, λ_t + 15.0, λ_t + 20.0]
        S = Puig_surviving(x_extreme, λ_t, k_t, T_t)
        @test all(isfinite, S)
        @test all(s -> s >= 0.0, S)
        # Monotonía decreciente
        @test issorted(S, rev=true)

        # 10.4 Continuidad: valores de supervivencia en la frontera del régimen asintótico
        # Para b - a justo en el límite (~4.0), la supervivencia debe ser continua
        x_boundary = [λ_t + 4.0 / sT]
        S_bnd = Puig_surviving(x_boundary, λ_t, k_t, T_t)
        @test isfinite(S_bnd[1])
        @test S_bnd[1] > 0.0
        @test S_bnd[1] < 1.0

        # 10.5 marcumq escalar vs vectorizado: deben ser consistentes
        x_vec_test = [λ_t + 5.0, λ_t + 8.0, λ_t + 12.0]
        mq_scalar = [marcumq(a_t, x_vec_test[i] * sT, k_t / 2) for i in 1:3]
        mq_vector = marcumq(a_t, x_vec_test .* sT, k_t / 2)
        @test mq_scalar ≈ mq_vector atol=1e-14

        # 10.6 Caso extremo b >> a: exponential decay monótono sin oscilaciones
        x_remote = [λ_t + 10.0, λ_t + 12.0, λ_t + 14.0, λ_t + 16.0, λ_t + 18.0]
        S_remote = Puig_surviving(x_remote, λ_t, k_t, T_t)
        @test all(isfinite, S_remote)
        @test issorted(S_remote, rev=true)

        # 10.7 Sin argumentos espurios: para ab < 30, fallback a pnchisq
        a_small, b_small = 1.0, 3.0  # ab = 3 < 30
        mq_small = marcumq(a_small, b_small, 1.0)
        @test isfinite(mq_small)
        @test mq_small > 0.0
        @test mq_small < 1.0

        # 10.8 marcumq_asymp con m=2 produce valores mayores que con m=1 (b > a)
        @test marcumq_asymp(10.0, 20.0, 2.0) > marcumq_asymp(10.0, 20.0, 1.0)
    end

end
