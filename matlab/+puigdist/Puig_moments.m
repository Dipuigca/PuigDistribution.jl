function mu = Puig_moments(lam, k, T, n)
    % PUIG_MOMENTS  Los n primeros momentos raw μ₁,…,μₙ (recurrencia de tres términos O(n)).
    %
    % Semillas: μ₁ (Laguerre real), μ₂ = k/T + λ², μ₃ (Laguerre), μ₄ = μ₂² + 2k/T² + 4λ²/T.
    % Recurrencia (n ≥ 5):
    %   μₙ = ((2n-4+k)/T + λ²)·μₙ₋₂ - (n-2)(n+k-4)/T²·μₙ₋₄
    %
    % Port de Puig_moments_raw / puig_moments de DistributionsPuig.jl.
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
        n   (1,1) double {mustBeInteger, mustBeGreaterThan(n, 0)}
    end
    n = int32(n);
    mu = zeros(1, n);
    T2 = T * T;
    mu(1) = puigdist.Puig_mean(lam, k, T);
    if n >= 2
        mu(2) = k / T + lam^2;
    end
    if n >= 3
        mu(3) = puigdist.Puig_moment3(lam, k, T);
    end
    if n >= 4
        mu(4) = mu(2)^2 + 2.0 * k / T2 + 4.0 * lam^2 / T;
    end
    for r = 5:n
        rr = double(r);
        mu(r) = ((2.0 * rr - 4.0 + k) / T + lam^2) * mu(r - 2) ...
            - (rr - 2.0) * (rr + k - 4.0) / T2 * mu(r - 4);
    end
end