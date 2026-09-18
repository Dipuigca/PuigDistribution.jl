function kur = Puig_kurtosis(lam, k, T)
    % PUIG_KURTOSIS  Kurtosis (4º momento estandarizado no restado) γ₂ = E[(X-μ)⁴]/σ⁴.
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end
    mu1 = puigdist.Puig_mean(lam, k, T);
    mu2 = k / T + lam^2;
    mu4 = mu2^2 + 2.0 * k / T^2 + 4.0 * lam^2 / T;
    mu3 = puigdist.Puig_moment3(lam, k, T);
    sigma4 = puigdist.Puig_var(lam, k, T)^2;
    kur = (mu4 - 4.0 * mu1 * mu3 + 6.0 * mu1^2 * mu2 - 3.0 * mu1^4) / sigma4;
end