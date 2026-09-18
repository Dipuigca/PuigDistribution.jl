function sk = Puig_skewness(lam, k, T)
    % PUIG_SKEWNESS  Asimetría estandarizada γ₁ = E[(X-μ)³]/σ³.
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end
    mu1 = puigdist.Puig_mean(lam, k, T);
    mu2 = puigdist.Puig_var(lam, k, T) + mu1^2;
    mu3 = puigdist.Puig_moment3(lam, k, T);
    sigma3 = puigdist.Puig_std(lam, k, T)^3;
    sk = (mu3 - 3.0 * mu1 * mu2 + 2.0 * mu1^3) / sigma3;
end