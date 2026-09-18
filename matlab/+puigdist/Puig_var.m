function v = Puig_var(lam, k, T)
    % PUIG_VAR  Varianza teórica E[X²] - (E[X])² con E[X²] = k/T + λ².
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end
    mu2 = k / T + lam^2;
    v = mu2 - puigdist.Puig_mean(lam, k, T)^2;
end