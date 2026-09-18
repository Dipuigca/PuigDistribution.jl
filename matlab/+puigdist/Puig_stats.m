function st = Puig_stats(lam, k, T)
    % PUIG_STATS  Estadísticos resumen: media, varianza, σ, skewness, kurtosis.
    mu1 = puigdist.Puig_mean(lam, k, T);
    sig2 = puigdist.Puig_var(lam, k, T);
    sig = sqrt(sig2);
    mu3 = puigdist.Puig_moment3(lam, k, T);
    mu2 = k / T + lam^2;
    mu4 = mu2^2 + 2.0 * k / T^2 + 4.0 * lam^2 / T;
    sk = (mu3 - 3.0 * mu1 * mu2 + 2.0 * mu1^3) / sig^3;
    kur = (mu4 - 4.0 * mu1 * mu3 + 6.0 * mu1^2 * mu2 - 3.0 * mu1^4) / sig^4;
    st = struct('mean', mu1, 'var', sig2, 'sig', sig, ...
        'skewness', sk, 'kurtosis', kur);
end