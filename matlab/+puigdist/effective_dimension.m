function k = effective_dimension(M)
    % EFFECTIVE_DIMENSION  Dimensión efectiva k = d²/(d + 2·Σ_{i<j} R_ij²) en O(d²).
    % M: matriz n×d (filas = observaciones). Port de src/puig_fit.jl.
    if size(M, 2) <= 1
        k = 1.0;
        return;
    end
    d = size(M, 2);
    R = corrcoef(M);
    iu = triu(true(d), 1);
    s_off = sum(R(iu).^2);
    denom = d + 2.0 * s_off;
    if denom == 0.0
        k = 1.0;
    else
        k = d^2 / denom;
    end
end