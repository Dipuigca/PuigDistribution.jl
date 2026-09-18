function out = Puig_pdf(x, lam, k, T, method)
    % PUIG_PDF  Densidad de probabilidad de la distribución de Puig.
    %
    %   f_P(x) = T·x^{k/2}/λ^{k/2-1}·exp(-T/2·(x-λ)²)·e^{-λ²T/2}·I_{k/2-1}(xλT)
    %
    % method: 'asymp' (por defecto) usa la expansión asintótica de 5 términos
    % cuando xλT ≥ 200 y besseli escalada en otro caso; 'arb' es una
    % referencia de precisión doble equivalente (besseli exacta para z < 200).
    % Para λ = 0 reduce a MB_pdf(x, k, 1/T).
    arguments
        x double {mustBeNonnegative}
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
        method string = "asymp"
    end

    if lam == 0
        out = puigdist.MB_pdf(x, k, 1.0 / T);
        return;
    end

    scalar = isscalar(x);
    x = x(:);
    nu = k / 2.0 - 1.0;
    half_T = 0.5 * T;
    pow_k_half = k / 2.0;
    lam_denom = lam^(pow_k_half - 1.0);

    out = zeros(size(x));
    for i = 1:numel(x)
        xi = x(i);
        if xi > 0
            arg_exp = -half_T * (xi - lam)^2;
            if arg_exp >= -740.0
                expo = exp(arg_exp);
                if expo > 0.0
                    z = xi * lam * T;
                    if z >= 200.0
                        bix = puigdist.besselix_asymp(nu, z);
                    else
                        bix = besseli(nu, z, 1);
                    end
                    factor = T * (xi^pow_k_half) / lam_denom;
                    val = factor * expo * bix;
                    if isfinite(val)
                        out(i) = val;
                    end
                end
            end
        end
    end
    if scalar
        out = out(1);
    end
end