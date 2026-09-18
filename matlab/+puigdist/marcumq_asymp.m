function v = marcumq_asymp(a, b, m)
    % MARCUMQ_ASYMP  Q_m(a, b) por expansión asintótica de 2 términos (Cantrell 1986).
    % Requiere b - a ≥ 4.0 y ab ≥ 30.0 (cola superior lejana). Evita
    % cancelaciones catastróficas y falsos suelos de ruido (~1e-50).
    % Port de marcumq_asymp de src/puig_pdf.jl (DistributionsPuig.jl).
    arguments
        a (1,1) double
        b (1,1) double
        m (1,1) double = 1.0
    end
    if (b - a) < 4.0
        error('marcumq_asymp requiere b - a >= 4.0');
    end
    if a * b < 30.0
        error('marcumq_asymp requiere ab >= 30.0');
    end

    diff = b - a;
    arg_exp = -0.5 * diff * diff;
    if arg_exp < -740.0
        v = 0.0;
        return;
    end

    rho = b / a;
    prefactor = rho^(m - 0.5) .* exp(arg_exp) ./ (sqrt(2.0 * pi) .* diff);
    num1 = b + a + (4.0 .* (m - 0.5).^2 - 1.0) ./ (4.0 .* diff);
    denom1 = 2.0 .* a .* b .* diff;
    corr = 1.0 - num1 ./ denom1;
    v = prefactor .* corr;
    v = max(v, 0.0);
    if ~isfinite(v)
        v = 0.0;
    end
end