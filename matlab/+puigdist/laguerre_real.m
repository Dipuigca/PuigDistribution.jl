function L = laguerre_real(n, alpha, x)
    % LAGUERRE_REAL  Función de Laguerre generalizada L_n^(α)(x) para n, α reales y x ≤ 0.
    %
    %   L_n^(α)(x) = Γ(n+α+1)/(Γ(n+1)Γ(α+1)) · e^{x}·M(α+1+n, α+1, -x)
    %
    % Con la transformación de Kummer M(-n, α+1, x) = e^{x}·M(α+1+n, α+1, -x).
    % En doble precisión:
    %   - y = -x < 20: serie de Kummer M(-n, b, -y) (serie directa estable).
    %   - y >= 20   : expansión asintótica (rama e^{-y}·M(a, b, y) ~ ...).
    % Port fiel a laguerre_real de src/puig_pdf.jl (con alta precisión en Julia).
    arguments
        n (1,1) double
        alpha (1,1) double
        x (1,1) double {mustBeLessThanOrEqual(x, 0)}
    end

    b = alpha + 1.0;
    y = -x;
    a = alpha + 1.0 + n;   % = b + n

    if y < 20.0
        f = kummer_series(-n, b, -y);
    else
        % e^{-y}·M(a, b, y) ≈ Γ(b)/Γ(a)·y^{a-b}·Σ_k (1-a)_k (b-a)_k / k! · y^{-k}
        f = asymp_series(a, b, y);
    end

    coef = gamma(n + alpha + 1.0) / (gamma(n + 1.0) * gamma(alpha + 1.0));
    L = coef * f;
end

function s = kummer_series(aa, bb, z)
    % M(a; b; z) por serie directa con criterio de parada relativo.
    term = 1.0;
    s = 1.0;
    maxk = 50000;
    for k = 1:maxk
        term = term * ((aa + k - 1.0) / (bb + k - 1.0)) * (z / k);
        s = s + term;
        if abs(term) < 1e-15 * max(abs(s), eps)
            break;
        end
    end
end

function s = asymp_series(a, b, y)
    % e^{-y}·M(a; b; y) ~ Γ(b)/Γ(a)·y^{a-b}·Σ_{k≥0} (1-a)_k (b-a)_k / k! · y^{-k}
    tot = 0.0;
    for k = 0:60
        c = pf(1.0 - a, k) * pf(b - a, k) / factorial(k);
        term = c * y^(-k);
        tot = tot + term;
        if k > 10 && abs(term) < 1e-16 * abs(tot)
            break;
        end
    end
    s = exp(gammaln(b) - gammaln(a)) * y^(a - b) * tot;
end

function r = pf(x, k)
    % Pochhammer (x)_k
    r = 1.0;
    for j = 0:(k - 1)
        r = r * (x + j);
    end
end