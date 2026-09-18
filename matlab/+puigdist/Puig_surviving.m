function out = Puig_surviving(x, lam, k, T)
    % PUIG_SURVIVING  Supervivencia S(x) = P(X ≥ x) = Q_{k/2}(λ√T, x√T).
    % Acepta (x, lam, k, T) o (x, dist) con dist un PuigDistribution.
    arguments
        x double {mustBeNonnegative}
        lam double
        k double = []
        T double = []
    end

    if isa(lam, 'puigdist.PuigDistribution')
        d = lam; lam = d.lam; k = d.k; T = d.T;
    elseif isa(lam, 'puigdist.MB')
        d = lam; lam = 0.0; k = d.k; T = 1.0 / d.B;
    end

    scalar = isscalar(x);
    x = x(:);

    if lam == 0
        out = puigdist.MB_surviving(x, k, 1.0 / T);
    else
        order = k / 2.0;
        sT = sqrt(T);
        a = lam * sT;
        b = x * sT;
        out = puigdist.marcumq(a, b, order);
        % Réplica del ajuste de robustez de Julia
        for i = 1:numel(out)
            if isnan(out(i))
                error('marcumq devolvió NaN en x=%g, a=%g, m=%g', x(i), a, order);
            end
            if out(i) <= 0.0 && x(i) < lam
                out(i) = 1.0;
            end
        end
    end
    if scalar
        out = out(1);
    end
end