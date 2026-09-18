function q = Puig_quantile(p, lam, k, T, tol, maxit)
    % PUIG_QUANTILE  Cuantil Q(p) tal que F(Q(p)) = p, por bisección.
    arguments
        p double
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
        tol (1,1) double = 1e-10
        maxit (1,1) double = 200
    end
    if any(p < 0) || any(p > 1)
        error('p debe estar en [0, 1]');
    end
    scalar = isscalar(p);
    p = p(:);
    q = zeros(size(p));
    for i = 1:numel(p)
        pi_val = p(i);
        if pi_val == 0
            q(i) = 0.0;
        elseif pi_val == 1
            q(i) = Inf;
        else
            lo = 0.0;
            hi = max(lam, 1.0) * 2.0;
            while puigdist.Puig_cumulative(hi, lam, k, T) < pi_val
                hi = hi * 2.0;
            end
            for iter = 1:maxit
                mid = 0.5 * (lo + hi);
                if puigdist.Puig_cumulative(mid, lam, k, T) < pi_val
                    lo = mid;
                else
                    hi = mid;
                end
                if (hi - lo) < tol
                    break;
                end
            end
            q(i) = 0.5 * (lo + hi);
        end
    end
    if scalar
        q = q(1);
    end
end