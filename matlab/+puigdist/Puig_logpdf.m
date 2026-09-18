function out = Puig_logpdf(x, lam, k, T)
    % PUIG_LOGPDF  Log-densidad log f_P(x) numéricamente estable (besseli escalada).
    arguments
        x double {mustBeNonnegative}
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end

    scalar = isscalar(x);
    x = x(:);
    out = zeros(size(x));
    for i = 1:numel(x)
        xi = x(i);
        if xi == 0
            out(i) = -Inf;
        elseif lam == 0
            B = 1.0 / T;
            kh = k / 2.0;
            out(i) = (1.0 - kh) * log(2.0) - kh * log(B) ...
                - gammaln(kh) + (k - 1.0) * log(xi) - xi^2 / (2.0 * B);
        else
            nu = k / 2.0 - 1.0;
            z = xi * lam * T;
            if z >= 200.0
                logbix = log(puigdist.besselix_asymp(nu, z));
            else
                logbix = log(besseli(nu, z, 1));
            end
            out(i) = log(T) + (k / 2.0) * log(xi) ...
                - (k / 2.0 - 1.0) * log(lam) ...
                - T / 2.0 * (xi - lam)^2 + logbix;
        end
    end
    if scalar
        out = out(1);
    end
end