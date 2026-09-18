function out = Puig_cumulative(x, lam, k, T)
    % PUIG_CUMULATIVE  Acumulada F(x) = P(X ≤ x) = 1 - S(x).
    % Acepta (x, lam, k, T) o (x, dist).
    arguments
        x double {mustBeNonnegative}
        lam double
        k double = []
        T double = []
    end
    if isa(lam, 'puigdist.PuigDistribution') || isa(lam, 'puigdist.MB')
        out = 1.0 - puigdist.Puig_surviving(x, lam);
    else
        out = 1.0 - puigdist.Puig_surviving(x, lam, k, T);
    end
end