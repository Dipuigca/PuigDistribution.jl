function s = Puig_std(lam, k, T)
    % PUIG_STD  Desviación estándar σ = √Var(X).
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end
    s = sqrt(puigdist.Puig_var(lam, k, T));
end