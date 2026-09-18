function x = Puig_rand(n, lam, k, T)
    % PUIG_RAND  Muestreo pseudoaleatorio exacto.
    %
    % Representación exacta: X = √(W/T) con W ~ NoncentralChisq(k, λ²T);
    % para λ = 0, X = √(W/T) con W ~ Gamma(k/2, 2/T).
    arguments
        n (1,1) {mustBePositive}
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end
    sz = [n 1];
    if lam == 0
        W = gamrnd(k / 2.0, 2.0 / T, sz);
    else
        W = ncx2rnd(k, lam^2 * T, sz);
    end
    x = sqrt(W / T);
end