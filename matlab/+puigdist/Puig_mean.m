function mu1 = Puig_mean(lam, k, T)
    % PUIG_MEAN  Media teórica E[X] = √(π/(2T))·L_{1/2}^{(k/2-1)}(-λ²T/2).
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end
    ls = lam * sqrt(T);
    mu1 = sqrt(pi / (2.0 * T)) * puigdist.laguerre_real(0.5, k / 2.0 - 1.0, -ls^2 / 2.0);
end