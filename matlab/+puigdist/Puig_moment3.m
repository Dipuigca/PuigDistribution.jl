function mu3 = Puig_moment3(lam, k, T)
    % PUIG_MOMENT3  Tercer momento raw μ₃ = (3/T)·√(π/(2T))·L_{3/2}^{(k/2-1)}(-λ²T/2).
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
    end
    ls = lam * sqrt(T);
    L = puigdist.laguerre_real(1.5, k / 2.0 - 1.0, -ls^2 / 2.0);
    mu3 = (3.0 / T) * sqrt(pi / (2.0 * T)) * L;
end