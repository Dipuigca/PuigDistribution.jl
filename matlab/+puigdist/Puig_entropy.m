function h = Puig_entropy(lam, k, T, N)
    % PUIG_ENTROPY  Entropía diferencial H(X) = -∫ f(x)·log f(x) dx (regla del trapecio).
    arguments
        lam (1,1) double {mustBeNonnegative}
        k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        T   (1,1) double {mustBePositive}
        N   (1,1) double = 400
    end
    ci = puigdist.Puig_ci996(lam, k, T);
    init = max(ci(1) / 10.0, 1e-6);
    fin = ci(2);
    x = linspace(init, fin, N);
    f = puigdist.Puig_pdf(x, lam, k, T);
    lf = puigdist.Puig_logpdf(x, lam, k, T);
    F = -f .* lf;
    hstep = (fin - init) / (N - 1);
    FH = F * hstep;
    FH(1) = FH(1) / 2.0;
    FH(end) = FH(end) / 2.0;
    h = sum(FH);
end