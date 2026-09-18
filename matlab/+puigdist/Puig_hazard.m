function out = Puig_hazard(x, lam, k, T)
    % PUIG_HAZARD  Función de riesgo h(x) = f(x) / S(x).
    f = puigdist.Puig_pdf(x, lam, k, T);
    S = puigdist.Puig_surviving(x, lam, k, T);
    out = f ./ S;
end