function ci = Puig_ci98(lam, k, T)
    % PUIG_CI98  Intervalo central al 98%: (P1, P99).
    q = puigdist.Puig_quantile([0.01 0.99], lam, k, T);
    ci = [q(1) q(2)];
end