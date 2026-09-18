function ci = Puig_ci996(lam, k, T)
    % PUIG_CI996  Intervalo central al 99.8%: (P0.1, P99.9).
    q = puigdist.Puig_quantile([0.001 0.999], lam, k, T);
    ci = [q(1) q(2)];
end