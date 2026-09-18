function st = calc_all_metrics(x, S_emp, S_mod)
    % CALC_ALL_METRICS  Métricas de bondad de ajuste sobre la supervivencia.
    % Devuelve struct con R2, MAE, RMSE, MaxAE, IAE (trapecio).
    S_e = S_emp(:); S_m = S_mod(:);
    denom = sum((S_e - mean(S_e)).^2);
    if denom == 0
        R2 = 1.0;
    else
        R2 = 1.0 - sum((S_e - S_m).^2) / denom;
    end
    MAE = mean(abs(S_e - S_m));
    RMSE = sqrt(mean((S_e - S_m).^2));
    MaxAE = max(abs(S_e - S_m));
    xv = x(:);
    diffs = abs(S_e - S_m);
    dx = diff(xv);
    IAE = sum(dx .* (diffs(1:end-1) + diffs(2:end)) / 2.0);
    st = struct('R2', R2, 'MAE', MAE, 'RMSE', RMSE, 'MaxAE', MaxAE, 'IAE', IAE);
end