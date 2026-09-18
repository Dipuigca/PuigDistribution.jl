function res = Puig_fit(Data, opts)
    % PUIG_FIT  Ajusta datos a la distribución de Puig.
    %
    %   res = Puig_fit(Data)
    %   res = Puig_fit(Data, 'method', ..., 'vars', ..., 'time_col', ..., 'k_fixed', ...)
    %
    % Entradas:
    %   table   → columnas numéricas; el indicador se calcula como √(Σ vars²) o
    %             se toma de time_col. opts.vars: cell/string; opts.time_col: string.
    %   m×d double → filas = observaciones, columnas = variables.
    %   vector 1×n / n×1 → serie de normas observadas.
    %
    % Métodos (matriz): 'dynamic' (k por correlación), 'fixed' (k = d), 'numerical' (2D).
    % Métodos (vector): 'dynamic'/'numerical' (3D Nelder-Mead) o 'fixed' con k_fixed.
    %
    % Devuelve struct con campos: params (PuigDistribution), indicator, stats (struct), method.
    arguments
        Data
        opts.vars          = []
        opts.time_col      = []
        opts.method string = "dynamic"
        opts.k_fixed       = []
    end

    if istable(Data)
        res = fit_table(Data, opts);
    elseif ismatrix(Data) && isnumeric(Data)
        if size(Data, 2) > 1
            res = fit_matrix(Data, opts);
        else
            res = fit_vector(Data(:), opts);
        end
    else
        error('Entrada no soportada');
    end
end

% -----------------------------------------------------------------------
% Ajuste desde table
% -----------------------------------------------------------------------
function res = fit_table(T, opts)
    var_names = opts.vars;
    time_col = opts.time_col;
    if isempty(var_names)
        if isempty(time_col)
            exclude = {};
        else
            exclude = {time_col};
        end
        var_names = T.Properties.VariableNames;
        keep = cellfun(@(c) isnumeric(T.(c)) && ~ismember(c, exclude), var_names);
        var_names = var_names(keep);
    end
    if isempty(var_names)
        error('Se requiere al menos una columna numérica para el ajuste');
    end
    mat = T{:, var_names};
    if isempty(time_col)
        times_vec = [];
    else
        times_vec = T{:, time_col};
    end
    res = fit_matrix(mat, struct('method', opts.method, 'times', times_vec));
end

% -----------------------------------------------------------------------
% Ajuste desde matrix (m×d)
% -----------------------------------------------------------------------
function res = fit_matrix(Data, opts)
    n_obs = size(Data, 1);
    d     = size(Data, 2);
    if n_obs < 10
        error('Se requieren al menos 10 observaciones para el ajuste (recibidas: %d)', n_obs);
    end

    if isfield(opts, 'times') && ~isempty(opts.times)
        t_raw = opts.times(:);
    else
        t_raw = sqrt(sum(Data.^2, 2));
    end
    if numel(t_raw) ~= n_obs
        error('times debe tener tantas observaciones como filas de Data');
    end

    valid = isfinite(t_raw) & (t_raw >= 0);
    t = t_raw(valid);
    if numel(t) < 10
        error('Se requieren al menos 10 observaciones válidas para el ajuste');
    end

    col_means = mean(Data(valid, :), 1);
    lam = sqrt(sum(col_means.^2));
    k_corr = puigdist.effective_dimension(Data(valid, :));

    min_t = min(t);
    max_t = max(t);
    x_opt = linspace(min_t, max_t, 150);
    S_emp_opt = puigdist.empirical_survival(t, x_opt);

    sdt = std(t);
    if isfinite(sdt) && sdt > 0
        B0 = sdt;
    else
        B0 = 1.0;
    end

    method = string(opts.method);
    k_est = k_corr;
    B_est = B0;

    if method == "dynamic" || method == "fixed"
        if method == "dynamic"
            k_target = k_corr;
        else
            k_target = double(d);
        end
        fun = @(B) err_l1(B, x_opt, S_emp_opt, lam, k_target);
        B_est = fminbnd(fun, 0.01, 5.0 * B0, optimset('TolX', 1e-12));
        k_est = k_target;
    elseif method == "numerical"
        fun = @(par) err_l2_matrix(par, x_opt, S_emp_opt, lam);
        init_par = [min(max(k_corr, 1.0), max(3.0, d)), B0];
        lb = [1.0, 0.001];
        ub = [max(3.0, double(d)), 10.0 * B0];
        opts_fms = optimset('MaxIter', 1000, 'TolX', 1e-8, 'TolFun', 1e-8);
        res_fms = fminsearch(fun, init_par, opts_fms);
        k_est = max(lb(1), min(ub(1), res_fms(1)));
        B_est = max(lb(2), min(ub(2), res_fms(2)));
    else
        error('Método desconocido: %s. Opciones: dynamic, fixed, numerical', method);
    end

    T_est = 1.0 / B_est^2;
    x_eval = linspace(min_t, max_t, 500);
    S_emp_eval = puigdist.empirical_survival(t, x_eval);
    S_mod_eval = puigdist.Puig_surviving(x_eval, lam, k_est, T_est);
    fit_metrics = puigdist.calc_all_metrics(x_eval, S_emp_eval, S_mod_eval);

    params = puigdist.PuigDistribution(lam, k_est, T_est);
    res = struct('params', params, 'indicator', t, 'stats', fit_metrics, 'method', method);
end

% -----------------------------------------------------------------------
% Ajuste desde vector 1D
% -----------------------------------------------------------------------
function res = fit_vector(Data, opts)
    valid = isfinite(Data) & (Data >= 0);
    t = Data(valid);
    n_obs = numel(t);
    if n_obs < 10
        error('Se requieren al menos 10 observaciones para el ajuste (recibidas: %d)', n_obs);
    end

    min_t = min(t);
    max_t = max(t);
    x_opt = linspace(min_t, max_t, 60);
    S_emp_opt = puigdist.empirical_survival(t, x_opt);

    m1 = mean(t);
    m2 = mean(t.^2);
    s = std(t);
    if isfinite(s) && s > 0
        B0 = s;
    else
        B0 = 1.0;
    end
    lam0 = sqrt(max(0.0, m1^2 - s^2));
    k0 = min(max((m2 - lam0^2) / B0^2, 1.0), 10.0);

    method = string(opts.method);
    k_fixed = opts.k_fixed;
    is_fixed_k = (method == "fixed") || ~isempty(k_fixed);

    if is_fixed_k
        if isempty(k_fixed)
            k_target = double(size(Data, 2));
        else
            k_target = double(k_fixed);
        end
        if k_target < 1.0
            error('k_fixed debe ser >= 1.0 (recibido: %g)', k_target);
        end
        init_p = [log(max(1e-3, lam0)), log(max(1e-3, B0))];
        fun = @(p) err_2d_fixed(p, x_opt, S_emp_opt, k_target);
        opts_fms = optimset('MaxIter', 500, 'TolX', 1e-8, 'TolFun', 1e-8);
        res_fms = fminsearch(fun, init_p, opts_fms);
        lam_est = exp(res_fms(1));
        B_est = exp(res_fms(2));
        k_est = k_target;
    elseif method == "dynamic" || method == "numerical"
        init_p = [log(max(1e-3, lam0)), log(max(1e-3, k0 - 1.0)), log(max(1e-3, B0))];
        fun = @(p) err_3d(p, x_opt, S_emp_opt);
        opts_fms = optimset('MaxIter', 600, 'TolX', 1e-8, 'TolFun', 1e-8);
        res_fms = fminsearch(fun, init_p, opts_fms);
        lam_est = exp(res_fms(1));
        k_est = 1.0 + exp(res_fms(2));
        B_est = exp(res_fms(3));
        if k_est > 10.0
            k_est = 10.0;
        end
    else
        error('Método desconocido: %s. Opciones: dynamic, numerical, fixed', method);
    end

    T_est = 1.0 / B_est^2;
    x_eval = linspace(min_t, max_t, 500);
    S_emp_eval = puigdist.empirical_survival(t, x_eval);
    S_mod_eval = puigdist.Puig_surviving(x_eval, lam_est, k_est, T_est);
    fit_metrics = puigdist.calc_all_metrics(x_eval, S_emp_eval, S_mod_eval);

    params = puigdist.PuigDistribution(lam_est, k_est, T_est);
    res = struct('params', params, 'indicator', t, 'stats', fit_metrics, 'method', method);
end

% -----------------------------------------------------------------------
% Funciones objetivo (subfunctions)
% -----------------------------------------------------------------------
function e = err_l1(B, x_opt, S_emp, lam, k_target)
    if B <= 0
        e = 1e12; return;
    end
    T_val = 1.0 / B^2;
    try
        S_mod = puigdist.Puig_surviving(x_opt, lam, k_target, T_val);
        e = sum(abs(S_emp - S_mod(:)));
    catch
        e = 1e12;
    end
end

function e = err_l2_matrix(par, x_opt, S_emp, lam)
    k_val = par(1); B_val = par(2);
    if k_val < 1.0 || B_val <= 0.0
        e = 1e12; return;
    end
    T_val = 1.0 / B_val^2;
    try
        S_mod = puigdist.Puig_surviving(x_opt, lam, k_val, T_val);
        if ~all(isfinite(S_mod))
            e = 1e12; return;
        end
        e = sum((S_emp - S_mod(:)).^2);
    catch
        e = 1e12;
    end
end

function e = err_2d_fixed(p, x_opt, S_emp, k_target)
    lam_val = exp(p(1));
    B_val = exp(p(2));
    T_val = 1.0 / B_val^2;
    try
        S_mod = puigdist.Puig_surviving(x_opt, lam_val, k_target, T_val);
        if ~all(isfinite(S_mod))
            e = 1e12; return;
        end
        e = sum((S_emp - S_mod(:)).^2);
    catch
        e = 1e12;
    end
end

function e = err_3d(p, x_opt, S_emp)
    lam_val = exp(p(1));
    k_val = 1.0 + exp(p(2));
    B_val = exp(p(3));
    if k_val > 10.0
        e = 1e12; return;
    end
    T_val = 1.0 / B_val^2;
    try
        S_mod = puigdist.Puig_surviving(x_opt, lam_val, k_val, T_val);
        if ~all(isfinite(S_mod))
            e = 1e12; return;
        end
        e = sum((S_emp - S_mod(:)).^2);
    catch
        e = 1e12;
    end
end