function figh = Puig_plot(dist, options)
    % PUIG_PLOT  Representaciones gráficas de la distribución de Puig.
    %
    %   Puig_plot(dist)                        → figura 2×2 con pdf, cdf, surv, datos.
    %   Puig_plot(dist, 'pdf')                 → solo densidad.
    %   Puig_plot(dist, ["pdf" "surv"])        → paneles seleccionados.
    %   Puig_plot(MB(k, B))                    → objeto MB también aceptado.
    %
    % Devuelve el handle de la figura.
    arguments
        dist
        options string = "all"
    end

    if isa(dist, 'puigdist.MB')
        lam = 0.0;
        k = dist.k;
        T = 1.0 / dist.B;
    elseif isa(dist, 'puigdist.PuigDistribution')
        lam = dist.lam;
        k = dist.k;
        T = dist.T;
    else
        error('dist debe ser un objeto PuigDistribution o MB');
    end

    if isequal(options, "all")
        options = ["pdf" "cdf" "surv" "data"];
    end
    known = ["pdf" "cdf" "surv" "data"];
    bad = setdiff(options, known);
    if ~isempty(bad)
        error('Variable de opción desconocida: %s', strjoin(bad, ', '));
    end

    ci98 = puigdist.Puig_quantile([0.01 0.99], lam, k, T);
    st = puigdist.Puig_stats(lam, k, T);
    qs = puigdist.Puig_quantile([0.25 0.5 0.75], lam, k, T);

    txt = sprintf(['\\lambda = %.4f   k = %.4f   T = %.4f\n\n' ...
        'media    = %.4f\nvarianza = %.4f\n\\sigma    = %.4f\n' ...
        'skewness = %.4f\nkurtosis = %.4f\n\nP1  = %.4f   P99 = %.4f\n\n' ...
        'Q1 (25%%) = %.4f\nQ2 (50%%) = %.4f\nQ3 (75%%) = %.4f'], ...
        lam, k, T, st.mean, st.var, st.sig, st.skewness, st.kurtosis, ...
        ci98(1), ci98(2), qs(1), qs(2), qs(3));

    figh = figure('Name', 'Distribución de Puig', 'NumberTitle', 'off');
    n = numel(options);
    for i = 1:n
        subplot(2, 2, i);
        o = options(i);
        if o == "data"
            text(0.05, 0.5, txt, 'FontName', 'Consolas', 'FontSize', 9);
            axis off;
        else
            ci = puigdist.Puig_ci996(lam, k, T);
            init = max(ci(1) / 10.0, 1e-6);
            fin = ci(2);
            x = linspace(init, fin, 400);
            switch o
                case "pdf"
                    y = puigdist.Puig_pdf(x, lam, k, T);
                    plot(x, y, 'LineWidth', 1.2);
                    title('Densidad (pdf)'); xlabel('x'); ylabel('f(x)');
                case "cdf"
                    y = puigdist.Puig_cumulative(x, lam, k, T);
                    plot(x, y, 'LineWidth', 1.2);
                    title('Acumulada (cdf)'); xlabel('x'); ylabel('F(x)');
                case "surv"
                    y = puigdist.Puig_surviving(x, lam, k, T);
                    plot(x, y, 'LineWidth', 1.2);
                    title('Supervivencia'); xlabel('x'); ylabel('S(x)');
            end
            grid on;
        end
    end
end