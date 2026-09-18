function S = empirical_survival(times, grid)
    % EMPIRICAL_SURVIVAL  S_emp(x) = #{t_i ≥ x} / N evaluada sobre la rejilla grid.
    t = times(:);
    g = grid(:);
    n = numel(t);
    if n == 0
        S = zeros(size(g));
        return;
    end
    S = zeros(size(g));
    for i = 1:numel(g)
        S(i) = sum(t >= g(i)) / n;
    end
end