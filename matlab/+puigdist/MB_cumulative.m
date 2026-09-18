function out = MB_cumulative(x, k, B)
    % MB_CUMULATIVE  Acumulada F(x) = 1 - S(x).
    arguments
        x double {mustBeNonnegative}
        k (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        B (1,1) double {mustBePositive}
    end
    out = 1.0 - puigdist.MB_surviving(x, k, B);
end