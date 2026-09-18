function out = MB_surviving(x, k, B)
    % MB_SURVIVING  Supervivencia S(x) = Γ(k/2, x²/(2B)) / Γ(k/2) (gamma regularizada superior).
    arguments
        x double {mustBeNonnegative}
        k (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        B (1,1) double {mustBePositive}
    end
    P = (x.^2) ./ (2.0 * B);
    out = gammainc(P, k / 2.0, 'upper');
end