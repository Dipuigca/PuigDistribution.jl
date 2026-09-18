function out = MB_pdf(x, k, B)
    % MB_PDF  Densidad del caso límite Maxwell-Boltzmann / chi generalizada.
    arguments
        x double {mustBeNonnegative}
        k (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
        B (1,1) double {mustBePositive}
    end
    kh = k / 2.0;
    E = exp(-(x.^2) / (2.0 * B));
    G = gamma(kh);
    C = 2.0^(1.0 - kh) * B^(-kh) / G;
    out = (x.^(k - 1.0)) .* E .* C;
end