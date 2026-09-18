classdef PuigDistribution
    % Distribución de Puig: chi no central generalizada de dimensión real k.
    %
    %   f_P(x; λ, k, T) = T·x^{k/2}/λ^{k/2-1}·exp(-T/2·(x²+λ²))·I_{k/2-1}(x·λ·T)
    %
    % Parámetros: λ ≥ 0 (norma de medias), k ≥ 1 (dimensión efectiva),
    % T > 0 (precisión / escala, T = 1/σ²). Port de DistributionsPuig.jl.
    properties
        lam (1,1) double
        k   (1,1) double
        T   (1,1) double
    end

    methods
        function obj = PuigDistribution(lam, k, T)
            % PuigDistribution(λ, k, T) construye la distribución.
            arguments
                lam (1,1) double {mustBeNonnegative}
                k   (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
                T   (1,1) double {mustBePositive}
            end
            obj.lam = lam;
            obj.k = k;
            obj.T = T;
        end
    end
end