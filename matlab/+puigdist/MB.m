classdef MB
    % Caso límite Maxwell-Boltzmann (λ → 0) de la distribución de Puig.
    %
    %   f(x; k, B) = 2^{1-k/2}·B^{-k/2}/Γ(k/2)·x^{k-1}·exp(-x²/(2B))
    %
    % Con B = 1/T (escala térmica / dispersión). Port de DistributionsPuig.jl.
    properties
        k (1,1) double
        B (1,1) double
    end

    methods
        function obj = MB(k, B)
            % MB(k, B) construye la distribución límite.
            arguments
                k (1,1) double {mustBeGreaterThanOrEqual(k, 1)}
                B (1,1) double {mustBePositive}
            end
            obj.k = k;
            obj.B = B;
        end
    end
end