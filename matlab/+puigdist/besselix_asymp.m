function res = besselix_asymp(nu, z)
    % BESSELIX_ASYMP  e^{-z}·I_ν(z) por expansión asintótica de 5 términos (DLMF 10.40.1).
    % Precisión relativa < 1e-12 para z ≥ 100 y < 1e-15 para z ≥ 200.
    % Port de besselix_asymp de src/puig_pdf.jl (DistributionsPuig.jl).
    mu = 4.0 * nu .* nu;
    res = ones(size(z));
    t1 = -(mu - 1.0) ./ (8.0 * z);
    res = res + t1;
    t2 = -t1 .* (mu - 9.0) ./ (16.0 * z);
    res = res + t2;
    t3 = -t2 .* (mu - 25.0) ./ (24.0 * z);
    res = res + t3;
    t4 = -t3 .* (mu - 49.0) ./ (32.0 * z);
    res = res + t4;
    t5 = -t4 .* (mu - 81.0) ./ (40.0 * z);
    res = res + t5;
    res = res ./ sqrt(2.0 * pi .* z);
end