function out = marcumq(a, b, m)
    % MARCUMQ  Función Marcum-Q generalizada Q_m(a, b) = ccdf(NoncentralChisq(2m, a²), b²).
    % b puede ser escalar o vector. Si commuta a la cola asintótica (b-a ≥ 4 y
    % ab ≥ 30) usa marcumq_asymp. Port de marcumq de src/puig_pdf.jl.
    arguments
        a (1,1) double
        b double
        m (1,1) double = 1.0
    end

    scalar = isscalar(b);
    b = b(:);
    out = zeros(size(b));
    for i = 1:numel(b)
        out(i) = marcumq_scalar(a, b(i), m);
    end
    if scalar
        out = out(1);
    end
end

function v = marcumq_scalar(a, b, m)
    if (b - a) >= 4.0 && a * b >= 30.0
        v = puigdist.marcumq_asymp(a, b, m);
    else
        v = 1.0 - ncx2cdf(b * b, 2.0 * m, a * a);
    end
end