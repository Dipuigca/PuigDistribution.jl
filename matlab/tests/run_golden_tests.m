% Validación del port MATLAB contra los golden files de Julia.
% Uso: matlab -batch "run('matlab/tests/run_golden_tests.m')"
% (debe ejecutarse desde la raíz del repositorio)

global PASS FAIL
PASS = 0; FAIL = 0;

function ok = check(label, v, gv, rel_tol, abs_tol)
  global PASS FAIL
  p = double(v); g = double(gv);
  if isinf(p) || isinf(g)
    ok = (p == g);
  elseif g == 0 && p == 0
    ok = true;
  elseif g == 0
    ok = abs(p) < abs_tol;
  elseif abs(g) > 1e-6
    ok = abs(p - g) / abs(g) < rel_tol;
  else
    ok = abs(p - g) < abs_tol;
  end
  if ok
    PASS = PASS + 1;
  else
    FAIL = FAIL + 1;
    fprintf('  FAIL %s: M=% .16g golden=% .16g err=%.2e\n', label, p, g, abs(p - g));
  end
end

gpath = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'golden');
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

fprintf('=== PDF ===\n');
pdf = readmatrix(fullfile(gpath, 'golden_pdf.csv'));
for i = 1:size(pdf, 1)
  lamv = pdf(i,1); kv = pdf(i,2); Tv = pdf(i,3); xv = pdf(i,4);
  pA = puigdist.Puig_pdf(xv, lamv, kv, Tv, 'asymp');
  check(sprintf('pdf_asymp lam=%g k=%g T=%g x=%g', lamv, kv, Tv, xv), pA, pdf(i,5), 1e-10, 1e-15);
  if lamv > 0
    pB = puigdist.Puig_pdf(xv, lamv, kv, Tv, 'arb');
    check(sprintf('pdf_arb lam=%g k=%g T=%g x=%g', lamv, kv, Tv, xv), pB, pdf(i,6), 1e-10, 1e-12);
  end
end

fprintf('=== SURV ===\n');
surv = readmatrix(fullfile(gpath, 'golden_surv.csv'));
for i = 1:size(surv, 1)
  lamv = surv(i,1); kv = surv(i,2); Tv = surv(i,3); xv = surv(i,4);
  p = puigdist.Puig_surviving(xv, lamv, kv, Tv);
  check(sprintf('surv lam=%g k=%g T=%g x=%g', lamv, kv, Tv, xv), p, surv(i,5), 5e-8, 1e-14);
end

fprintf('=== CDF ===\n');
cdf = readmatrix(fullfile(gpath, 'golden_cdf.csv'));
for i = 1:size(cdf, 1)
  lamv = cdf(i,1); kv = cdf(i,2); Tv = cdf(i,3); xv = cdf(i,4);
  p = puigdist.Puig_cumulative(xv, lamv, kv, Tv);
  check(sprintf('cdf lam=%g k=%g T=%g x=%g', lamv, kv, Tv, xv), p, cdf(i,5), 5e-8, 1e-14);
end

fprintf('=== MOMENTS ===\n');
mom = readmatrix(fullfile(gpath, 'golden_moments.csv'));
for i = 1:size(mom, 1)
  lamv = mom(i,1); kv = mom(i,2); Tv = mom(i,3); nv = int32(mom(i,4));
  mu = puigdist.Puig_moments(lamv, kv, Tv, nv);
  check(sprintf('moments lam=%g k=%g T=%g n=%d', lamv, kv, Tv, nv), mu(nv), mom(i,5), 1e-10, 1e-12);
end

fprintf('=== QUANTILE ===\n');
qg = readmatrix(fullfile(gpath, 'golden_quantile.csv'));
for i = 1:size(qg, 1)
  lamv = qg(i,1); kv = qg(i,2); Tv = qg(i,3); pv = qg(i,4);
  q = puigdist.Puig_quantile(pv, lamv, kv, Tv);
  check(sprintf('quantile lam=%g k=%g T=%g p=%g', lamv, kv, Tv, pv), q, qg(i,5), 1e-7, 1e-10);
end

fprintf('=== MB ===\n');
mb = readmatrix(fullfile(gpath, 'golden_mb.csv'));
for i = 1:size(mb, 1)
  kv = mb(i,1); Bv = mb(i,2); xv = mb(i,3);
  check(sprintf('mb_pdf k=%g B=%g x=%g', kv, Bv, xv), ...
    puigdist.MB_pdf(xv, kv, Bv), mb(i,4), 1e-10, 1e-15);
  check(sprintf('mb_surv k=%g B=%g x=%g', kv, Bv, xv), ...
    puigdist.MB_surviving(xv, kv, Bv), mb(i,5), 1e-10, 1e-15);
  check(sprintf('mb_cdf k=%g B=%g x=%g', kv, Bv, xv), ...
    puigdist.MB_cumulative(xv, kv, Bv), mb(i,6), 1e-10, 1e-15);
end

fprintf('=== MB MOMENTS ===\n');
mbm = readmatrix(fullfile(gpath, 'golden_mb_moments.csv'));
for i = 1:size(mbm, 1)
  kv = mbm(i,1); Bv = mbm(i,2); nv = int32(mbm(i,3));
  mu = puigdist.Puig_moments(0.0, kv, 1.0 / Bv, nv);
  check(sprintf('mb_mom k=%g B=%g n=%d', kv, Bv, nv), mu(nv), mbm(i,4), 1e-10, 1e-12);
end

fprintf('\n%s\nRESULTADOS: %d PASS, %d FAIL\n', repmat('=', 1, 60), PASS, FAIL);
if FAIL > 0
  error('Tests golden fallados');
end
fprintf('Todos los tests pasaron (OK)\n');