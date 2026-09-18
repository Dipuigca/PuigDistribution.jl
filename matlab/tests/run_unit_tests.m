% Tests unitarios del port MATLAB (port de R/tests/run_unit_tests.R).
% Uso: matlab -batch "run('matlab/tests/run_unit_tests.m')"
% (debe ejecutarse desde la raíz del repositorio)

global PASS FAIL
PASS = 0; FAIL = 0;

function ok = check(cond, label)
  global PASS FAIL
  if isscalar(cond) && islogical(cond) && cond
    PASS = PASS + 1;
  else
    FAIL = FAIL + 1;
    fprintf('  FAIL: %s\n', label);
  end
end

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

ok = @(cond, label) check(cond, label); % alias de check

% ---- 0. Tipos ----
d = puigdist.PuigDistribution(3.0, 2.0, 1.0);
check(d.lam == 3.0 && d.k == 2.0 && d.T == 1.0, 'campos d');
mb = puigdist.MB(2.0, 1.0);
check(mb.k == 2.0 && mb.B == 1.0, 'campos MB');

% ---- 1. effective_dimension ----
ok(puigdist.effective_dimension(randn(50, 1)) == 1.0, 'effd 1D');
rng(123);
M_uncorr = randn(30000, 3);
e1 = puigdist.effective_dimension(M_uncorr);
ok(e1 > 2.9 && e1 < 3.1, sprintf('effd incorreladas ~3 (%g)', e1));
xb = randn(100, 1);
M_corr = [xb, xb, xb];
e2 = puigdist.effective_dimension(M_corr);
ok(e2 < 1.0001, sprintf('effd correladas ~1 (%g)', e2));
k_exacto = 3^2 / (3 + 2 * (0.8^2 + 0.6^2 + 0.7^2));
ok(abs(k_exacto - 1.5050167224) < 1e-6, 'k formula exacta');

% ---- 2. empirical_survival ----
t = [1 2 3 4 5]; grid = [0 1 3 5 6];
s = puigdist.empirical_survival(t, grid);
ok(all([s(1) == 1, s(2) == 1, abs(s(3) - 0.6) < 1e-12, abs(s(4) - 0.2) < 1e-12, s(5) == 0]), ...
   'emp_survival valores');

% ---- 3. Métricas ----
S_true = [1 0.8 0.5 0.2 0]; xg = 1:5;
st = puigdist.calc_all_metrics(xg, S_true, S_true);
ok(abs(st.R2 - 1) < 1e-12 && st.MAE == 0 && st.RMSE == 0 && st.MaxAE == 0 && st.IAE == 0, ...
   'metricas perfectas');
st2 = puigdist.calc_all_metrics(xg, S_true, S_true + 0.05);
ok(abs(st2.MAE - 0.05) < 1e-12 && abs(st2.MaxAE - 0.05) < 1e-12 && st2.R2 < 1, ...
   'metricas discrepancia');

% ---- 4. Puig_fit (table / matrix / vector) ----
rng(42);
N = 150;
X1 = 150 + randn(N, 1) * 20;
X2 = 100 + 0.8 * X1 + randn(N, 1) * 10;
X3 = 200 + 0.5 * X1 + randn(N, 1) * 15;
df = table(X1, X2, X3, 'VariableNames', {'Squat', 'Bench', 'Deadlift'});

res_dyn = puigdist.Puig_fit(df, 'method', 'dynamic');
ok(isstruct(res_dyn) && isa(res_dyn.params, 'puigdist.PuigDistribution') ...
   && isfield(res_dyn, 'indicator') && isfield(res_dyn, 'stats') ...
   && isfield(res_dyn, 'method'), 'res es struct con params/indicator/stats');
ok(numel(res_dyn.indicator) == N, 'len indicador');
ok(res_dyn.params.lam > 0 && res_dyn.params.T > 0, 'lam,T > 0');
ok(res_dyn.params.k >= 1 && res_dyn.params.k <= 3, 'k en [1,3]');
ok(res_dyn.stats.R2 > 0.95, 'R2 dynamic > 0.95');

res_fix = puigdist.Puig_fit(df, 'method', 'fixed');
ok(abs(res_fix.params.k - 3) < 1e-9 && res_fix.stats.R2 > 0.90, 'fixed k=3');
res_num = puigdist.Puig_fit(df, 'method', 'numerical');
ok(res_num.params.k >= 1 && res_num.params.k <= 3 && res_num.stats.R2 > 0.95, ...
   'numerical');

res_sub = puigdist.Puig_fit(df, 'vars', {'Squat', 'Bench'});
ok(res_sub.params.k <= 2 && numel(res_sub.indicator) == N, 'vars subset');

df2 = df;
df2.Total = sqrt(X1.^2 + X2.^2 + X3.^2);
res_time = puigdist.Puig_fit(df2, 'vars', {'Squat', 'Bench', 'Deadlift'}, ...
   'time_col', 'Total');
ok(max(abs(res_time.indicator - df2.Total)) < 1e-9, 'time_col indicador');

res_mat = puigdist.Puig_fit(df{:, {'Squat', 'Bench', 'Deadlift'}}, 'method', 'dynamic');
ok(abs(res_mat.params.lam - res_dyn.params.lam) / res_dyn.params.lam < 1e-6 && ...
     abs(res_mat.params.k - res_dyn.params.k) < 1e-6, 'matrix == dynamic');

res_vec = puigdist.Puig_fit(df2.Total);
ok(res_vec.params.lam > 0 && res_vec.params.k >= 1 && res_vec.params.T > 0 && ...
     res_vec.stats.R2 > 0.95, 'vector 1D');

rng(77);
muestra = puigdist.Puig_rand(500, 40.0, 2.5, 0.05);
res_1d = puigdist.Puig_fit(muestra);
ok(abs(res_1d.params.lam - 40) < 5 && res_1d.params.k >= 1 && res_1d.stats.R2 > 0.98, ...
   sprintf('1D sintetico lam=%.2f R2=%.4f', res_1d.params.lam, res_1d.stats.R2));
res_1d_kf = puigdist.Puig_fit(muestra, 'k_fixed', 2.5);
ok(abs(res_1d_kf.params.k - 2.5) < 1e-6 && res_1d_kf.stats.R2 > 0.98, '1D k_fixed');

% ---- 5. Errores ----
try, puigdist.Puig_fit(rand(5, 2)); err_ = false; catch, err_ = true; end
ok(err_, 'fit n<10');
try, puigdist.Puig_fit(rand(5, 1)); err_ = false; catch, err_ = true; end
ok(err_, 'fit vec n<10');
try, puigdist.Puig_fit(df, 'method', 'metodo_inexistente'); err_ = false; catch, err_ = true; end
ok(err_, 'metodo desconocido');
try, puigdist.Puig_fit(rand(20, 1), 'k_fixed', 0.5); err_ = false; catch, err_ = true; end
ok(err_, 'k_fixed < 1');
TblText = table({'a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'i'; 'j'; 'k'});
try, puigdist.Puig_fit(TblText); err_ = false; catch, err_ = true; end
ok(err_, 'sin columnas numeric');

% ---- 6. Interfaz básica ----
ok(puigdist.Puig_pdf(0, 3, 2, 1) == 0, 'pdf(0)=0');
ok(isinf(puigdist.Puig_logpdf(0, 3, 2, 1)), 'logpdf(0)=-Inf');
ok(puigdist.Puig_cumulative(0, 3, 2, 1) == 0 && puigdist.Puig_surviving(0, 3, 2, 1) == 1, ...
   'cdf/surv(0)');
ok(abs(puigdist.Puig_cumulative(100, 3, 2, 1) - 1) < 1e-6, 'cdf(100)~1');
ok(abs(puigdist.Puig_surviving(100, 3, 2, 1)) < 1e-6, 'surv(100)~0');
ok(abs(puigdist.Puig_cumulative(5, 3, 2, 1) + puigdist.Puig_surviving(5, 3, 2, 1) - 1) < 1e-10, ...
   'cdf+surv=1');
ok(puigdist.Puig_quantile(0.01, 3, 2, 1) >= 0, 'quantile>=0');
m = puigdist.Puig_mean(3, 2, 1); v = puigdist.Puig_var(3, 2, 1);
st3 = puigdist.Puig_std(3, 2, 1);
ok(abs(v - st3^2) < 1e-6 && m > 0 && v > 0, 'mean/var/std coherentes');
rng(42);
samples = puigdist.Puig_rand(100, 3, 2, 1);
ok(all(samples >= 0) && abs(mean(samples) - m) < 1.0, 'rand coherente');

% ---- 7. Robustez asymp vs arb ----
for lamv = [10 50 100 250 500 1000]
  for fr = [0.7 0.9 1.0 1.1 1.3]
    xv = fr * lamv;
    va = puigdist.Puig_pdf(xv, lamv, 2.4, 0.01, 'asymp');
    vb = puigdist.Puig_pdf(xv, lamv, 2.4, 0.01, 'arb');
    if vb > 1e-15
      test_ok = abs(va - vb) / vb < 1e-12;
    else
      test_ok = abs(va - vb) < 1e-15;
    end
    ok(test_ok, sprintf('asymp/arb lam=%g frac=%g', lamv, fr));
  end
end
for kv = [1.0001 1.2 2.0 3.0 5.0 10.0]
  for fr = [0.8 1.0 1.2]
    xv = fr * 50;
    va = puigdist.Puig_pdf(xv, 50, kv, 0.05, 'asymp');
    vb = puigdist.Puig_pdf(xv, 50, kv, 0.05, 'arb');
    if vb > 1e-15
      test_ok = abs(va - vb) / vb < 1e-12;
    else
      test_ok = abs(va - vb) < 1e-15;
    end
    ok(test_ok, sprintf('asymp/arb k=%g frac=%g', kv, fr));
  end
end
try, puigdist.Puig_pdf(-5, 100, 2.5, 0.05); err_ = false; catch, err_ = true; end
ok(err_, 'pdf negativo error');

% ---- 8. Momentos ----
t4 = puigdist.Puig_moments(2.5, 3.2, 1.4, 4);
v8 = puigdist.Puig_moments(2.5, 3.2, 1.4, 8);
ok(all(abs(v8(1:4) - t4) < 1e-12), 'mom4 primeros');
ok(numel(v8) == 8, 'mom8 len');
mu1v = v8(1); mu2v = v8(2); mu3v = v8(3); mu4v = v8(4);
lamv = 2.5; kv = 3.2; Tv = 1.4;
mu5_rec = ((10 - 4 + kv) / Tv + lamv^2) * mu3v - 3 * (1 + kv) / Tv^2 * mu1v;
ok(abs(v8(5) - mu5_rec) < 1e-10, 'mu5 recurrencia');
mu6_rec = ((12 - 4 + kv) / Tv + lamv^2) * mu4v - 4 * (2 + kv) / Tv^2 * mu2v;
ok(abs(v8(6) - mu6_rec) < 1e-10, 'mu6 recurrencia');

Bmb = 1.0; kmb = 2.0;
vec_mb = puigdist.Puig_moments(0, kmb, 1 / Bmb, 8);
for nn = 1:8
  mu_an = (2 * Bmb)^(nn / 2) * gamma((kmb + nn) / 2) / gamma(kmb / 2);
  ok(abs(vec_mb(nn) - mu_an) < 1e-10, sprintf('MB mom %d', nn));
end

% Caso regulador: momentos con λ grande (fix de la serie de Kummer).
ml = puigdist.Puig_moments(1000, 3.0, 0.05, 7);
ok(abs(ml(1) - 1000.02) < 0.5 && abs(ml(3) - 1.00012e9) / 1.00012e9 < 1e-6 ...
   && ml(1) > 0, 'mu1(1000,3,0.05)>0 y ~lam');

try, puigdist.Puig_moments(2.5, 3.2, 1.4, 0); err_ = false; catch, err_ = true; end
ok(err_, 'n=0 no valido');
sts = puigdist.Puig_stats(2.5, 3.2, 1.4);
ok(abs(sts.mean - t4(1)) < 1e-12 && abs(sts.var - (t4(2) - t4(1)^2)) < 1e-12, ...
   'stats coherentes');

% ---- 9. Marcum-Q ----
mq = puigdist.marcumq_asymp(10, 15, 1);
ok(mq > 0 && mq < 1, 'marcumq_asymp(10,15,1) en (0,1)');
sTv = sqrt(1); a_t = 10 * sTv;
mq_lo = puigdist.marcumq(a_t, 25 * sTv, 3.2 / 2);
ok(mq_lo < 1e-30, 'cola lejana < 1e-30');
Sext = puigdist.Puig_surviving([20 25 30], 10, 3.2, 1);
ok(all(isfinite(Sext)) && all(Sext >= 0) && issorted(fliplr(Sext(:)')), 'colas monot');
Sbnd = puigdist.Puig_surviving(10 + 4, 10, 3.2, 1);
ok(isfinite(Sbnd) && Sbnd > 0 && Sbnd < 1, 'frontera asymp');
ms = puigdist.marcumq(1, 3, 1);
ok(isfinite(ms) && ms > 0 && ms < 1, 'ab<30 fallback ncx2');
ok(puigdist.marcumq_asymp(10, 20, 2) > puigdist.marcumq_asymp(10, 20, 1), ...
   'm mayor -> mayor Q');

fprintf('\n%s\nRESULTADOS: %d PASS, %d FAIL\n', repmat('=', 1, 60), PASS, FAIL);
if FAIL > 0
  error('Tests unitarios fallados');
end
fprintf('Todos los tests pasaron (OK)\n');