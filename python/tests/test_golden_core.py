"""Verificación de paridad del núcleo Python contra los golden files de Julia.

Uso:  python tests/test_golden_core.py
"""

import sys
import os
import csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from puigdist.core import (
    Puig_pdf, Puig_surviving, Puig_cumulative, Puig_moments, Puig_quantile,
)
from puigdist.mb import MB_pdf, MB_surviving, MB_cumulative

GOLDEN = os.path.join(os.path.dirname(__file__), "..", "..", "golden")

PASS = 0
FAIL = 0


def check(label, python_val, golden_val, rel_tol=1e-9, abs_tol=1e-14):
    global PASS, FAIL
    g = float(golden_val)
    p = float(python_val)
    if g == 0 and p == 0:
        PASS += 1
        return
    if g == 0:
        err = abs(p - g)
        ok = err < abs_tol
    elif abs(g) > 1e-6:
        err = abs(p - g) / abs(g)
        ok = err < rel_tol
    else:
        err = abs(p - g)
        ok = err < abs_tol
    if ok:
        PASS += 1
    else:
        FAIL += 1
        print("  FAIL %s: python=%.16g, golden=%.16g, err=%.2e" % (label, p, g, err))


def rows_of(name):
    with open(os.path.join(GOLDEN, name), encoding="utf-8") as f:
        for row in csv.DictReader(f):
            yield row


def test_pdf():
    print("=== PDF ===")
    for row in rows_of("golden_pdf.csv"):
        lam = float(row["\u03bb"].strip())
        kv = float(row["k"])
        Tv = float(row["T"])
        xv = float(row["x"])
        g_asymp = float(row["pdf_asymp"])
        g_arb = float(row["pdf_arb"])
        p_asymp = Puig_pdf(xv, lam, kv, Tv, method="asymp")
        check("pdf_asymp lam=%g k=%g T=%g x=%g" % (lam, kv, Tv, xv),
              p_asymp, g_asymp, rel_tol=1e-10, abs_tol=1e-15)
        if lam > 0:
            p_arb = Puig_pdf(xv, lam, kv, Tv, method="arb")
            check("pdf_arb lam=%g k=%g T=%g x=%g" % (lam, kv, Tv, xv),
                  p_arb, g_arb, abs_tol=1e-12)


def test_surv():
    print("=== SURV ===")
    for row in rows_of("golden_surv.csv"):
        lam = float(row["\u03bb"].strip())
        kv = float(row["k"])
        Tv = float(row["T"])
        xv = float(row["x"])
        g = float(row["su"])
        p = Puig_surviving([xv], lam, kv, Tv)[0]
        # 5e-8: margen para diferencias scipy ncx2.sf vs Distributions.jl
        check("surv lam=%g k=%g T=%g x=%g" % (lam, kv, Tv, xv),
              p, g, rel_tol=5e-8, abs_tol=1e-14)


def test_cdf():
    print("=== CDF ===")
    for row in rows_of("golden_cdf.csv"):
        lam = float(row["\u03bb"].strip())
        kv = float(row["k"])
        Tv = float(row["T"])
        xv = float(row["x"])
        g = float(row["cu"])
        p = Puig_cumulative([xv], lam, kv, Tv)[0]
        check("cdf lam=%g k=%g T=%g x=%g" % (lam, kv, Tv, xv),
              p, g, rel_tol=5e-8, abs_tol=1e-14)


def test_moments():
    print("=== MOMENTS ===")
    for row in rows_of("golden_moments.csv"):
        lam = float(row["\u03bb"].strip())
        kv = float(row["k"])
        Tv = float(row["T"])
        n = int(row["n"])
        g = float(row["mu"])
        mu = Puig_moments(lam, kv, Tv, n)
        p = float(mu[n - 1])
        check("moments lam=%g k=%g T=%g n=%d" % (lam, kv, Tv, n),
              p, g, rel_tol=1e-10)


def test_quantile():
    print("=== QUANTILE ===")
    for row in rows_of("golden_quantile.csv"):
        lam = float(row["\u03bb"].strip())
        kv = float(row["k"])
        Tv = float(row["T"])
        pval = float(row["p"])
        g = float(row["q"])
        q = Puig_quantile(pval, lam, kv, Tv)
        check("quantile lam=%g k=%g T=%g p=%g" % (lam, kv, Tv, pval),
              q, g, rel_tol=1e-7, abs_tol=1e-10)


def test_mb():
    print("=== MB ===")
    for row in rows_of("golden_mb.csv"):
        kv = float(row["k"])
        Bv = float(row["B"])
        xv = float(row["x"])
        gp, gs, gc = float(row["mb_pdf"]), float(row["mb_surv"]), float(row["mb_cdf"])
        pp = MB_pdf([xv], kv, Bv)[0]
        ps = MB_surviving([xv], kv, Bv)[0]
        pc = MB_cumulative([xv], kv, Bv)[0]
        check("mb_pdf k=%g B=%g x=%g" % (kv, Bv, xv), pp, gp, rel_tol=1e-10)
        check("mb_surv k=%g B=%g x=%g" % (kv, Bv, xv), ps, gs, rel_tol=1e-10)
        check("mb_cdf k=%g B=%g x=%g" % (kv, Bv, xv), pc, gc, rel_tol=1e-10)
    print("=== MB MOMENTS ===")
    for row in rows_of("golden_mb_moments.csv"):
        kv = float(row["k"])
        Bv = float(row["B"])
        n = int(row["n"])
        g = float(row["mb_mu"])
        mu = Puig_moments(0.0, kv, 1.0 / Bv, n)
        p = float(mu[n - 1])
        check("mb_mom k=%g B=%g n=%d" % (kv, Bv, n), p, g, rel_tol=1e-10)


if __name__ == "__main__":
    test_pdf()
    test_surv()
    test_cdf()
    test_moments()
    test_quantile()
    test_mb()
    print("=" * 60)
    print("RESULTADOS: %d PASS, %d FAIL" % (PASS, FAIL))
    if FAIL:
        sys.exit(1)
    print("Todos los tests pasaron (OK)")