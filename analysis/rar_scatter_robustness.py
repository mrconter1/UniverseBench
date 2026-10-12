"""How robust is the observed intrinsic between-galaxy RAR scatter (0.058 dex,
quality-1, analysis/rar_age_v2.py)? Same estimator (ML tau with each galaxy's
offset error budget), varying the assumptions one at a time:
error-budget scalings, stellar M/L uncertainty, a 10% baryonic-velocity
systematic (as tested by Desmond 2023, arXiv:2303.11314), a0, M/L values,
and subsamples."""
import sys
sys.path.insert(0, ".")
import numpy as np
from scipy.optimize import minimize_scalar
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import rar_nu

KPC = 3.0857e13
T1 = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    if len(f) != 19 or not f[1].isdigit():
        continue
    T1[f[0]] = dict(D=float(f[2]), eD=float(f[3]), fD=int(f[4]), inc=float(f[5]), einc=float(f[6]),
                    L=float(f[7]), MHI=float(f[13]), Q=int(f[17]))
GALS = [g for g in load_sparc() if g.name in T1 and T1[g.name]["inc"] >= 30]


def rows(a0=1.2e-10, ud=0.5, ub=0.7, kD=1.0, kI=1.0, eml=0.1, vsys=0.0, kstat=1.0):
    A0 = a0 * KPC
    out = []
    for g in GALS:
        t = T1[g.name]
        vb2 = v_baryon_sq(g, ud, ub)
        vstar2 = ud * g.Vdisk * np.abs(g.Vdisk) + ub * g.Vbul * np.abs(g.Vbul)
        ok = (vb2 > 0) & (g.Vobs > 0)
        if ok.sum() < 5:
            continue
        gb, go = vb2[ok] / g.R[ok], g.Vobs[ok] ** 2 / g.R[ok]
        off = np.median(np.log10(go) - np.log10(gb * rar_nu(gb / A0)))
        stat = kstat * np.median(2 * g.sigma[ok] / g.Vobs[ok] / np.log(10)) / np.sqrt(ok.sum())
        e_dist = kD * t["eD"] / t["D"] / np.log(10)
        e_inc = kI * 2 / np.tan(np.radians(t["inc"])) * np.radians(t["einc"]) / np.log(10)
        fstar = np.median(np.clip(vstar2[ok] / vb2[ok], 0, 1))
        slope = 1 - 0.5 / (1 + np.sqrt(np.median(gb / A0)))
        e_ml = eml * fstar * slope
        e_vb = 2 * vsys / np.log(10) * slope          # 10% in V_bar -> 0.087 dex in g_bar
        out.append(dict(off=off, err=np.sqrt(stat**2 + e_dist**2 + e_inc**2 + e_ml**2 + e_vb**2),
                        Q=t["Q"], fD=t["fD"], name=g.name))
    return out


def tau(sel):
    o = np.array([r["off"] for r in sel]); e = np.array([r["err"] for r in sel])
    def nll(tau):
        v = e**2 + tau**2
        mu = np.sum(o / v) / np.sum(1 / v)
        return 0.5 * np.sum((o - mu) ** 2 / v + np.log(v))
    return minimize_scalar(nll, bounds=(0, 0.5), method="bounded").x


def boot(sel, n=400):
    rng = np.random.default_rng(0)
    t = [tau([sel[i] for i in rng.integers(0, len(sel), len(sel))]) for _ in range(n)]
    return np.percentile(t, [16, 84])


q1 = lambda R: [r for r in R if r["Q"] == 1]
base = q1(rows())
lo, hi = boot(base)
print(f"baseline, quality 1 (n={len(base)}): tau = {tau(base):.3f} dex (68% bootstrap {lo:.3f}-{hi:.3f})")
cases = [("distance errors x0.5", dict(kD=0.5)), ("distance errors x1.5", dict(kD=1.5)),
         ("inclination errors x0.5", dict(kI=0.5)), ("inclination errors x1.5", dict(kI=1.5)),
         ("M/L uncertainty 0.05 dex", dict(eml=0.05)), ("M/L uncertainty 0.15 dex", dict(eml=0.15)),
         ("+10% baryonic-velocity systematic", dict(vsys=0.10)), ("+5% baryonic-velocity systematic", dict(vsys=0.05)),
         ("velocity errors x1.5", dict(kstat=1.5)),
         ("a0 = 1.0e-10", dict(a0=1.0e-10)), ("a0 = 1.4e-10", dict(a0=1.4e-10)),
         ("M/L disc 0.4, bulge 0.56", dict(ud=0.4, ub=0.56)), ("M/L disc 0.6, bulge 0.84", dict(ud=0.6, ub=0.84))]
for lab, kw in cases:
    print(f"  {lab:36s}: tau = {tau(q1(rows(**kw))):.3f}")
allq = rows()
print(f"all qualities (n={len(allq)}): tau = {tau(allq):.3f}")
direct = [r for r in allq if r["fD"] in (2, 3, 4, 5)]
print(f"direct distances (TRGB, Cepheids, cluster, SN; n={len(direct)}): tau = {tau(direct):.3f}")
dq1 = [r for r in direct if r["Q"] == 1]
print(f"quality 1 AND direct distances (n={len(dq1)}): tau = {tau(dq1):.3f} (68% {boot(dq1)[0]:.3f}-{boot(dq1)[1]:.3f})")
q2 = [r for r in allq if r["Q"] == 2]
print(f"quality 2 only (n={len(q2)}): tau = {tau(q2):.3f}")
