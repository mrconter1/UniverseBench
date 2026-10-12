"""Candidate A, decisive version: between-galaxy scatter of the RAR with each
galaxy's own distance, inclination and M/L uncertainties (SPARC Table 1)
included in its offset error budget (Li et al. 2018 logic, linearised):
  distance D -> D(1+d): g_obs ~ 1/D, g_bar ~ D^0   => offset error e_D/D / ln10
  inclination i:        g_obs ~ 1/sin^2 i          => 2 cot(i) e_i / ln10
  stellar M/L 0.1 dex:  shifts g_bar on the stellar share; the RAR response
                        d log g / d log g_bar (0.5 deep, 1 high) maps it to g.
Intrinsic between-galaxy scatter = ML estimate of tau in
offset_j ~ N(mu, sqrt(err_j^2 + tau^2)). Compared for halves split by
evolution proxies. A relaxation attractor predicts larger tau for less
evolved galaxies; a force law predicts equal tau."""
import sys
sys.path.insert(0, '.')
import numpy as np
from scipy.optimize import minimize_scalar
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G, rar_nu

A0 = 1.2e-10 * 3.0857e13
T1 = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    # data rows: name, T, D, e_D, f_D, Inc, e_Inc, L, e_L, Reff, SBeff, Rdisk, SBdisk,
    #            MHI, RHI, Vflat, e_Vflat, Q, Ref  (whitespace-separated; robust to spacing)
    if len(f) != 19 or not f[1].isdigit():
        continue
    try:
        T1[f[0]] = dict(T=int(f[1]), D=float(f[2]), eD=float(f[3]), inc=float(f[5]), einc=float(f[6]),
                        L=float(f[7]), SBd=float(f[12]), MHI=float(f[13]), Q=int(f[17]))
    except ValueError:
        continue
assert len(T1) == 175, len(T1)

rows = []
for g in load_sparc():
    t = T1.get(g.name)
    if t is None or t["inc"] < 30:            # Li+2018 / McGaugh+2016 cut
        continue
    vb2 = v_baryon_sq(g, 0.5, 0.7)
    vstar2 = 0.5 * g.Vdisk * np.abs(g.Vdisk) + 0.7 * g.Vbul * np.abs(g.Vbul)
    ok = (vb2 > 0) & (g.Vobs > 0)
    if ok.sum() < 5:
        continue
    gb = vb2[ok] / g.R[ok]
    go = g.Vobs[ok] ** 2 / g.R[ok]
    res = np.log10(go) - np.log10(gb * rar_nu(gb / A0))
    off = np.median(res)
    stat = np.median(2 * g.sigma[ok] / g.Vobs[ok] / np.log(10)) / np.sqrt(ok.sum())
    e_dist = t["eD"] / t["D"] / np.log(10)
    e_inc = 2 / np.tan(np.radians(t["inc"])) * np.radians(t["einc"]) / np.log(10)
    fstar = np.median(np.clip(vstar2[ok] / vb2[ok], 0, 1))
    y = np.median(gb / A0)
    slope = 1 - 0.5 / (1 + np.sqrt(y))        # d log g / d log g_bar, 0.5 deep -> 1 high
    e_ml = 0.1 * fstar * slope
    err = np.sqrt(stat ** 2 + e_dist ** 2 + e_inc ** 2 + e_ml ** 2)
    mgas = 1.33 * t["MHI"]
    fgas = mgas / (mgas + 0.5 * t["L"])
    rows.append(dict(name=g.name, off=off, err=err, fgas=fgas, lsb=np.log10(max(t["SBd"], 1e-3)),
                     lm=np.log10(1e9 * (0.5 * t["L"] + mgas)), T=t["T"], Q=t["Q"],
                     e_dist=e_dist, e_inc=e_inc))

def tau_ml(sel):
    o = np.array([r["off"] for r in sel]); e = np.array([r["err"] for r in sel])
    def nll(tau):
        v = e ** 2 + tau ** 2
        mu = np.sum(o / v) / np.sum(1 / v)
        return 0.5 * np.sum((o - mu) ** 2 / v + np.log(v))
    return minimize_scalar(nll, bounds=(0, 0.5), method="bounded").x

def boot(sel, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    return np.array([tau_ml([sel[i] for i in rng.integers(0, len(sel), len(sel))]) for _ in range(n)])

print(f"{len(rows)} galaxies (inclination >= 30 deg). Median offset error budget: "
      f"distance {np.median([r['e_dist'] for r in rows]):.3f}, inclination {np.median([r['e_inc'] for r in rows]):.3f}, "
      f"total {np.median([r['err'] for r in rows]):.3f} dex")
print(f"all galaxies: intrinsic between-galaxy scatter tau = {tau_ml(rows):.3f} dex  (raw scatter of offsets {np.std([r['off'] for r in rows]):.3f})\n")
for key, lab, less_evolved_is_high in [("fgas", "gas fraction", True), ("lsb", "disc central surface brightness", False),
                                       ("lm", "baryonic mass", False), ("T", "Hubble type (late = high)", True)]:
    v = np.array([r[key] for r in rows]); med = np.median(v)
    lo = [r for r in rows if r[key] < med]; hi = [r for r in rows if r[key] >= med]
    less, more = (hi, lo) if less_evolved_is_high else (lo, hi)
    tl, tm = tau_ml(less), tau_ml(more)
    d = boot(less) - boot(more, seed=1)
    print(f"{lab:32s}: less evolved tau {tl:.3f} (n={len(less)}) | more evolved tau {tm:.3f} (n={len(more)}) | "
          f"difference {np.median(d):+.3f}, 95% CI {np.percentile(d,2.5):+.3f}..{np.percentile(d,97.5):+.3f}")
q1 = [r for r in rows if r["Q"] == 1]
print(f"\nquality-1 galaxies only (n={len(q1)}): tau = {tau_ml(q1):.3f} dex")
print("\n=== Same splits within quality-1 galaxies only (removes the data-quality confounder) ===")
for key, lab, less_evolved_is_high in [("fgas", "gas fraction", True), ("lsb", "disc central surface brightness", False),
                                       ("lm", "baryonic mass", False), ("T", "Hubble type (late = high)", True)]:
    v = np.array([r[key] for r in q1]); med = np.median(v)
    lo = [r for r in q1 if r[key] < med]; hi = [r for r in q1 if r[key] >= med]
    less, more = (hi, lo) if less_evolved_is_high else (lo, hi)
    d = boot(less) - boot(more, seed=1)
    print(f"{lab:32s}: less evolved tau {tau_ml(less):.3f} (n={len(less)}) | more evolved tau {tau_ml(more):.3f} (n={len(more)}) | "
          f"difference {np.median(d):+.3f}, 95% CI {np.percentile(d,2.5):+.3f}..{np.percentile(d,97.5):+.3f}")
from collections import Counter
print("\nquality flag by mass half:", "low mass", Counter(r["Q"] for r in rows if r["lm"] < np.median([x['lm'] for x in rows])),
      "| high mass", Counter(r["Q"] for r in rows if r["lm"] >= np.median([x['lm'] for x in rows])))
