"""Candidate A: is the galaxy relation (RAR) a relaxation attractor?
If so, its INTRINSIC scatter should be smaller for more evolved galaxies
(gas-poor, high surface brightness) than for less evolved ones (gas-rich,
low surface brightness). A fundamental force law predicts no difference.
Uses all 165 SPARC galaxies, the McGaugh+2016 RAR with a0 = 1.2e-10, and
stellar M/L 0.5 (disc) / 0.7 (bulge). Intrinsic scatter = observed scatter
minus measurement scatter (velocity errors + 5% floor), in quadrature."""
import sys
sys.path.insert(0, '.')
import numpy as np
from scipy.stats import spearmanr
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G, rar_nu

A0 = 1.2e-10 * 3.0857e13
rows = []
for g in load_sparc():
    vb2 = v_baryon_sq(g, 0.5, 0.7)
    ok = (vb2 > 0) & (g.Vobs > 0)
    if ok.sum() < 5:
        continue
    gb = vb2[ok] / g.R[ok]
    go = g.Vobs[ok] ** 2 / g.R[ok]
    res = np.log10(go) - np.log10(gb * rar_nu(gb / A0))
    err = 2 * g.sigma[ok] / g.Vobs[ok] / np.log(10)      # dex in g_obs
    Ld, Lb = g.luminosity()
    ms = 0.5 * Ld + 0.7 * Lb
    mg = max(float(g.Vgas[-1] * abs(g.Vgas[-1]) * g.R[-1] / G), 0.0)
    fgas = mg / (mg + ms)
    sb0 = np.max(g.SBdisk[:3]) if len(g.SBdisk) else np.nan      # central disc surface brightness
    rows.append(dict(name=g.name, res=res, err=err, off=np.median(res), fgas=fgas,
                     lsb=np.log10(max(sb0, 1e-3)), lm=np.log10(ms + mg)))

def intrinsic(sel):
    """Galaxy-level: scatter of per-galaxy offsets, minus their errors;
    point-level: residual scatter minus measurement scatter."""
    res = np.concatenate([r["res"] for r in sel]); err = np.concatenate([r["err"] for r in sel])
    tot = np.var(res); meas = np.mean(err ** 2)
    offs = np.array([r["off"] for r in sel])
    e_off = np.array([np.median(r["err"]) / np.sqrt(len(r["res"])) for r in sel])
    gtot = np.var(offs); gmeas = np.mean(e_off ** 2)
    return (np.sqrt(max(tot - meas, 0)), np.sqrt(tot), np.sqrt(meas),
            np.sqrt(max(gtot - gmeas, 0)), len(sel))

def boot(sel, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    return np.array([intrinsic([sel[i] for i in rng.integers(0, len(sel), len(sel))])[3] for _ in range(n)])

print(f"{len(rows)} galaxies\n")
for key, lab in [("fgas", "gas fraction"), ("lsb", "log central surface brightness"), ("lm", "log baryonic mass")]:
    v = np.array([r[key] for r in rows]); med = np.median(v)
    lo = [r for r in rows if r[key] < med]; hi = [r for r in rows if r[key] >= med]
    a, b = intrinsic(lo), intrinsic(hi)
    ba, bb = boot(lo), boot(hi)
    d = bb - ba
    print(f"split by {lab} (median {med:.2f}):")
    print(f"   low  half: galaxy-offset intrinsic scatter {a[3]:.3f} dex (68%: {np.percentile(ba,16):.3f}-{np.percentile(ba,84):.3f}); point scatter total {a[1]:.3f}, measurement {a[2]:.3f}, intrinsic {a[0]:.3f}")
    print(f"   high half: galaxy-offset intrinsic scatter {b[3]:.3f} dex (68%: {np.percentile(bb,16):.3f}-{np.percentile(bb,84):.3f}); point scatter total {b[1]:.3f}, measurement {b[2]:.3f}, intrinsic {b[0]:.3f}")
    print(f"   difference (high - low): {np.median(d):+.3f} dex, 95% CI {np.percentile(d,2.5):+.3f}..{np.percentile(d,97.5):+.3f}")
    rho, p = spearmanr(v, [abs(r['off']) for r in rows])
    print(f"   |galaxy offset| vs {lab}: Spearman rho {rho:+.2f} (p {p:.3f})\n")

print("=== Within-galaxy scatter (each galaxy's own offset removed): robust to galaxy-wide distance/inclination errors ===")
def within(sel):
    res = np.concatenate([r["res"] - r["off"] for r in sel]); err = np.concatenate([r["err"] for r in sel])
    n_eff = np.concatenate([np.full(len(r["res"]), (len(r["res"]) - 1) / len(r["res"])) for r in sel])
    tot = np.mean(res ** 2); meas = np.mean(err ** 2 * n_eff)
    return np.sqrt(max(tot - meas, 0)), np.sqrt(tot), np.sqrt(meas)
def bootw(sel, n=2000, seed=1):
    rng = np.random.default_rng(seed)
    return np.array([within([sel[i] for i in rng.integers(0, len(sel), len(sel))])[0] for _ in range(n)])
for key, lab in [("fgas", "gas fraction"), ("lsb", "central surface brightness"), ("lm", "baryonic mass")]:
    v = np.array([r[key] for r in rows]); med = np.median(v)
    lo = [r for r in rows if r[key] < med]; hi = [r for r in rows if r[key] >= med]
    a, b = within(lo), within(hi); d = bootw(hi) - bootw(lo)
    print(f"{lab:28s}: low half intrinsic {a[0]:.3f} (total {a[1]:.3f}, meas {a[2]:.3f}) | high half intrinsic {b[0]:.3f} (total {b[1]:.3f}, meas {b[2]:.3f}) | diff {np.median(d):+.3f}, 95% CI {np.percentile(d,2.5):+.3f}..{np.percentile(d,97.5):+.3f}")
