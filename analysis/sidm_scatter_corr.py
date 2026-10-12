"""Does correlated halo-baryon scatter rescue the tightness of the RAR?
Idea 1 drew halo concentration offsets at random (0.11 dex) and got 0.10 dex
of relation scatter against ~0.06-0.07 observed. Here each galaxy's
concentration offset is correlated (coefficient rho) with its disc
compactness: the residual of log R_disk (SPARC Table 1) at fixed baryonic
mass, sign chosen so that rho > 0 means compact discs sit in MORE
concentrated halos (the simulation claim). Halo-mass scatter 0.15 dex kept
independent. Run for the frozen cross-section and for 3x it."""
import sys, importlib.util
sys.path.insert(0, '.')
import numpy as np
from scipy.optimize import minimize_scalar
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G

spec = importlib.util.spec_from_file_location("csidm", "submissions/solver_claude/my_model.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
T1 = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    if len(f) == 19 and f[1].isdigit():
        T1[f[0]] = dict(inc=float(f[5]), Q=int(f[17]), L=float(f[7]), Rd=float(f[11]), MHI=float(f[13]))
GALS = [g for g in load_sparc() if g.name in T1 and T1[g.name]["Q"] == 1 and T1[g.name]["inc"] >= 30]
u = m.UPS_DISK
MB = np.array([m._lum_total(g)[0] * u + 1.4 * u * m._lum_total(g)[1] + m._gas_mass(g) for g in GALS])
lR = np.log10([T1[g.name]["Rd"] for g in GALS])
lM = np.log10(MB)
coef = np.polyfit(lM, lR, 1)
dR = lR - np.polyval(coef, lM)                       # + : more extended than typical
comp = -(dR - dR.mean()) / dR.std()                  # standardised compactness (+ = compact)
print(f"{len(GALS)} galaxies; size-mass slope {coef[0]:.2f}, size scatter at fixed mass {dR.std():.2f} dex")

def run(rho, seed, sig_c=0.11, sig_m=0.15):
    rng = np.random.default_rng(seed)
    gb_all, gm_all, gi = [], [], []
    for j, g in enumerate(GALS):
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        mh = m.halo_mass(MB[j]) * 10 ** (sig_m * rng.standard_normal())
        dc = sig_c * (rho * comp[j] + np.sqrt(1 - rho ** 2) * rng.standard_normal())
        c = m.concentration(mh) * 10 ** dc
        h = m.SIDMHalo(mh, c, m._mb_profile(g, u))
        v2 = np.clip(vb2 + G * h.mass(g.R) / g.R, 0, None)
        ok = vb2 > 0
        gb_all += list(vb2[ok] / g.R[ok]); gm_all += list(v2[ok] / g.R[ok]); gi += [j] * int(ok.sum())
    gb, gm, gi = np.array(gb_all), np.array(gm_all), np.array(gi)
    rar = lambda a0: gb / (-np.expm1(-np.sqrt(gb / a0)))
    a0 = 10 ** minimize_scalar(lambda la: np.mean((np.log10(gm) - np.log10(rar(10 ** la))) ** 2),
                               bounds=(-1, 4), method="bounded").x
    res = np.log10(gm) - np.log10(rar(a0))
    return np.std([np.median(res[gi == j]) for j in range(len(GALS))])

base = m.SIGMA0
for lab, s0 in [("frozen cross-section", base), ("cross-section x3", base * 3)]:
    m.SIGMA0 = s0
    print(f"\n{lab}:")
    for rho in (-0.9, -0.5, 0.0, 0.5, 0.9, 1.0):
        vals = [run(rho, s) for s in (11, 12)]
        print(f"   rho {rho:+.1f}: relation scatter {np.mean(vals):.3f} dex  (draws {', '.join(f'{v:.3f}' for v in vals)})")
    vals = [run(1.0, s, sig_m=0.0) for s in (11, 12)]
    print(f"   rho +1.0, no halo-mass scatter: {np.mean(vals):.3f} dex")
    vals = [run(0.0, s, sig_c=0.0, sig_m=0.0) for s in (11,)]
    print(f"   no halo scatter at all:         {np.mean(vals):.3f} dex")
print("\ntarget (observed intrinsic): 0.058 dex all quality-1; 0.070 with direct distances")
