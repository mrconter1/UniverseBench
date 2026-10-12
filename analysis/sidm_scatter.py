"""Idea 1: can collapsing SIDM reproduce the observed tightness of the RAR?
Observed: intrinsic between-galaxy scatter tau = 0.058 dex for SPARC
quality-1 galaxies (analysis/rar_age_v2.py). Here the frozen cSIDM model is
evaluated on the same galaxies with REALISTIC halo diversity: per-galaxy
log-normal scatter in halo mass at fixed baryonic mass (0.15 dex) and in
concentration (0.11 dex, Dutton & Maccio 2014). For each realisation we fit
one RAR a0 to the model's accelerations and measure the scatter of
per-galaxy offsets. Same draws for plain NFW (cold dark matter) halos, to see
how much self-interactions erase halo diversity."""
import sys, importlib.util
sys.path.insert(0, '.')
import numpy as np
from scipy.optimize import minimize_scalar
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G, nfw_r200

spec = importlib.util.spec_from_file_location("csidm", "submissions/solver_claude/my_model.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
T1 = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    if len(f) == 19 and f[1].isdigit():
        T1[f[0]] = dict(inc=float(f[5]), Q=int(f[17]))
GALS = [g for g in load_sparc() if g.name in T1 and T1[g.name]["Q"] == 1 and T1[g.name]["inc"] >= 30]
M = m.CollapsingSIDM()
u = m.UPS_DISK

def mb_total(g):
    Ld, Lb = m._lum_total(g)
    return u * Ld + 1.4 * u * Lb + m._gas_mass(g)

def run(dm_scatter, dc_scatter, seed, sidm=True):
    rng = np.random.default_rng(seed)
    gb_all, gm_all, gal_idx = [], [], []
    for j, g in enumerate(GALS):
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        mh = m.halo_mass(mb_total(g)) * 10 ** (dm_scatter * rng.standard_normal())
        c = m.concentration(mh) * 10 ** (dc_scatter * rng.standard_normal())
        if sidm:
            h = m.SIDMHalo(mh, c, m._mb_profile(g, u))
            mdm = h.mass(g.R)
        else:
            mdm = m.nfw_m(g.R, mh, c, nfw_r200(mh))
        v2 = np.clip(vb2 + G * mdm / g.R, 0, None)
        ok = vb2 > 0
        gb_all += list(vb2[ok] / g.R[ok]); gm_all += list(v2[ok] / g.R[ok]); gal_idx += [j] * int(ok.sum())
    gb, gm, gi = np.array(gb_all), np.array(gm_all), np.array(gal_idx)
    rar = lambda a0: gb / (-np.expm1(-np.sqrt(gb / a0)))
    a0 = 10 ** minimize_scalar(lambda la: np.mean((np.log10(gm) - np.log10(rar(10 ** la))) ** 2),
                               bounds=(-1, 4), method="bounded").x       # in (km/s)^2/kpc units
    res = np.log10(gm) - np.log10(rar(a0))
    offs = np.array([np.median(res[gi == j]) for j in range(len(GALS)) if np.any(gi == j)])
    within = np.sqrt(np.mean([np.mean((res[gi == j] - np.median(res[gi == j])) ** 2) for j in range(len(GALS)) if np.any(gi == j)]))
    run.pooled = float(np.std(res))      # point-by-point pooled width (what simulation papers quote)
    return np.std(offs), within, a0 * 1e6 / 3.0857e19     # (km/s)^2/kpc -> m/s^2

print(f"{len(GALS)} quality-1 galaxies (inclination >= 30)")
print("observed (data): between-galaxy intrinsic scatter 0.058 dex\n")
for lab, sidm in [("cSIDM", True), ("plain NFW (CDM)", False)]:
    s0, w0, a0 = run(0, 0, 0, sidm)
    print(f"{lab:16s} no halo scatter:        between {s0:.3f} dex, within-galaxy {w0:.3f}, a0 {a0:.2e}")
    vals, pooled = [], []
    for s in (1, 2, 3):
        vals.append(run(0.15, 0.11, s, sidm)); pooled.append(run.pooled)
    print(f"{lab:16s} realistic halo scatter: pooled point-by-point width {np.mean(pooled):.3f} dex")
    print(f"{lab:16s} realistic halo scatter: between {np.mean([v[0] for v in vals]):.3f} dex "
          f"(realisations {', '.join(f'{v[0]:.3f}' for v in vals)}), within-galaxy {np.mean([v[1] for v in vals]):.3f}")
    vals = [run(0.0, 0.11, s, sidm) for s in (1, 2)]
    print(f"{lab:16s} concentration only:     between {np.mean([v[0] for v in vals]):.3f} dex\n")
