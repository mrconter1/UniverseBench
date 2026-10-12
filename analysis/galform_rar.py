"""What radial acceleration relation does solver model 2 (cold dark matter +
DC14 galaxy-formation halos) produce on SPARC? Noise-free accelerations from
the model with its hidden-run constants, quality-1 galaxies (inc >= 30):
fitted a0 (McGaugh form) and between-galaxy scatter, without and with
realistic halo diversity (0.15 dex halo mass, 0.11 dex concentration at fixed
baryons). Compare: observed intrinsic 0.058 dex (0.03-0.08 under reasonable
error models), a0 = 1.2e-10; cSIDM 0.056 / 0.10; NFW 0.068 / 0.10 (sidm_scatter.py)."""
import sys, importlib.util
sys.path.insert(0, ".")
import numpy as np
from scipy.optimize import minimize_scalar
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G, nfw_concentration

spec = importlib.util.spec_from_file_location("gf", "submissions/solver_galform/my_model.py")
gf = importlib.util.module_from_spec(spec); spec.loader.exec_module(gf)
GP = {"ups_disk": 0.42078, "halo_A": 11.55989, "halo_B": 0.8545, "halo_lfloor": 9.8397, "dlogc": -0.11012}
T1 = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    if len(f) == 19 and f[1].isdigit():
        T1[f[0]] = dict(inc=float(f[5]), Q=int(f[17]))
GALS = [g for g in load_sparc() if g.name in T1 and T1[g.name]["Q"] == 1 and T1[g.name]["inc"] >= 30]
u = GP["ups_disk"]


def run(dm, dc, seed, contract=True):
    rng = np.random.default_rng(seed)
    gb_all, gm_all, gi = [], [], []
    for j, g in enumerate(GALS):
        Ld, Lb = gf._lum_total(g)
        mstar = u * Ld + 1.4 * u * Lb
        mb = mstar + gf._gas_mass(g)
        mh = (10 ** GP["halo_A"] * (mb / 1e10) ** GP["halo_B"] + 10 ** GP["halo_lfloor"]) * 10 ** (dm * rng.standard_normal())
        c = nfw_concentration(mh) * 10 ** GP["dlogc"] * 10 ** (dc * rng.standard_normal())
        X = np.log10(mstar / mh) if contract else -4.5
        h = gf.DC14Halo(mh, c, X)
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        ok = vb2 > 0
        v2 = vb2 + G * h.mass(g.R) / g.R
        gb_all += list(vb2[ok] / g.R[ok]); gm_all += list(v2[ok] / g.R[ok]); gi += [j] * int(ok.sum())
    gb, gm, gi = map(np.array, (gb_all, gm_all, gi))
    rar = lambda a0: gb / (-np.expm1(-np.sqrt(gb / a0)))
    la = minimize_scalar(lambda la: np.mean(np.log10(gm / rar(10 ** la)) ** 2), bounds=(-1, 4), method="bounded").x
    res = np.log10(gm / rar(10 ** la))
    offs = np.array([np.median(res[gi == j]) for j in range(len(GALS))])
    return offs.std(), res.std(), 10 ** la * 1e6 / 3.0857e19


print(f"{len(GALS)} quality-1 galaxies; observed: a0 = 1.2e-10, intrinsic between-galaxy scatter 0.058 (0.03-0.08)")
for lab, kw in [("DC14 (galaxy formation)", {}), ("same halos, plain NFW shape", {"contract": False})]:
    b, p, a0 = run(0, 0, 0, **kw)
    print(f"{lab:30s} no halo scatter:   a0 = {a0:.2e}, between {b:.3f}, pooled {p:.3f}")
    v = [run(0.15, 0.11, s, **kw) for s in (1, 2, 3)]
    print(f"{lab:30s} realistic scatter: a0 = {np.mean([x[2] for x in v]):.2e}, between {np.mean([x[0] for x in v]):.3f}, "
          f"pooled {np.mean([x[1] for x in v]):.3f}")
    v = [run(0.08, 0.11, s, **kw) for s in (1, 2, 3)]
    print(f"{lab:30s} 0.08 dex halo mass scatter: between {np.mean([x[0] for x in v]):.3f}")


def offsets_vs_X():
    rows = []
    for g in GALS:
        Ld, Lb = gf._lum_total(g)
        mstar = u * Ld + 1.4 * u * Lb
        mb = mstar + gf._gas_mass(g)
        mh = 10 ** GP["halo_A"] * (mb / 1e10) ** GP["halo_B"] + 10 ** GP["halo_lfloor"]
        c = nfw_concentration(mh) * 10 ** GP["dlogc"]
        X = np.log10(mstar / mh)
        vb2 = v_baryon_sq(g, u, 1.4 * u); ok = vb2 > 0
        gb = vb2[ok] / g.R[ok]
        out = []
        for Xi in (X, -4.5):
            gm = (vb2[ok] + G * gf.DC14Halo(mh, c, Xi).mass(g.R[ok]) / g.R[ok]) / g.R[ok]
            out.append(np.median(np.log10(gm / (gb / (-np.expm1(-np.sqrt(gb / (1.2e-10 * 3.0857e13))))))))
        rows.append((X, out[0], out[1], np.log10(mb)))
    r = np.array(rows)
    print("\nOffset from the observed RAR (a0 = 1.2e-10), noise-free, by X = log M*/M_halo:")
    for lo, hi in [(-6, -3.5), (-3.5, -2.7), (-2.7, -2.2), (-2.2, -1.8), (-1.8, -1.5), (-1.5, 0)]:
        s = (r[:, 0] >= lo) & (r[:, 0] < hi)
        if s.any():
            print(f"  X in [{lo:5.1f},{hi:5.1f}): n={s.sum():3d}  DC14 offset {np.median(r[s,1]):+.3f}  NFW offset {np.median(r[s,2]):+.3f}  "
                  f"log Mb {np.median(r[s,3]):.1f}")


offsets_vs_X()
