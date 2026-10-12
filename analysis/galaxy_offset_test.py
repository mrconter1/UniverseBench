"""A sharper rotation-curve test: whole-galaxy offsets with each galaxy's own
error budget. Test 01 treats points as independent with a 5% floor, but the
dominant errors (distance, inclination, M/L) move a galaxy's WHOLE curve, so
every model scores ~55-58 there and the test barely discriminates.
Here: for each galaxy, offset = median over radii of log10(g_obs / g_model),
error = statistical + distance + inclination + 0.1 dex M/L (as in
rar_age_v2.py). z = offset / error; % of galaxies with |z| < 2. A perfect
model scores ~95. No per-galaxy freedom."""
import sys, importlib.util
sys.path.insert(0, ".")
import numpy as np
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import rar_nu

T1 = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    if len(f) == 19 and f[1].isdigit():
        T1[f[0]] = dict(D=float(f[2]), eD=float(f[3]), inc=float(f[5]), einc=float(f[6]), Q=int(f[17]))
GALS = [g for g in load_sparc() if g.name in T1 and T1[g.name]["inc"] >= 30]


def load(path):
    spec = importlib.util.spec_from_file_location(path, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


gf, cs = load("submissions/solver_galform/my_model.py"), load("submissions/solver_claude/my_model.py")
GP_GF = {**{k: v[0] for k, v in gf.MODEL.GLOBAL_PARAMS.items()},
         "ups_disk": 0.42078, "halo_A": 11.55989, "halo_B": 0.8545, "halo_lfloor": 9.8397, "dlogc": -0.11012}
A0 = 1.2e-10 * 3.0857e13


def models():
    yield "MOND (a0 1.2e-10, M/L 0.5)", 0.5, lambda g, u: np.sqrt(np.clip(v_baryon_sq(g, u, 1.4 * u), 1e-9, None) / g.R
                                                                  * rar_nu(np.clip(v_baryon_sq(g, u, 1.4 * u), 1e-9, None) / g.R / A0) * g.R)
    yield "CDM + galaxy formation (DC14)", GP_GF["ups_disk"], lambda g, u: gf.MODEL.predict(g, g.R, GP_GF, {"ups_scale": 1.0})
    yield "cSIDM (frozen)", cs.UPS_DISK, lambda g, u: cs.MODEL.predict(g, g.R, {}, {"ups_scale": 1.0})


def score(fn, u, sel):
    zs = []
    for g in sel:
        t = T1[g.name]
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        ok = (vb2 > 0) & (g.Vobs > 0)
        if ok.sum() < 5:
            continue
        v = np.asarray(fn(g, u))
        off = np.median(np.log10(g.Vobs[ok] ** 2 / v[ok] ** 2))
        vstar2 = u * g.Vdisk * np.abs(g.Vdisk) + 1.4 * u * g.Vbul * np.abs(g.Vbul)
        fstar = np.median(np.clip(vstar2[ok] / vb2[ok], 0, 1))
        stat = np.median(2 * g.sigma[ok] / g.Vobs[ok] / np.log(10)) / np.sqrt(ok.sum())
        e = np.sqrt(stat ** 2 + (t["eD"] / t["D"] / np.log(10)) ** 2
                    + (2 / np.tan(np.radians(t["inc"])) * np.radians(t["einc"]) / np.log(10)) ** 2
                    + (0.1 * fstar * 0.75) ** 2)
        zs.append(off / e)
    zs = np.array(zs)
    return 100 * np.mean(np.abs(zs) < 2), np.median(zs), len(zs)


for lab, sel in [("quality 1", [g for g in GALS if T1[g.name]["Q"] == 1]),
                 ("quality 1+2", [g for g in GALS if T1[g.name]["Q"] <= 2])]:
    print(f"--- {lab}")
    for name, u, fn in models():
        s, med, n = score(fn, u, sel)
        print(f"  {name:32s}: {s:5.1f}% of {n} galaxies within 2 sigma (median z {med:+.2f})", flush=True)
