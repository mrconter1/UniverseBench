"""
Galaxy Benchmark, phase 1: can a theory of gravity / dark matter
predict how fast stars and gas orbit in 175 real galaxies (SPARC)?

Data: SPARC (Lelli, McGaugh & Schombert 2016, AJ 152, 157).
Each galaxy file gives, per radius: observed rotation speed Vobs with
error, and the speeds that the *visible* gas, stellar disk and bulge
alone would produce under Newtonian gravity (disk/bulge for a
mass-to-light ratio of 1).

A model sees only the visible matter and must say how fast things orbit.

Two tests, both scored ONLY on hidden galaxies the model never saw:

  A. PREDICT     - global parameters only (tuned on public galaxies).
                   Pure prediction from visible matter, no per-galaxy knobs.
  B. EXTRAPOLATE - the model may tune its declared per-galaxy parameters
                   on the inner 60% of each hidden galaxy's points, then
                   must predict the outer 40%. Extra knobs only help if
                   they capture real physics, so overfitting is punished
                   automatically.

Per-point uncertainty: sigma^2 = errV^2 + (0.05 * Vobs)^2. The 5% floor
stands in for distance and inclination uncertainties that SPARC's
quoted errors leave out.

  C. CLUSTERS    - total mass of hidden galaxy clusters from their gas and
                   stars (see clusters.py).

Headline per test: percentage of hidden points predicted within 2 sigma
(a perfect model scores about 95). Total = mean of the tests.
Also reported: median reduced chi^2 and the public-vs-hidden gap.

Global constants are fitted jointly on the public data of ALL tests: a
theory gets one value of each constant for galaxies and clusters alike.
"""
import glob
import hashlib
import os

import numpy as np
from scipy.optimize import minimize

import clusters

HERE = os.path.dirname(os.path.abspath(__file__))
FLOOR = 0.05
MIN_POINTS = 6
# The secret that decides which galaxies are hidden. Set BENCH_SALT to
# your own private string before giving an agent its public bundle.
SPLIT_SALT = os.environ.get("BENCH_SALT", "demo-salt")


class Galaxy:
    def __init__(self, path):
        self.name = os.path.basename(path).replace("_rotmod.dat", "")
        with open(path) as f:
            self.distance_mpc = float(f.readline().split("=")[1].split()[0])
        d = np.loadtxt(path, comments="#")
        (self.R, self.Vobs, self.errV, self.Vgas, self.Vdisk,
         self.Vbul, self.SBdisk, self.SBbul) = d.T[:8]
        self.sigma = np.sqrt(self.errV**2 + (FLOOR * self.Vobs) ** 2)

    def luminosity(self):
        """Rough 3.6um luminosities [Lsun] by integrating surface
        brightness over the measured radii (misses the faint outskirts)."""
        R_pc = self.R * 1000.0
        dR = np.gradient(R_pc)
        ring = 2 * np.pi * R_pc * dR
        return float(np.sum(ring * self.SBdisk)), float(np.sum(ring * self.SBbul))

    def subset(self, mask):
        g = object.__new__(Galaxy)
        g.__dict__.update(self.__dict__)
        for k in ("R", "Vobs", "errV", "Vgas", "Vdisk", "Vbul",
                  "SBdisk", "SBbul", "sigma"):
            setattr(g, k, getattr(self, k)[mask])
        return g


def load_all(folder=None):
    folder = folder or os.path.join(HERE, "data", "all")
    gals = [Galaxy(p) for p in sorted(glob.glob(os.path.join(folder, "*_rotmod.dat")))]
    return [g for g in gals if len(g.R) >= MIN_POINTS]


def is_hidden(name, salt=SPLIT_SALT):
    return int(hashlib.sha256((salt + name).encode()).hexdigest(), 16) % 2 == 1


def split(gals, salt=SPLIT_SALT):
    pub = [g for g in gals if not is_hidden(g.name, salt)]
    hid = [g for g in gals if is_hidden(g.name, salt)]
    return pub, hid


# ---------------------------------------------------------------- fitting
def _vec_to_dict(vec, spec):
    return {k: float(np.clip(v, lo, hi)) for v, (k, (_, lo, hi)) in zip(vec, spec.items())}


def _chi2(model, g, gp, lp):
    v = model.predict(g, g.R, gp, lp)
    if not np.all(np.isfinite(v)):
        return 1e12
    return float(np.sum(((v - g.Vobs) / g.sigma) ** 2))


def fit_global(model, gals, cl=None):
    """Tune global constants on PUBLIC data of ALL tests at once: one set of
    constants must serve galaxies and clusters alike. Each test contributes
    its mean chi^2 per point, so tests weigh equally regardless of size."""
    spec = model.GLOBAL_PARAMS
    if not spec:
        return {}
    x0 = [v[0] for v in spec.values()]
    lp_default = {k: v[0] for k, v in model.LOCAL_PARAMS.items()}
    n_rot = sum(len(g.R) for g in gals)

    def f(x):
        gp = _vec_to_dict(x, spec)
        total = sum(_chi2(model, g, gp, lp_default) for g in gals) / n_rot
        if cl:
            z = clusters.residuals(model, cl, gp)
            if z is not None:
                total += float(np.mean(z**2)) if np.all(np.isfinite(z)) else 1e12
        return total

    res = minimize(f, x0, method="Nelder-Mead",
                   options={"maxiter": 6000, "xatol": 1e-5, "fatol": 1e-6})
    return _vec_to_dict(res.x, spec)


def fit_local(model, g, gp):
    spec = model.LOCAL_PARAMS
    if not spec:
        return {}
    best = None
    # a few starts to avoid bad local minima
    starts = [[v[0] for v in spec.values()]]
    rng = np.random.default_rng(0)
    for _ in range(4):
        starts.append([rng.uniform(lo, hi) for (_, lo, hi) in spec.values()])
    for x0 in starts:
        res = minimize(lambda x: _chi2(model, g, gp, _vec_to_dict(x, spec)),
                       x0, method="Nelder-Mead",
                       options={"maxiter": 3000, "xatol": 1e-4, "fatol": 1e-4})
        if best is None or res.fun < best.fun:
            best = res
    return _vec_to_dict(best.x, spec)


# ---------------------------------------------------------------- scoring
def _summ(z, n_gal):
    z = np.concatenate(z) if z else np.array([])
    return {
        "within_2sigma_pct": round(100 * float(np.mean(np.abs(z) < 2)), 1) if len(z) else 0.0,
        "reduced_chi2_median": round(float(np.median(z**2)), 2) if len(z) else None,
        "reduced_chi2_mean": round(float(np.mean(z**2)), 2) if len(z) else None,
        "n_points": int(len(z)),
        "n_galaxies": n_gal,
    }


def test_predict(model, gals, gp):
    zs = []
    for g in gals:
        lp = {k: v[0] for k, v in model.LOCAL_PARAMS.items()}
        zs.append((model.predict(g, g.R, gp, lp) - g.Vobs) / g.sigma)
    return _summ(zs, len(gals)), zs


def test_extrapolate(model, gals, gp, inner_frac=0.6):
    zs, n = [], 0
    for g in gals:
        cut = int(np.ceil(inner_frac * len(g.R)))
        n_local = len(model.LOCAL_PARAMS)
        if cut < max(3, n_local + 1) or len(g.R) - cut < 2:
            continue
        inner = g.subset(np.arange(len(g.R)) < cut)
        outer = g.subset(np.arange(len(g.R)) >= cut)
        lp = fit_local(model, inner, gp)
        zs.append((model.predict(outer, outer.R, gp, lp) - outer.Vobs) / outer.sigma)
        n += 1
    return _summ(zs, n), zs


def split_clusters(pts, salt=SPLIT_SALT):
    """Both radii of a cluster always land on the same side."""
    return ([c for c in pts if not is_hidden(c.name, salt)],
            [c for c in pts if is_hidden(c.name, salt)])


def test_clusters(model, pts, gp):
    z = clusters.residuals(model, pts, gp)
    if z is None:  # model has no cluster prediction: scores zero
        return _summ([np.full(len(pts), np.inf)], len({c.name for c in pts}))
    out = _summ([z], len({c.name for c in pts}))
    out["mean_ratio_pred_over_obs"] = round(float(np.mean(
        [model.predict_cluster_mass(c, gp) / c.M_tot for c in pts])), 3)
    return out


TESTS = ("A_predict", "B_extrapolate", "C_clusters")


def run(model, gals=None, cl=None, salt=SPLIT_SALT, verbose=True):
    gals = gals or load_all()
    cl = cl if cl is not None else clusters.load_clusters()
    pub, hid = split(gals, salt)
    cpub, chid = split_clusters(cl, salt)
    gp = fit_global(model, pub, cpub)
    res = {
        "A_predict": test_predict(model, hid, gp)[0],
        "B_extrapolate": test_extrapolate(model, hid, gp)[0],
        "C_clusters": test_clusters(model, chid, gp),
    }
    total = round(float(np.mean([res[t]["within_2sigma_pct"] for t in TESTS])), 1)
    out = {
        "model": model.NAME,
        "global_params": {k: round(v, 4) for k, v in gp.items()},
        "n_global_params": len(model.GLOBAL_PARAMS),
        "n_local_params": len(model.LOCAL_PARAMS),
        **{t + "_hidden": res[t] for t in TESTS},
        "A_predict_public": test_predict(model, pub, gp)[0],
        "C_clusters_public": test_clusters(model, cpub, gp),
        "TOTAL": total,
    }
    if verbose:
        s = "  ".join(f"{t.split('_')[0]}={res[t]['within_2sigma_pct']:5.1f}" for t in TESTS)
        print(f"{model.NAME:30s} {s}  TOTAL={total:5.1f}   {out['global_params']}")
    return out
