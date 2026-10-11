"""
Tests 01 and 02: galaxy rotation curves (SPARC).

Data: Lelli, McGaugh & Schombert (2016), AJ 152, 157. For each galaxy and
radius: observed rotation speed Vobs (km/s) with error, and the speeds the
visible gas, stellar disk and bulge would produce alone under Newtonian
gravity (disk and bulge for a mass-to-light ratio of 1). 165 galaxies with
at least 6 points are used.

Uncertainty per point: sigma^2 = errV^2 + (0.05 Vobs)^2; the 5% floor stands
in for distance and inclination errors SPARC's quoted errors leave out.

01 PREDICT: global constants only; predict hidden galaxies' whole curves.
02 EXTRAPOLATE: per-galaxy knobs (LOCAL_PARAMS) are tuned on the inner 60%
   of each hidden galaxy; predict the outer 40%.

Model method: predict(g, R, gp, lp) -> rotation speed [km/s] at radii R [kpc]
"""
import glob
import os

import numpy as np
from scipy.optimize import minimize

from ..core import Test, data_path, default_locals, is_hidden, _vec_to_dict

FLOOR = 0.05
MIN_POINTS = 6


class Galaxy:
    ups_fixed = False

    def __init__(self, path):
        self.name = os.path.basename(path).replace("_rotmod.dat", "")
        with open(path) as f:
            self.distance_mpc = float(f.readline().split("=")[1].split()[0])
        d = np.loadtxt(path, comments="#")
        (self.R, self.Vobs, self.errV, self.Vgas, self.Vdisk,
         self.Vbul, self.SBdisk, self.SBbul) = d.T[:8]
        self.sigma = np.sqrt(self.errV**2 + (FLOOR * self.Vobs) ** 2)

    def luminosity(self):
        """Rough 3.6um luminosities [Lsun] by integrating surface brightness
        over the measured radii (misses the faint outskirts)."""
        R_pc = self.R * 1000.0
        ring = 2 * np.pi * R_pc * np.gradient(R_pc)
        return float(np.sum(ring * self.SBdisk)), float(np.sum(ring * self.SBbul))

    def stellar_mass(self, ups_disk):
        Ld, Lb = self.luminosity()
        return ups_disk * Ld + 1.4 * ups_disk * Lb

    def subset(self, mask):
        g = object.__new__(Galaxy)
        g.__dict__.update(self.__dict__)
        for k in ("R", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul", "sigma"):
            setattr(g, k, getattr(self, k)[mask])
        return g


def load_sparc(folder=None):
    folder = folder or data_path("sparc")
    gals = [Galaxy(p) for p in sorted(glob.glob(os.path.join(folder, "*_rotmod.dat")))]
    return [g for g in gals if len(g.R) >= MIN_POINTS]


def _z(model, g, gp, lp):
    v = model.predict(g, g.R, gp, lp)
    return (np.asarray(v) - g.Vobs) / g.sigma


def fit_local(model, g, gp):
    spec = getattr(model, "LOCAL_PARAMS", {})
    if not spec:
        return {}
    chi2 = lambda x: float(np.sum(_z(model, g, gp, _vec_to_dict(x, spec)) ** 2))
    starts = [[v[0] for v in spec.values()]]
    rng = np.random.default_rng(0)
    starts += [[rng.uniform(lo, hi) for (_, lo, hi) in spec.values()] for _ in range(4)]
    best = None
    for x0 in starts:
        res = minimize(lambda x: chi2(x) if np.isfinite(chi2(x)) else 1e12, x0,
                       method="Nelder-Mead",
                       options={"maxiter": 3000, "xatol": 1e-4, "fatol": 1e-4})
        if best is None or res.fun < best.fun:
            best = res
    return _vec_to_dict(best.x, spec)


class _SparcTest(Test):
    def load(self):
        self.galaxies = load_sparc()

    def part(self, part, salt):
        if part == "all":
            return self.galaxies
        want = part == "hidden"
        return [g for g in self.galaxies if is_hidden(g.name, salt) == want]


class RotationPredict(_SparcTest):
    ID = "01_rotation_predict"
    LABEL = "Rotation curves: predict"

    def export_public(self, dest_data, salt):
        import shutil
        os.makedirs(os.path.join(dest_data, "sparc"), exist_ok=True)
        for g in self.part("public", salt):
            shutil.copy(data_path("sparc", f"{g.name}_rotmod.dat"), os.path.join(dest_data, "sparc"))

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict"):
            return None
        lp = default_locals(model)
        return np.concatenate([_z(model, g, gp, lp) for g in self.part(part, salt)])


class RotationExtrapolate(_SparcTest):
    ID = "02_rotation_extrapolate"
    LABEL = "Rotation curves: extrapolate"
    USED_IN_FIT = False   # the predict test already uses these galaxies

    def residuals(self, model, gp, part, salt, inner_frac=0.6):
        if not hasattr(model, "predict"):
            return None
        zs = []
        n_local = len(getattr(model, "LOCAL_PARAMS", {}))
        for g in self.part(part, salt):
            cut = int(np.ceil(inner_frac * len(g.R)))
            if cut < max(3, n_local + 1) or len(g.R) - cut < 2:
                continue
            idx = np.arange(len(g.R))
            inner, outer = g.subset(idx < cut), g.subset(idx >= cut)
            zs.append(_z(model, outer, gp, fit_local(model, inner, gp)))
        return np.concatenate(zs)
