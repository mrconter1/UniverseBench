"""
Test 09: the Bullet Cluster (1E0657-558).

Two clusters collided. Their hot gas, which holds most of the visible
matter, was slowed by the collision and sits between them, while the
galaxies passed through. Weak lensing (Clowe et al. 2006, ApJ 648, L109)
shows the mass peaks lying with the galaxies, offset from the gas at about
8 sigma.

Simplified test of that argument. For each of the two sub-clusters we ask
where its lensing mass peak lies on the line from its galaxies (brightest
cluster galaxy, BCG) to its gas peak, as a fraction: 0 = at the galaxies,
1 = at the gas. Observed: 0 for both, with uncertainty 1/8 of the
separation (the 8 sigma offset from the gas). Inputs per sub-cluster from
Table 2 (100 kpc apertures): gas mass at the plasma peak and stellar mass
at the BCG.

Caveat: real analyses fit full lensing maps; a mass-weighted position is
a crude stand-in, kept deliberately simple.

Model method: predict_bullet_offset(s, gp) -> fraction in [0, 1], given
s.M_gas, s.M_star [Msun], s.sep_kpc.
"""
import csv
import math

import numpy as np

from ..core import Test, data_path

KPC_PER_ARCSEC = 4.41          # at z = 0.296 (flat LCDM, h = 0.7)
OBS, OBS_ERR = 0.0, 1 / 8.0


def _radec(ra, dec):
    h, m, s = (float(x) for x in ra.split(":"))
    d, dm, ds = (abs(float(x)) for x in dec.split(":"))
    sign = -1 if dec.strip().startswith("-") else 1
    return 15 * (h + m / 60 + s / 3600), sign * (d + dm / 60 + ds / 3600)


class SubCluster:
    def __init__(self, name, bcg, plasma):
        self.name = name
        self.M_gas = float(plasma["M_gas"]) * 1e12
        self.M_star = float(bcg["M_star"]) * 1e12
        (a1, d1), (a2, d2) = _radec(bcg["ra"], bcg["dec"]), _radec(plasma["ra"], plasma["dec"])
        sep_arcsec = 3600 * math.hypot((a1 - a2) * math.cos(math.radians(d1)), d1 - d2)
        self.sep_kpc = sep_arcsec * KPC_PER_ARCSEC


class BulletCluster(Test):
    ID = "09_bullet_cluster"
    LABEL = "Bullet Cluster"
    SIMPLIFIED = ("mass-weighted peak position from aperture masses, not a fit to the "
                  "lensing map; a faithful model should compute where its lensing mass sits")
    HIDDEN_SPLIT = False
    DATA_DIRS = ("bullet_cluster",)

    def load(self):
        path = data_path("bullet_cluster", "clowe2006.csv")
        rows = {r["component"]: r for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))}
        self.subs = [SubCluster("main", rows["main_BCG"], rows["main_plasma"]),
                     SubCluster("sub", rows["sub_BCG"], rows["sub_plasma"])]

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_bullet_offset"):
            return None
        return np.array([(model.predict_bullet_offset(s, gp) - OBS) / OBS_ERR for s in self.subs])
