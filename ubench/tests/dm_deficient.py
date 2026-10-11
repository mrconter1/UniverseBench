"""
Test 10: galaxies apparently lacking dark matter (CONTESTED).

NGC 1052-DF2 and DF4 are faint, puffy ("ultra-diffuse") galaxies whose
stars move so slowly that they seem to contain little or no dark matter
(Danieli et al. 2019; van Dokkum et al. 2019). That is awkward for standard
dark matter (galaxies should sit in halos) and for MOND (the boost should
apply to everyone, though the host's external field softens it). The
distances, and with them the masses, are disputed, so the test is scored
but never used to fit constants.

The galaxies are passed as dwarf-like objects, so models reuse
predict_dwarf_dispersion from test 05. Host field: g_ext = V^2 / D with
V ~ 200 km/s for NGC 1052 (approximate) and the projected distances (lower
limits on the true 3D distances).

Uncertainty: the dispersion error (asymmetric for DF4) plus the stellar
mass error propagated through each model's prediction.
"""
import copy
import csv

import numpy as np

from ..core import Test, data_path

V_HOST_NGC1052 = 200.0


class UDG:
    def __init__(self, r):
        self.name, self.host = r["name"], "NGC1052"
        self.M_bar = float(r["M_star"])
        self.log_Mb_err = float(r["log_M_err"])
        self.r_half = float(r["r_half_kpc"])
        self.D_host = float(r["D_host_kpc"])
        self.g_ext = V_HOST_NGC1052**2 / self.D_host
        self.sigma_obs = float(r["sigma"])
        self.err_lo, self.err_hi = float(r["sigma_err_lo"]), float(r["sigma_err_hi"])


class DMDeficientGalaxies(Test):
    ID = "10_dm_deficient"
    LABEL = "Galaxies lacking dark matter"
    STATUS = "contested"
    HIDDEN_SPLIT = False
    USED_IN_FIT = False
    DATA_DIRS = ("dm_deficient",)

    def load(self):
        path = data_path("dm_deficient", "df2_df4.csv")
        self.galaxies = [UDG(r) for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_dwarf_dispersion"):
            return None
        z = []
        for g in self.galaxies:
            pred = model.predict_dwarf_dispersion(g, gp)
            gu = copy.copy(g)
            gu.M_bar = g.M_bar * 10 ** g.log_Mb_err
            dm = model.predict_dwarf_dispersion(gu, gp) - pred
            err = g.err_hi if pred > g.sigma_obs else g.err_lo
            z.append((pred - g.sigma_obs) / np.hypot(err, dm))
        return np.array(z)
