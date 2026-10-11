"""
Test 06: wide binary stars (CONTESTED).

Pairs of stars thousands of AU apart orbit each other under extremely weak
gravity, around the MOND acceleration scale. Two analyses of the same Gaia
DR3 data disagree:
  * Chae 2023: gravity is boosted by 1.43 +/- 0.06 at g_N = 10^-10.15 m/s^2,
    as MOND (AQUAL) with the Galaxy's external field predicts.
  * Banik et al. 2024: alpha_grav = -0.021 (+0.065/-0.045), i.e. Newtonian;
    MOND excluded at 16-19 sigma. Converted to a boost of 0.99 +/- 0.02
    (see data/wide_binaries/measurements.csv for the conversion).

No theory can match both, so this test caps everyone at 50 until the
dispute is resolved. As a contested test it is scored but never used to
fit constants.

The binaries sit in the Milky Way's field at the Sun: g_ext = V0^2 / R0 with
V0 = 229 km/s, R0 = 8.122 kpc (Eilers+ 2019), about 2.1e-10 m/s^2.

Model method: predict_wide_binary_boost(gN_si, g_ext_si, gp) -> g_obs / g_N
"""
import csv

import numpy as np

from ..core import Test, data_path
from ..physics import KPC_M

G_EXT_SUN_SI = 229.0**2 / 8.122 * 1e6 / KPC_M      # (km/s)^2/kpc -> m/s^2


class WideBinaries(Test):
    ID = "06_wide_binaries"
    LABEL = "Wide binary stars"
    STATUS = "contested"
    HIDDEN_SPLIT = False
    USED_IN_FIT = False
    DATA_DIRS = ("wide_binaries",)

    def load(self):
        path = data_path("wide_binaries", "measurements.csv")
        self.rows = list(csv.DictReader(l for l in open(path) if not l.startswith("#")))

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_wide_binary_boost"):
            return None
        return np.array([(model.predict_wide_binary_boost(float(r["gN_si"]), G_EXT_SUN_SI, gp)
                          - float(r["boost"])) / float(r["boost_err"]) for r in self.rows])
