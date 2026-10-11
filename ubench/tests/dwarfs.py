"""
Test 05: dwarf spheroidal galaxies.

Data: 25 dwarf spheroidal satellites of the Milky Way (19) and Andromeda
(6) from McGaugh & Wolf 2010, ApJ 722, 248, Tables 1-2: distance from the
host, line-of-sight velocity dispersion (V_c / sqrt 3), baryonic (stellar)
mass and 3D half-light radius. Verified against the table's own checksum
column F_b = M_b / (45 V_c^4).

These are the most dark-matter-dominated galaxies known. In MOND they are
also where the "external field effect" bites: the host galaxy's pull can
suppress the MOND boost. To let any model use it, each dwarf carries the
host's gravitational acceleration at its distance, g_ext = V_host^2 / D,
with observed host circular speeds of ~190 km/s (Milky Way) and ~230 km/s
(Andromeda) at these distances (approximate).

Uncertainty: the dispersion error plus the effect of the 0.4 dex stellar
mass uncertainty, propagated through each model's own prediction (the test
re-predicts with the mass raised by its error).

Caveat: several ultra-faint dwarfs may be tidally disturbed (McGaugh & Wolf
discuss this), which affects both camps' predictions.

Model method: predict_dwarf_dispersion(d, gp) -> line-of-sight velocity
dispersion [km/s], given d.M_bar [Msun], d.r_half [kpc, 3D],
d.D_host [kpc], d.g_ext [(km/s)^2/kpc], d.host ("MW" or "M31").
"""
import copy
import csv

import numpy as np

from ..core import Test, data_path, is_hidden

V_HOST = {"MW": 190.0, "M31": 230.0}


class Dwarf:
    def __init__(self, r):
        self.name, self.host = r["name"], r["host"]
        self.D_host = float(r["D"])
        self.M_bar = 10 ** float(r["log_Mb"])
        self.log_Mb_err = float(r["log_Mb_err"])
        self.r_half = float(r["r_half"]) / 1000.0
        self.sigma_obs = float(r["Vc"]) / np.sqrt(3)
        self.sigma_err = float(r["Vc_err"]) / np.sqrt(3)
        self.g_ext = V_HOST[self.host] ** 2 / self.D_host


def load_dwarfs():
    path = data_path("dwarfs", "mcgaugh_wolf2010.csv")
    return [Dwarf(r) for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]


class DwarfSpheroidals(Test):
    ID = "05_dwarfs"
    LABEL = "Dwarf galaxies"

    def load(self):
        self.dwarfs = load_dwarfs()

    def export_public(self, dest_data, salt):
        import os
        os.makedirs(os.path.join(dest_data, "dwarfs"), exist_ok=True)
        with open(data_path("dwarfs", "mcgaugh_wolf2010.csv")) as fi, \
                open(os.path.join(dest_data, "dwarfs", "mcgaugh_wolf2010.csv"), "w") as fo:
            for line in fi:
                if line.startswith(("#", "name,")) or not is_hidden(line.split(",")[0], salt):
                    fo.write(line)

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_dwarf_dispersion"):
            return None
        ds = self.dwarfs if part == "all" else \
            [d for d in self.dwarfs if is_hidden(d.name, salt) == (part == "hidden")]
        z = []
        for d in ds:
            pred = model.predict_dwarf_dispersion(d, gp)
            d_up = copy.copy(d)
            d_up.M_bar = d.M_bar * 10 ** d.log_Mb_err
            pred_up = model.predict_dwarf_dispersion(d_up, gp)
            sig = np.hypot(d.sigma_err, pred_up - pred)
            z.append((pred - d.sigma_obs) / sig)
        return np.array(z)
