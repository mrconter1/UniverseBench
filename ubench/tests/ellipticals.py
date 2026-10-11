"""
Test 08: massive elliptical galaxies via strong gravitational lensing.

Data: 70 early-type lens galaxies from the SLACS survey (Auger et al. 2009,
ApJ 705, 1099, Table 4; lenses with complete entries). Each bends light from
a background galaxy into a ring of radius r_Ein, which precisely measures the
projected total mass inside it, M_Ein. Verified with the table's internal
Salpeter/Chabrier ratio (~1.78 in every row).

Input: projected stellar mass inside r_Ein for a Chabrier stellar
population (f*_Chab x M_Ein: photometry times stellar mass-to-light, not a
dynamical quantity) and the galaxy's total stellar mass (Chabrier). The
stellar population of giant ellipticals may be heavier than Chabrier
(Salpeter-like), so models should carry a global stellar-mass factor for
ellipticals; baselines fit one between 0.8 and 2.0.

Uncertainty: stellar-mass error propagated through each model's prediction,
plus a 5% floor on the lensing mass.

Model method: predict_lens_mass(L, gp) -> projected total mass [Msun] inside
L.r_Ein [kpc], given L.Mstar_proj_chab (stars inside r_Ein), L.Mstar_chab
(total stars), L.mstar_rel_err.
"""
import copy
import csv

import numpy as np

from ..core import Test, data_path, is_hidden

FLOOR = 0.05


class Lens:
    def __init__(self, r):
        self.name = r["name"]
        self.r_Ein = float(r["r_Ein_kpc"])
        self.M_Ein = 10 ** float(r["log_M_Ein"])
        self.Mstar_proj_chab = float(r["fstar_Chab"]) * self.M_Ein
        self.Mstar_chab = 10 ** float(r["log_Mstar_Chab"])
        self.mstar_rel_err = float(r["fstar_Chab_err"]) / float(r["fstar_Chab"])


class Ellipticals(Test):
    ID = "08_ellipticals"
    LABEL = "Elliptical galaxies (lensing)"

    def load(self):
        path = data_path("ellipticals", "auger2009_slacs.csv")
        self.lenses = [Lens(r) for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]

    def export_public(self, dest_data, salt):
        import os
        os.makedirs(os.path.join(dest_data, "ellipticals"), exist_ok=True)
        with open(data_path("ellipticals", "auger2009_slacs.csv")) as fi, \
                open(os.path.join(dest_data, "ellipticals", "auger2009_slacs.csv"), "w") as fo:
            for line in fi:
                if line.startswith(("#", "name,")) or not is_hidden(line.split(",")[0], salt):
                    fo.write(line)

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_lens_mass"):
            return None
        ls = self.lenses if part == "all" else \
            [l for l in self.lenses if is_hidden(l.name, salt) == (part == "hidden")]
        z = []
        for L in ls:
            m = model.predict_lens_mass(L, gp)
            Lu = copy.copy(L)
            Lu.Mstar_proj_chab *= 1 + L.mstar_rel_err
            Lu.Mstar_chab *= 1 + L.mstar_rel_err
            dm = model.predict_lens_mass(Lu, gp) - m
            z.append((m - L.M_Ein) / np.sqrt(dm**2 + (FLOOR * L.M_Ein) ** 2))
        return np.array(z)
