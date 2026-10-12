"""
Test 05: dwarf spheroidal galaxies.

Data: 78 dispersion-supported dwarfs from the Local Volume Database
(Pace 2024, OJAp 7, 111; CC0): 40 satellites of the Milky Way, 26 of
Andromeda and 12 isolated field dwarfs. For each: stellar line-of-sight
velocity dispersion (or a 95% upper limit for 9 of them), stellar mass,
3D half-light radius, distance to the host (or nearest giant galaxy), and
HI mass where detected. Built by data/dwarfs/build_lvdb.py, whose
docstring lists the selection: tidally disrupted systems are dropped,
upper limits are kept, and there is no cut on the dispersion itself.

(Until Oct 2026 this test used the 25 dwarfs of McGaugh & Wolf 2010, kept
in data/dwarfs/archive/. Five of those dispersions differ from modern
values by more than 2 sigma.)

These are the most dark-matter-dominated galaxies known, and where MOND's
"external field effect" bites: the host's pull can suppress the MOND
boost. Each dwarf carries the host's acceleration at its distance,
g_ext = V_host^2 / D, with flat host circular speeds of 190 km/s (Milky
Way), 230 km/s (Andromeda) and 67 km/s (NGC 3109, for Antlia B);
field dwarfs use their nearest giant.

Uncertainty: the measured dispersion error (upper or lower side, as
appropriate) plus the effect of a 0.2 dex stellar-mass uncertainty,
propagated through each model's own prediction. Upper limits are one-sided:
a prediction at or below the limit passes; above it, z = 2 x pred / limit.

Caveats: single-epoch dispersions of the faintest systems may be inflated
by binary stars (upward only). Gas-rich field dwarfs (WLM, Leo A, IC 1613,
Pegasus, Sag dIrr) also rotate somewhat; the test uses their dispersions.

Model method: predict_dwarf_dispersion(d, gp) -> line-of-sight velocity
dispersion [km/s], given
    d.M_bar   stellar mass [Msun] (as before: stars only)
    d.M_gas   HI mass x 1.33 for helium [Msun], 0 if none detected; mostly
              outside the half-light radius, so a model decides how much
              to include
    d.r_half  3D half-light radius [kpc]
    d.D_host  distance to host / nearest giant [kpc]
    d.g_ext   host's acceleration [(km/s)^2/kpc]
    d.host    "MW", "M31" or "field"
"""
import copy
import csv
import os

import numpy as np

from ..core import Test, data_path, is_hidden

V_HOST = {"MW": 190.0, "M31": 230.0, "ngc_3109": 67.0}
LOG_MSTAR_ERR = 0.2
FILE = "lvdb_dwarfs.csv"


def _f(x):
    return float(x) if x not in ("", None) else None


class Dwarf:
    def __init__(self, r):
        self.name, self.host = r["name"], r["host"]
        self.nearest_giant = r["nearest_giant"]
        self.D_host = float(r["D_host"])
        self.M_bar = 10 ** float(r["log_Mstar"])
        self.M_gas = 1.33 * 10 ** float(r["log_MHI"]) if r["log_MHI"] else 0.0
        self.log_Mb_err = LOG_MSTAR_ERR
        self.r_half = float(r["r_half_pc"]) / 1000.0
        self.sigma_obs = _f(r["sigma"])
        self.sigma_em, self.sigma_ep = _f(r["sigma_em"]), _f(r["sigma_ep"])
        self.sigma_ul = _f(r["sigma_ul"])
        self.is_limit = self.sigma_obs is None
        # symmetric error kept for models/tools that expect one number
        self.sigma_err = None if self.is_limit else 0.5 * (self.sigma_em + self.sigma_ep)
        self.g_ext = V_HOST[self.nearest_giant] ** 2 / self.D_host


def load_dwarfs():
    path = data_path("dwarfs", FILE)
    return [Dwarf(r) for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]


class DwarfSpheroidals(Test):
    ID = "05_dwarfs"
    LABEL = "Dwarf galaxies"

    def load(self):
        self.dwarfs = load_dwarfs()

    def n_bound(self):
        # upper limits pass automatically for a perfect theory; about half
        # of them land in the hidden part
        return int(round(sum(d.is_limit for d in self.dwarfs) / 2))

    def export_public(self, dest_data, salt):
        os.makedirs(os.path.join(dest_data, "dwarfs"), exist_ok=True)
        with open(data_path("dwarfs", FILE)) as fi, \
                open(os.path.join(dest_data, "dwarfs", FILE), "w") as fo:
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
            if d.is_limit:
                z.append(0.0 if pred <= d.sigma_ul else 2.0 * pred / d.sigma_ul)
                continue
            d_up = copy.copy(d)
            d_up.M_bar = d.M_bar * 10 ** d.log_Mb_err
            pred_up = model.predict_dwarf_dispersion(d_up, gp)
            e_obs = d.sigma_ep if pred > d.sigma_obs else d.sigma_em
            sig = np.hypot(e_obs, pred_up - pred)
            z.append((pred - d.sigma_obs) / sig)
        return np.array(z)
