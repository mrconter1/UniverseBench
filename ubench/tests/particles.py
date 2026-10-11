"""
Stage 4 particle guard rails (tests 24, 25, 36).

Many dark matter proposals come with a particle that should show up in
experiments. These tests check a model's declared particle, from
model.dm_particle(gp), against null results:

24 Direct detection: LZ's limit on spin-independent scattering off nuclei.
25 Indirect detection: Fermi-LAT's dwarf-galaxy gamma-ray limits exclude a
   thermal-relic annihilation rate for dark matter lighter than ~100 GeV.
36 Cosmic antimatter: no antihelium in cosmic rays (BESS-Polar), ruling
   out large primordial antimatter domains. Models may predict a ratio with
   predict_antihelium_ratio(gp); otherwise standard baryogenesis (no
   antimatter domains, ratio 0) is assumed.

dm_particle returns a dict (mass_GeV, sigma_SI_cm2, sigma_v_cm3s) or None
for "no dark matter particle" (e.g. MOND), which predicts no signal and
passes. A model without the method makes no prediction ("—").
"""
import numpy as np

from ..core import Test
from .cosmo import _rows


class _Particle(Test):
    STAGE = 4
    HIDDEN_SPLIT = False
    DATA_DIRS = ("particles",)
    NUMBER = ""

    def load(self):
        self.row = [r for r in _rows("particles", "limits.csv") if r["test"] == self.NUMBER][0]


class DirectDetection(_Particle):
    ID, LABEL, NUMBER = "24_direct_detection", "Direct detection", "24"

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "dm_particle"):
            return None
        p = model.dm_particle(gp)
        if not p or not p.get("sigma_SI_cm2"):
            return np.array([0.0])
        m = p.get("mass_GeV") or float(self.row["ref_mass_GeV"])
        lim = float(self.row["limit"]) * max(1.0, m / float(self.row["ref_mass_GeV"]))
        return np.array([0.0 if p["sigma_SI_cm2"] <= lim else 2.0 * p["sigma_SI_cm2"] / lim])


class IndirectDetection(_Particle):
    ID, LABEL, NUMBER = "25_indirect_detection", "Indirect detection", "25"

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "dm_particle"):
            return None
        p = model.dm_particle(gp)
        if not p or not p.get("sigma_v_cm3s") or (p.get("mass_GeV") or 1e9) > float(self.row["ref_mass_GeV"]):
            return np.array([0.0])
        lim = float(self.row["limit"])
        return np.array([0.0 if p["sigma_v_cm3s"] < lim else 2.0 * p["sigma_v_cm3s"] / lim])


class CosmicAntimatter(_Particle):
    ID, LABEL, NUMBER = "36_cosmic_antimatter", "Cosmic antimatter", "36"

    def residuals(self, model, gp, part, salt):
        r = model.predict_antihelium_ratio(gp) if hasattr(model, "predict_antihelium_ratio") else 0.0
        lim = float(self.row["limit"])
        return np.array([0.0 if r <= lim else 2.0 * r / lim])
