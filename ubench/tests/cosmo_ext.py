"""
Stage 3b: extended cosmology (tests 28-31, 34 and open question 29b).

Data: data/extended/constraints.csv (each row checked against its paper's
abstract) and data/growth/ (SDSS DR16 growth-rate files from Cobaya).
Predictions come from model.cosmology(gp), plus predict_gravity for the
number of spacetime dimensions.

28 Inflation: primordial gravitational waves (r) and non-Gaussianity (f_NL).
29 Neutrinos: number of relativistic species, KATRIN's upper limit and the
   oscillation lower limit on the sum of neutrino masses.
29b (CONTESTED): the cosmological upper limits on the neutrino mass sum,
   which now press against the oscillation minimum.
30 Curvature and topology: Omega_K (Planck + BAO) and the minimum size of a
   finite (torus) universe.
31 Growth rate: f sigma8 from SDSS/eBOSS redshift-space distortions.
34 Standard sirens: GW170817's Hubble constant and the number of spacetime
   dimensions (gravity leaking into extra dimensions would weaken the
   gravitational waves relative to light).

Bounds score 0 inside the allowed region and grow beyond it (z = 2 at the
bound times 1 + the fractional excess).
"""
import csv

import numpy as np

from ..core import Test, data_path
from .cosmo import _asym_z, _cosmo, _rows


def bound_z(pred, row):
    v = float(row["value"])
    if row["kind"] == "upper":
        return 0.0 if pred <= v else 2.0 * pred / v
    if row["kind"] == "lower":
        return 0.0 if pred >= v else 2.0 * v / max(pred, 1e-12)
    return _asym_z(pred, v, float(row["err_plus"]), float(row["err_minus"]))


def _get(c, model, gp, q):
    """Prediction for quantity q from the cosmology object (or gravity)."""
    if q == "tensor_to_scalar_r":
        return c.r
    if q == "f_NL_local":
        return c.f_NL
    if q == "N_eff":
        return c.Neff
    if q == "sum_mnu_eV":
        return c.mnu
    if q == "Omega_K":
        return c.omk
    if q == "topology_Ri_over_chi_rec":
        return c.topology_Ri_over_chi_rec
    if q == "H0_siren":
        return c.H0
    if q == "spacetime_dimensions":
        return model.predict_gravity("spacetime_dimensions", gp) if hasattr(model, "predict_gravity") else None
    raise KeyError(q)


class _ExtTable(Test):
    STAGE = 3
    GUARD_RAIL = True    # standard cosmologies pass these automatically; they catch exotic physics
    HIDDEN_SPLIT = False
    DATA_DIRS = ("extended",)
    NUMBER = ""

    def load(self):
        self.rows = [r for r in _rows("extended", "constraints.csv") if r["test"] == self.NUMBER]

    def n_bound(self):
        return sum(r["kind"] != "gauss" for r in self.rows)

    def residuals(self, model, gp, part, salt):
        c = _cosmo(model, gp)
        if c is None:
            return None
        if c == "bad":
            return np.full(len(self.rows), 1e6)
        z = []
        for r in self.rows:
            v = _get(c, model, gp, r["quantity"])
            if v is None:
                return None
            z.append(bound_z(v, r))
        return np.array(z)


class Inflation(_ExtTable):
    ID, LABEL, NUMBER = "28_inflation", "Inflation signatures", "28"


class Neutrinos(_ExtTable):
    ID, LABEL, NUMBER = "29_neutrinos", "Neutrino masses and species", "29"


class NeutrinoSqueeze(_ExtTable):
    ID, LABEL, NUMBER = "29b_neutrino_squeeze", "Neutrino mass squeeze", "29b"
    GUARD_RAIL = False
    STATUS = "contested"
    USED_IN_FIT = False


class Curvature(_ExtTable):
    ID, LABEL, NUMBER = "30_curvature", "Curvature and topology", "30"


class StandardSirens(_ExtTable):
    ID, LABEL, NUMBER = "34_standard_sirens", "Standard sirens", "34"


class GrowthRate(Test):
    ID, LABEL = "31_growth_rate", "Growth rate of structure"
    STAGE = 3
    HIDDEN_SPLIT = False
    DATA_DIRS = ("growth",)
    FILES = (("sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8.dat", "sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8_covtot.txt"),
             ("sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8.dat", "sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8_covtot.txt"))

    def load(self):
        self.z, self.val, self.Linv = [], [], []
        for mean, cov in self.FILES:
            rows = [l.split() for l in open(data_path("growth", mean)) if l.strip() and not l.startswith("#")]
            C = np.loadtxt(data_path("growth", cov))
            idx = [i for i, r in enumerate(rows) if r[2].startswith("f_sigma8")]
            sub = C[np.ix_(idx, idx)]
            self.z.append(np.array([float(rows[i][0]) for i in idx]))
            self.val.append(np.array([float(rows[i][1]) for i in idx]))
            self.Linv.append(np.linalg.inv(np.linalg.cholesky(sub)))

    def residuals(self, model, gp, part, salt):
        c = _cosmo(model, gp)
        if c is None:
            return None
        if c == "bad":
            return np.full(sum(len(v) for v in self.val), 1e6)
        out = []
        for z, v, L in zip(self.z, self.val, self.Linv):
            pred = np.array([c.fsigma8(zz) for zz in z])
            out.append(L @ (pred - v))
        return np.concatenate(out)
