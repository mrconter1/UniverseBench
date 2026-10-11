"""
Open questions in cosmology (CONTESTED; reported, never scored or fitted).

32 Cosmological principle: the dipole in the number counts of distant
   quasars is about twice what our motion relative to the CMB implies
   (Secrest+ 2021, 4.9 sigma). Any FLRW (homogeneous, isotropic) model
   predicts the kinematic value ~0.0072; models may override with
   predict_number_count_dipole(gp).
33 Dark flow: Planck finds no large-scale bulk flow (< 254 km/s at 95% on
   Gpc scales), against earlier claims of a large coherent flow. Models
   give the rms linear-theory bulk velocity in a 2 Gpc sphere.
35 Early supermassive black holes: a 1.6e9 Msun black hole only ~0.7 Gyr
   after the Big Bang (J0313-1806, z = 7.642). Yardstick: can it grow from
   a 100 Msun seed formed at z = 30, accreting at the Eddington limit
   (e-folding time 50 Myr)? Depends on the model's cosmic clock (and on
   astrophysics; heavy seeds or super-Eddington growth would ease it).
   Models may override with predict_max_bh_mass(z, gp).
"""
import numpy as np

from ..core import Test
from .cosmo import _cosmo, _rows

KINEMATIC_DIPOLE = 0.0072
SEED_MSUN, SEED_Z, EFOLD_GYR = 100.0, 30.0, 0.050


def bulk_flow_rms(c, R_mpc):
    """rms bulk velocity [km/s] in a top-hat sphere of radius R (linear theory, z = 0)."""
    lnk = np.linspace(np.log(1e-5), np.log(10.0), 3000)
    k = np.exp(lnk)
    x = k * R_mpc
    W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
    integ = c.delta2_linear(k, 0.0) / k**2 * W**2           # Delta^2/k^2 W^2 dlnk
    return float(c.growth_rate(0.0) * c.H0 * np.sqrt(np.trapezoid(integ, lnk)))


class _Open(Test):
    STAGE = 3
    STATUS = "contested"
    HIDDEN_SPLIT = False
    USED_IN_FIT = False
    DATA_DIRS = ("open_questions",)
    NUMBER = ""

    def load(self):
        self.row = [r for r in _rows("open_questions", "cosmic.csv") if r["test"] == self.NUMBER][0]


class CosmicDipole(_Open):
    ID, LABEL, NUMBER = "32_cosmic_dipole", "Cosmic dipole (isotropy)", "32"

    def residuals(self, model, gp, part, salt):
        if hasattr(model, "predict_number_count_dipole"):
            d = model.predict_number_count_dipole(gp)
        else:
            c = _cosmo(model, gp)
            if c is None:
                return None
            d = KINEMATIC_DIPOLE
        return np.array([(d - float(self.row["value"])) / float(self.row["err"])])


class DarkFlow(_Open):
    ID, LABEL, NUMBER = "33_dark_flow", "Dark flow (bulk motions)", "33"
    N_BOUND = 1

    def residuals(self, model, gp, part, salt):
        c = _cosmo(model, gp)
        if c is None:
            return None
        if c == "bad":
            return np.array([1e6])
        v = bulk_flow_rms(c, 2000.0)
        lim = float(self.row["value"])
        return np.array([0.0 if v <= lim else 2.0 * v / lim])


class EarlyBlackHoles(_Open):
    ID, LABEL, NUMBER = "35_early_black_holes", "Early supermassive black holes", "35"
    N_BOUND = 1

    def residuals(self, model, gp, part, salt):
        z_q = 7.642
        if hasattr(model, "predict_max_bh_mass"):
            m_max = model.predict_max_bh_mass(z_q, gp)
        else:
            c = _cosmo(model, gp)
            if c is None:
                return None
            if c == "bad":
                return np.array([1e6])
            dt = c.age_gyr(z_q) - c.age_gyr(SEED_Z)
            m_max = SEED_MSUN * np.exp(dt / EFOLD_GYR)
        need = float(self.row["value"])
        return np.array([max(0.0, (need - np.log10(m_max)) / float(self.row["err"]))])
