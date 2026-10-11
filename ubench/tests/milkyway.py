"""
Test 04: the Milky Way.

Data: circular velocity of our galaxy from 5 to 25 kpc, measured from
~23,000 red giant stars with Gaia, APOGEE, WISE and 2MASS (Eilers et al.
2019). Quoted errors are statistical; the paper estimates 2-5% systematic
uncertainty, so a 3% floor is added.

Visible matter: the McMillan (2017) Milky Way mass model: bulge, thin and
thick stellar discs, atomic (HI) and molecular (H2) gas discs. Caveat: its
stellar disc masses were fitted alongside a dark matter halo, so they are
not entirely independent of the dark matter hypothesis.

This is a pure out-of-sample test. The Milky Way is never used to fit
constants, and no per-galaxy tuning is allowed: a model must predict our
galaxy with the constants it learned from the other tests.

The Milky Way is handed to models as a galaxy object, just like SPARC
galaxies, except that ups_fixed = True: its stellar speeds already
correspond to real masses, so mass-to-light ratios must not be applied.
"""
import numpy as np

from ..core import Test, data_path, default_locals
from ..physics import G, v2_exponential_disc, v2_thin_disc

FLOOR = 0.03

# McMillan 2017, Table 3 (best fit) and Table 1 (gas)
BULGE_MASS = 9.23e9
THIN = dict(S0=896e6, Rd=2.50)            # Msun/kpc^2, kpc
THICK = dict(S0=183e6, Rd=3.02)
HI = dict(S0=53.1e6, Rd=7.0, Rm=4.0)
H2 = dict(S0=2180e6, Rd=1.5, Rm=12.0)
STELLAR_MASS = 5.43e10


def gas_sigma(p):
    return lambda a: p["S0"] * np.exp(-p["Rm"] / a - a / p["Rd"])


class MilkyWay:
    name = "MilkyWay"
    ups_fixed = True

    def __init__(self, path=None):
        path = path or data_path("milkyway", "eilers2019.csv")
        lines = [l for l in open(path) if l[:1].isdigit()]
        d = np.loadtxt(lines, delimiter=",")
        self.R, self.Vobs = d[:, 0], d[:, 1]
        self.errV = (d[:, 2] + d[:, 3]) / 2
        self.sigma = np.sqrt(self.errV**2 + (FLOOR * self.Vobs) ** 2)
        v2_disk = v2_exponential_disc(self.R, **THIN) + v2_exponential_disc(self.R, **THICK)
        v2_gas = v2_thin_disc(self.R, gas_sigma(HI)) + v2_thin_disc(self.R, gas_sigma(H2))
        v2_bul = G * BULGE_MASS / self.R          # all bulge mass lies inside 5 kpc
        self.Vdisk, self.Vgas, self.Vbul = (np.sqrt(v2_disk), np.sqrt(v2_gas),
                                            np.sqrt(v2_bul))

    def stellar_mass(self, ups=None):
        return STELLAR_MASS


class MilkyWayTest(Test):
    ID = "04_milkyway"
    LABEL = "Milky Way"
    HIDDEN_SPLIT = False
    USED_IN_FIT = False     # pure out-of-sample
    DATA_DIRS = ("milkyway",)

    def load(self):
        self.mw = MilkyWay()

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict"):
            return None
        mw = self.mw
        return (np.asarray(model.predict(mw, mw.R, gp, default_locals(model))) - mw.Vobs) / mw.sigma
