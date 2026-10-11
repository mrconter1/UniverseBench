"""
Test D: the Milky Way.

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
import os

import numpy as np
from scipy.special import ellipk, i0, i1, k0, k1

HERE = os.path.dirname(os.path.abspath(__file__))
G = 4.30091e-6        # kpc (km/s)^2 / Msun
FLOOR = 0.03

# McMillan 2017, Table 3 (best fit) and Table 1 (gas)
BULGE_MASS = 9.23e9
THIN = dict(S0=896e6, Rd=2.50)            # Msun/kpc^2, kpc
THICK = dict(S0=183e6, Rd=3.02)
HI = dict(S0=53.1e6, Rd=7.0, Rm=4.0)
H2 = dict(S0=2180e6, Rd=1.5, Rm=12.0)
STELLAR_MASS = 5.43e10


def v2_exponential(R, S0, Rd):
    """Exact speed^2 of a razor-thin exponential disc (Freeman 1970)."""
    y = R / (2 * Rd)
    return 4 * np.pi * G * S0 * Rd * y**2 * (i0(y) * k0(y) - i1(y) * k1(y))


def v2_numeric(R, sigma, a_max=80.0, n=40000):
    """Speed^2 of any razor-thin axisymmetric disc with surface density
    sigma(a), by summing the potentials of thin rings and differentiating."""
    a = (np.arange(n) + 0.5) * (a_max / n)
    dm = 2 * np.pi * a * sigma(a) * (a_max / n)

    def phi(r):
        m = 4 * r * a / (r + a) ** 2
        ok = m < 1 - 1e-12        # skip the one ring sitting exactly at r
        return -np.sum(2 * G * dm[ok] * ellipk(m[ok]) / (np.pi * (r + a[ok])))

    out = []
    for r in np.atleast_1d(R):
        h = 1.0137e-3 * r          # odd step so r +/- h avoids ring radii
        out.append(r * (phi(r + h) - phi(r - h)) / (2 * h))
    return np.array(out)


def gas_sigma(p):
    return lambda a: p["S0"] * np.exp(-p["Rm"] / a - a / p["Rd"])


class MilkyWay:
    name = "MilkyWay"
    ups_fixed = True

    def __init__(self, path=None):
        path = path or os.path.join(HERE, "data", "milkyway", "eilers2019.csv")
        lines = [l for l in open(path) if l[:1].isdigit()]
        d = np.loadtxt(lines, delimiter=",")
        self.R, self.Vobs = d[:, 0], d[:, 1]
        self.errV = (d[:, 2] + d[:, 3]) / 2
        self.sigma = np.sqrt(self.errV**2 + (FLOOR * self.Vobs) ** 2)
        v2_disk = v2_exponential(self.R, **THIN) + v2_exponential(self.R, **THICK)
        v2_gas = v2_numeric(self.R, gas_sigma(HI)) + v2_numeric(self.R, gas_sigma(H2))
        v2_bul = G * BULGE_MASS / self.R          # all bulge mass lies inside 5 kpc
        self.Vdisk, self.Vgas, self.Vbul = (np.sqrt(v2_disk), np.sqrt(v2_gas),
                                            np.sqrt(v2_bul))

    def stellar_mass(self, ups=None):
        return STELLAR_MASS


def residuals(model, mw, gp):
    lp = {k: v[0] for k, v in model.LOCAL_PARAMS.items()}
    return (model.predict(mw, mw.R, gp, lp) - mw.Vobs) / mw.sigma
