"""Candidate 9: cosmologically coupled black holes as dark energy
(Farrah et al. 2023; Croker et al. 2024, arXiv:2405.12282).

With coupling strength k = 3, a black hole's mass grows as a^3 while black
holes dilute as a^-3, so each generation contributes a CONSTANT physical
density from its birth on, a cosmological constant switched on gradually as
stars die. The dark-energy density is then
    rho_DE(a) = f * Int_0^a psi(a') / a'^3 dt'
with psi the cosmic star-formation rate density (Madau & Dickinson 2014,
ARA&A 52, 415: psi = 0.015 (1+z)^2.7 / (1 + ((1+z)/2.9)^5.6) Msun/yr/Mpc^3)
and f fixed by today's dark-energy density: the SHAPE of the dark-energy
history has no free constant. No cosmological constant. The black holes'
birth mass comes out of the baryons (Croker et al. 2024, eq. 2.2), so the
matter density falls at late times: with this star-formation history 46% of
the baryons (Omega = 0.023) have been converted by today. omch2/ombh2 are
the early-universe values.
* Fast engine: exact rho_DE(a) shape (tabulated), growth ODE with it.
* CLASS (test 15): the CPL (w0, wa) whose H(z) matches the exact history
  to < 0.3% for z < 4 (dark energy is negligible at recombination, so the
  CMB mostly sees the distance to last scattering).
Galaxies, clusters: the LCDM baseline. Action: not checkable (the coupling
is a statement about black-hole interiors, outside Horndeski).
"""
import numpy as np
import pathlib

from scipy.integrate import cumulative_trapezoid, solve_ivp

from ubench.cosmology import Cosmology
from ubench.models.lcdm import LCDM


_ROOT = pathlib.Path(__file__).resolve().parents[1]


def _shape():
    """rho_DE(a)/rho_DE(1), the comoving mass that has gone into black holes
    M(a)/M(1), and the ratio R = (baryons converted today)/(rho_DE today),
    from the star-formation history (H from LCDM for dt; the shapes are
    insensitive to that choice)."""
    g0, g1, g2, g3 = np.loadtxt(_ROOT / "data/star_formation/madau_dickinson2014.csv",
                                delimiter=",", comments="#", skiprows=6)
    z = np.linspace(0, 20, 4001)
    psi = g0 * (1 + z) ** g1 / (1 + ((1 + z) / g2) ** g3)
    Hz = np.sqrt(0.315 * (1 + z) ** 3 + 0.685)
    dt_dz = 1 / ((1 + z) * Hz)
    rev = lambda f: -cumulative_trapezoid(f[::-1], z[::-1], initial=0)[::-1]   # Int_z^inf
    de = rev(psi * (1 + z) ** 3 * dt_dz)                  # psi / a^3 dt: constant after birth
    mb = rev(psi * dt_dz)                                 # comoving mass converted
    a = 1 / (1 + z)
    return a[::-1], (de / de[0])[::-1], (mb / mb[0])[::-1], mb[0] / de[0]


A_TAB, S_TAB, M_TAB, R_MISSING = _shape()


def de_shape(a):
    return np.interp(a, A_TAB, S_TAB, left=0.0, right=1.0)


def missing_shape(a):
    """Converted baryon density / rho_DE(today), physical, at scale factor a."""
    return R_MISSING * np.interp(a, A_TAB, M_TAB, left=0.0, right=1.0) * np.asarray(a, float) ** -3


def net_shape(a):
    """What a matter-conserving analysis sees as dark energy: the black-hole
    density minus the baryons converted into them. Normalized to 1 today:
    the black holes carry rho_DE/(1 - R) and the converted baryons
    R rho_DE/(1 - R), with Omega_DE the net value."""
    return (de_shape(a) - missing_shape(a)) / (1 - R_MISSING)


def cpl_fit():
    """CPL (w0, wa) whose H(z) best matches the exact history for z < 4
    (the naive fit to rho_DE shape is poor: the true w(z) is phantom at all
    z, -1.00 today, -1.08 at z = 1, -1.30 at z = 2)."""
    from scipy.optimize import least_squares
    a = np.linspace(0.2, 1, 400)
    E_true = np.sqrt(0.315 * a ** -3 + 0.685 * net_shape(a))
    def res(p):
        cpl = a ** (-3 * (1 + p[0] + p[1])) * np.exp(-3 * p[1] * (1 - a))
        return np.sqrt(0.315 * a ** -3 + 0.685 * cpl) / E_true - 1
    return least_squares(res, [-1.0, 0.0]).x


W0, WA = cpl_fit()


class CCBHCosmology(Cosmology):
    def E(self, z):
        z = np.asarray(z, dtype=float)
        a = 1.0 / (1.0 + z)
        base = super().E(z) ** 2 - self.ode * a ** (-3 * (1 + self.w0 + self.wa)) * np.exp(-3 * self.wa * (1 - a))
        return np.sqrt(base + self.ode * net_shape(a))

    def growth(self, z):
        if not hasattr(self, "_growth_a"):
            om = self.om
            Eg = lambda a: np.sqrt(om * a ** -3 + (1 - om) * net_shape(a))

            def rhs(lna, y):
                a = np.exp(lna)
                E = Eg(a)
                dlnE = (np.log(Eg(a * np.exp(1e-5))) - np.log(E)) / 1e-5
                return [y[1], -(2 + dlnE) * y[1] + 1.5 * (om * a ** -3 - (1 - om) * missing_shape(a) / (1 - R_MISSING)) / E ** 2 * y[0]]
            self._growth_a = solve_ivp(rhs, [np.log(1e-3), 0.0], [1e-3, 1e-3], dense_output=True,
                                       rtol=1e-8, atol=1e-12)
        lna = np.log(1.0 / (1.0 + np.asarray(z, dtype=float)))
        return self._growth_a.sol(lna)[0]


class CoupledBH(LCDM):
    NAME = "Cosmologically coupled black holes (k = 3)"
    ACTION = {"type": "other", "note": "black-hole interior coupling, outside Horndeski"}

    def class_params(self, gp):
        p = super().class_params(gp)
        if p is not None:
            p.update({"Omega_Lambda": 0.0, "w0_fld": float(W0), "wa_fld": float(WA), "use_ppf": "yes"})
        return p

    def cosmology(self, gp):
        keys = ("H0", "ombh2", "omch2", "ns", "ln10As", "mnu")
        args = {k: gp[k] for k in keys if k in gp}
        key = tuple(sorted(args.items()))
        cache = self.__dict__.setdefault("_cosmo_cache", {})
        if key not in cache:
            if len(cache) > 8:
                cache.clear()
            cache[key] = CCBHCosmology(**args)
        return cache[key]


MODEL = CoupledBH()
