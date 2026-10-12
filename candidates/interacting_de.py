"""Candidate B, dark matter interacting with dark energy (the phenomenological
"interacting dark energy" form of coupled quintessence).

Energy flows between the two dark components at a rate Q = xi * H * rho_DE
(xi > 0: dark energy decays into dark matter). With w = -1 the background is
exact:
    rho_DE(a) = rho_DE0 a^-xi
    rho_DM(a) = rho_DM,early a^-3 + xi/(3 - xi) rho_DE0 a^-xi
so the dark-matter density the early universe sees (set by omch2, which
fixes the sound horizon and recombination) differs from today's. xi = 0 is
exactly LCDM.

Growth: dark matter feels the transfer as a drag, D'' + (2 + dlnE/dlna + q) D'
= 1.5 Omega_m(a) D, with q = Q/(H rho_DM). Dark-energy perturbations and the
q' source term are neglected (approximation). The compressed CMB shift
parameter uses the early matter density.

Galaxies, clusters and everything below cosmological scales: the LCDM
baseline (an interaction this weak is irrelevant inside halos).
"""
import numpy as np
from scipy.integrate import solve_ivp

from ubench.cosmology import Cosmology, _nu_ratio, T_NU_EV, C_KMS
from ubench.models.lcdm import LCDM


class IDECosmology(Cosmology):
    def __init__(self, xi=0.0, **kw):
        self.xi, self.dm_extra = 0.0, 0.0          # parent builds a LCDM table first
        super().__init__(**kw)
        self.xi = xi
        self.om_early = self.om                          # a^-3 part today (early-universe density)
        leftover = self.ode                              # what LCDM would assign to dark energy
        self.ode = leftover * (3 - xi) / 3               # flatness with the transferred matter
        self.dm_extra = xi / (3 - xi) * self.ode         # extra dark matter today
        self.om = self.om_early + self.dm_extra          # total matter today
        self._build_distance_table()

    def E(self, z):
        z = np.asarray(z, dtype=float)
        a = 1.0 / (1.0 + z)
        rad = (self.og + self.onu_massless) * a ** -4
        nu = self.onu_rel_each * a ** -4 * (_nu_ratio(self.mnu * a / T_NU_EV) if self.mnu > 0 else 1.0)
        dark = (self.ode + self.dm_extra) * a ** -self.xi
        return np.sqrt(rad + nu + self.ocb * a ** -3 + self.omk * a ** -2 + dark)

    def cmb_distance_priors(self):
        # R must use the matter density at recombination, not today's
        om_today, self.om = self.om, self.om_early
        try:
            return super().cmb_distance_priors()
        finally:
            self.om = om_today

    def growth(self, z):
        if not hasattr(self, "_growth_a"):
            ob = self.ocb - self.oc
            oc_e, ode, x = self.oc + self.omnuh2 / self.h ** 2, self.ode, self.xi

            def parts(a):
                rde = ode * a ** -x
                rdm = oc_e * a ** -3 + x / (3 - x) * ode * a ** -x
                return rde, rdm

            def Eg(a):
                rde, rdm = parts(a)
                return np.sqrt(ob * a ** -3 + rdm + rde + self.omk * a ** -2 - self.omk)

            def rhs(lna, y):
                a = np.exp(lna)
                E = Eg(a)
                eps = 1e-5
                dlnE = (np.log(Eg(a * np.exp(eps))) - np.log(E)) / eps
                rde, rdm = parts(a)
                om_a = (ob * a ** -3 + rdm) / E ** 2
                q = x * rde / rdm
                D, dD = y
                return [dD, -(2 + dlnE + q) * dD + 1.5 * om_a * D]
            a0 = 1e-3
            self._growth_a = solve_ivp(rhs, [np.log(a0), 0.0], [a0, a0], dense_output=True,
                                       rtol=1e-8, atol=1e-12)
        lna = np.log(1.0 / (1.0 + np.asarray(z, dtype=float)))
        return self._growth_a.sol(lna)[0]


class InteractingDE(LCDM):
    NAME = "Interacting dark energy (dark matter <-> dark energy)"
    GLOBAL_PARAMS = {**LCDM.GLOBAL_PARAMS, "xi": (0.0, -0.6, 0.6)}

    def class_params(self, gp):
        """Vanilla CLASS cannot express a Q = xi H rho_DE energy exchange, so
        this model makes no prediction for test 15 (full CMB spectra) until a
        modified CLASS (an interacting-dark-energy fork) is wired in. Returning
        the LCDM spectra instead would silently score it as LCDM."""
        return None

    def cosmology(self, gp):
        keys = ("H0", "ombh2", "omch2", "ns", "ln10As", "mnu", "xi")
        args = {k: gp[k] for k in keys if k in gp}
        key = tuple(sorted(args.items()))
        cache = self.__dict__.setdefault("_cosmo_cache", {})
        if key not in cache:
            if len(cache) > 8:
                cache.clear()
            cache[key] = IDECosmology(**args)
        return cache[key]


MODEL = InteractingDE()
