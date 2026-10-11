"""
A small, self-contained cosmology engine for the Stage 3 tests.

Flat Friedmann universe with photons, 3 neutrino species (one massive, sum
of masses mnu, total N_eff), baryons, cold dark matter and dark energy with
a CPL equation of state w(a) = w0 + wa (1 - a). Provides:

  * expansion history E(z) = H(z)/H0 and distances (D_M, D_H, D_V, D_L)
  * recombination z* (Hu & Sugiyama 1996) and the sound horizon r_s(z*)
  * the BAO drag-epoch sound horizon r_d (Aubourg+ 2015 fitting formula,
    0.02% accurate for standard physics)
  * CMB distance priors R and l_A (as defined by Chen, Huang & Wang 2019)
  * linear growth D(z) and f(z), the linear matter power spectrum from the
    Eisenstein & Hu (1998) no-wiggle transfer function, sigma8 and S8
  * the Sheth-Tormen halo mass function
  * Big Bang nucleosynthesis yields from the Steigman fitting formulas

Validated against Planck 2018 best-fit values (see validate()).

Units: distances in Mpc, H in km/s/Mpc, masses in Msun.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp

C_KMS = 299792.458
T_CMB = 2.7255
OMEGA_GAMMA_H2 = 2.47282e-5 * (T_CMB / 2.7255) ** 4
NU_MASS_CONV = 93.14          # sum m_nu [eV] / Omega_nu h^2 for non-relativistic neutrinos
RHO_CRIT_H2 = 2.775e11        # critical density today, Msun / Mpc^3 / h^2
Z_STAR_CALIB = 1089.92 / 1091.9081   # see Cosmology.z_star

# Massive-neutrino energy density relative to its massless value, as a
# function of y = m_nu c^2 / (k T_nu): F(y) = int q^2 sqrt(q^2+y^2)/(e^q+1) / int q^3/(e^q+1)
_Q = np.linspace(1e-4, 40, 4000)
_DEN = np.trapezoid(_Q**3 / (np.exp(_Q) + 1), _Q)
_Y = np.logspace(-4, 5, 600)
_F = np.array([np.trapezoid(_Q**2 * np.sqrt(_Q**2 + y**2) / (np.exp(_Q) + 1), _Q) for y in _Y]) / _DEN
T_NU_EV = 8.617333e-5 * T_CMB * (4 / 11) ** (1 / 3)


def _nu_ratio(y):
    y = np.asarray(y, dtype=float)
    return np.where(y < _Y[0], 1.0, np.where(y > _Y[-1], 0.180309 * y,
                    np.exp(np.interp(np.log(np.clip(y, _Y[0], _Y[-1])), np.log(_Y), np.log(_F)))))


class Cosmology:
    def __init__(self, H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, Neff=3.046,
                 ns=0.9649, ln10As=3.044, w0=-1.0, wa=0.0):
        self.H0, self.h = H0, H0 / 100.0
        self.ombh2, self.omch2, self.mnu, self.Neff = ombh2, omch2, mnu, Neff
        self.ns, self.ln10As, self.w0, self.wa = ns, ln10As, w0, wa
        self.As = np.exp(ln10As) * 1e-10
        h2 = self.h**2
        self.og = OMEGA_GAMMA_H2 / h2
        # one massive species carrying all of mnu, plus (Neff - Neff/3) massless
        self.onu_rel_each = (7 / 8) * (4 / 11) ** (4 / 3) * self.og * Neff / 3
        self.onu_massless = 2 * self.onu_rel_each
        self.omnuh2 = mnu / NU_MASS_CONV
        self.ob, self.oc = ombh2 / h2, omch2 / h2
        self.ocb = self.ob + self.oc
        self.om = self.ocb + self.omnuh2 / h2        # matter today (neutrinos non-relativistic)
        onu_massive_today = self.onu_rel_each * float(_nu_ratio(mnu / T_NU_EV)) if mnu > 0 \
            else self.onu_rel_each
        self.ode = 1.0 - self.ocb - self.og - self.onu_massless - onu_massive_today
        self._build_distance_table()

    # ------------------------------------------------------------ background
    def E(self, z):
        z = np.asarray(z, dtype=float)
        a = 1.0 / (1.0 + z)
        rad = (self.og + self.onu_massless) * a**-4
        if self.mnu > 0:
            nu = self.onu_rel_each * a**-4 * _nu_ratio(self.mnu * a / T_NU_EV)
        else:
            nu = self.onu_rel_each * a**-4
        de = self.ode * a ** (-3 * (1 + self.w0 + self.wa)) * np.exp(-3 * self.wa * (1 - a))
        return np.sqrt(rad + nu + self.ocb * a**-3 + de)

    def H(self, z):
        return self.H0 * self.E(z)

    def _build_distance_table(self):
        z = np.concatenate([np.linspace(0, 3, 3001)[:-1], np.logspace(np.log10(3), np.log10(1500), 3000)])
        integrand = 1.0 / self.E(z)
        self._z = z
        self._dc = C_KMS / self.H0 * cumulative_trapezoid(integrand, z, initial=0.0)

    def D_M(self, z):
        return np.interp(z, self._z, self._dc)

    def D_H(self, z):
        return C_KMS / self.H(z)

    def D_V(self, z):
        return (z * self.D_M(z) ** 2 * self.D_H(z)) ** (1 / 3)

    def D_L(self, z, z_hel=None):
        return (1 + (z if z_hel is None else z_hel)) * self.D_M(z)

    def distance_modulus(self, z, z_hel=None):
        return 5 * np.log10(self.D_L(z, z_hel)) + 25

    # ------------------------------------------------------------ sound horizon / CMB
    def z_star(self):
        """Recombination redshift: Hu & Sugiyama (1996, Eq. E1), rescaled by
        Z_STAR_CALIB so that the Planck 2018 best fit gives Planck's own
        z* = 1089.92 (the raw formula gives 1091.91 there). The acoustic
        angle theta* is measured to 0.03%, which this 0.18% bias would spoil."""
        ob, om = self.ombh2, self.ombh2 + self.omch2 + self.omnuh2
        g1 = 0.0783 * ob**-0.238 / (1 + 39.5 * ob**0.763)
        g2 = 0.560 / (1 + 21.1 * ob**1.81)
        return Z_STAR_CALIB * 1048 * (1 + 0.00124 * ob**-0.738) * (1 + g1 * om**g2)

    def theta_star(self):
        """Angular size of the sound horizon at recombination, r_s / D_M."""
        zs = self.z_star()
        return self.r_s(zs) / float(self.D_M(zs))

    def r_s(self, z):
        """Comoving sound horizon at redshift z (Mpc)."""
        Rb = 3 * self.ombh2 / (4 * OMEGA_GAMMA_H2)
        f = lambda a: 1.0 / (a**2 * self.E(1 / a - 1) * np.sqrt(3 * (1 + Rb * a)))
        return C_KMS / self.H0 * quad(f, 1e-8, 1 / (1 + z), limit=200, epsabs=0, epsrel=1e-8)[0]

    def r_d(self):
        """Sound horizon at the baryon drag epoch (Aubourg+ 2015, Eq. 16)."""
        wcb, wb, wnu = self.ombh2 + self.omch2, self.ombh2, self.omnuh2
        return 55.154 * np.exp(-72.3 * (wnu + 0.0006) ** 2) / (wcb**0.25351 * wb**0.12807)

    def cmb_distance_priors(self):
        zs = self.z_star()
        DM = float(self.D_M(zs)) if zs <= self._z[-1] else None
        R = np.sqrt(self.om) * self.H0 * DM / C_KMS
        lA = np.pi * DM / self.r_s(zs)
        return {"R": R, "l_A": lA, "omega_b": self.ombh2, "n_s": self.ns}

    # ------------------------------------------------------------ growth and power
    def growth(self, z):
        """Linear growth factor D(z), normalised to D = a deep in matter
        domination, from the growth ODE (valid for any w0, wa)."""
        if not hasattr(self, "_growth_a"):
            # matter (cdm + baryons + massive neutrinos) and dark energy only:
            # the 'D = a in matter domination' normalisation assumes no radiation
            om_tot = self.om
            de = lambda a: (1 - om_tot) * a ** (-3 * (1 + self.w0 + self.wa)) * np.exp(-3 * self.wa * (1 - a))
            Eg = lambda a: np.sqrt(om_tot * a**-3 + de(a))

            def rhs(lna, y):
                a = np.exp(lna)
                E = Eg(a)
                eps = 1e-5
                dlnE = (np.log(Eg(a * np.exp(eps))) - np.log(E)) / eps
                om_a = om_tot * a**-3 / E**2
                D, dD = y
                return [dD, -(2 + dlnE) * dD + 1.5 * om_a * D]
            a0 = 1e-3
            sol = solve_ivp(rhs, [np.log(a0), 0.0], [a0, a0], dense_output=True, rtol=1e-8, atol=1e-12)
            self._growth_a = sol
        lna = np.log(1.0 / (1.0 + np.asarray(z, dtype=float)))
        return self._growth_a.sol(lna)[0]

    def growth_rate(self, z):
        lna = np.log(1.0 / (1.0 + np.asarray(z, dtype=float)))
        self.growth(0.0)
        D, dD = self._growth_a.sol(lna)
        return dD / D

    def transfer_eh_nowiggle(self, k):
        """Eisenstein & Hu 1998 no-wiggle transfer function; k in 1/Mpc."""
        om0h2 = (self.ombh2 + self.omch2)
        fb = self.ombh2 / om0h2
        theta = T_CMB / 2.7
        s = 44.5 * np.log(9.83 / om0h2) / np.sqrt(1 + 10 * self.ombh2**0.75)
        alpha = 1 - 0.328 * np.log(431 * om0h2) * fb + 0.38 * np.log(22.3 * om0h2) * fb**2
        gamma_eff = om0h2 / self.h * (alpha + (1 - alpha) / (1 + (0.43 * k * s) ** 4))
        q = k / self.h * theta**2 / gamma_eff            # EH98 Eq. 28
        L0 = np.log(2 * np.e + 1.8 * q)
        C0 = 14.2 + 731 / (1 + 62.5 * q)
        return L0 / (L0 + C0 * q**2)

    def delta2_linear(self, k, z=0.0):
        """Dimensionless linear matter power Delta^2(k, z); k in 1/Mpc."""
        k = np.asarray(k, dtype=float)
        h0_inv_mpc = self.H0 / C_KMS
        prim = self.As * (k / 0.05) ** (self.ns - 1)
        # massive neutrinos free-stream and suppress small-scale power by
        # ~8 f_nu (Hu, Eisenstein & Tegmark 1998); smooth step above k_fs
        f_nu = self.omnuh2 / (self.ombh2 + self.omch2 + self.omnuh2)
        k_fs = 0.02 * self.h
        nu_supp = 1 - 8 * f_nu * (k / k_fs) ** 2 / (1 + (k / k_fs) ** 2)
        return (4 / 25) * prim * (k / h0_inv_mpc) ** 4 * (self.growth(z) / self.om) ** 2 \
            * self.transfer_eh_nowiggle(k) ** 2 * nu_supp

    def P_linear(self, k, z=0.0):
        k = np.asarray(k, dtype=float)
        return 2 * np.pi**2 * self.delta2_linear(k, z) / k**3

    def sigma_R(self, R_mpc, z=0.0):
        lnk = np.linspace(np.log(1e-5), np.log(1e3), 3000)
        k = np.exp(lnk)
        R = np.atleast_1d(R_mpc)[:, None]
        x = k[None, :] * R
        W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
        out = np.sqrt(np.trapezoid(self.delta2_linear(k, z)[None, :] * W**2, lnk, axis=1))
        return out if np.ndim(R_mpc) else float(out[0])

    def sigma8(self, z=0.0):
        return self.sigma_R(8.0 / self.h, z)

    def S8(self, alpha=0.5):
        return self.sigma8() * (self.om / 0.3) ** alpha

    # ------------------------------------------------------------ halos
    def halo_mass_function(self, M, z):
        """Sheth & Tormen (1999) dn/dlnM [1/Mpc^3] for masses M [Msun]."""
        M = np.atleast_1d(np.asarray(M, dtype=float))
        rho_m = RHO_CRIT_H2 * self.h**2 * self.ocb
        R = (3 * M / (4 * np.pi * rho_m)) ** (1 / 3)
        sig = self.sigma_R(R, z)
        dlns_dlnM = np.gradient(np.log(sig), np.log(M))
        nu = 1.686 / sig
        A, a, p = 0.3222, 0.707, 0.3
        f = A * np.sqrt(2 * a / np.pi) * (1 + (a * nu**2) ** -p) * nu * np.exp(-a * nu**2 / 2)
        return f * rho_m / M * np.abs(dlns_dlnM)

    # ------------------------------------------------------------ BBN
    def bbn(self):
        """Primordial abundances from the Steigman fitting formulas
        (Steigman 2012, arXiv:1208.0032): helium-4 mass fraction Yp,
        deuterium D/H and lithium-7 Li/H."""
        eta10 = 273.9 * self.ombh2
        S = np.sqrt(1 + 7 * (self.Neff - 3.046) / 43)
        Yp = 0.2485 + 0.0016 * ((eta10 - 6) + 100 * (S - 1))
        etaD = eta10 - 6 * (S - 1)
        DH = 2.55e-5 * (6 / etaD) ** 1.6
        etaLi = eta10 - 3 * (S - 1)
        LiH = 4.82e-10 * (etaLi / 6) ** 2
        return {"Yp": Yp, "D_H": DH, "Li_H": LiH}


PLANCK_2018 = {"z_star": 1089.92, "r_s_star": 144.43, "r_d": 147.09, "theta100": 1.04110,
               "R": 1.7502, "sigma8": 0.8111, "Omega_m": 0.3153}


def validate():
    """Engine vs Planck 2018 (TT,TE,EE+lowE+lensing, A&A 641, A6, Table 2)
    at the Planck best-fit parameters. R is from Chen, Huang & Wang 2019."""
    c = Cosmology()
    p = c.cmb_distance_priors()
    ours = {"z_star": c.z_star(), "r_s_star": c.r_s(c.z_star()), "r_d": c.r_d(),
            "theta100": 100 * c.theta_star(), "R": p["R"], "sigma8": c.sigma8(), "Omega_m": c.om}
    return {k: (round(float(ours[k]), 5), PLANCK_2018[k],
                f"{100 * (float(ours[k]) / PLANCK_2018[k] - 1):+.3f}%") for k in ours}


def lya_compressed(c, z=3.0, k_skm=0.009):
    """Linear power amplitude Delta^2_L and logarithmic slope n_eff at the
    Lyman-alpha pivot (k = 0.009 s/km, z = 3), as used by McDonald+ 2005 and
    the PRIYA analysis. Converts k from s/km to 1/Mpc with H(z)/(1+z)."""
    k = k_skm * float(c.H(z)) / (1 + z)
    d2 = float(c.delta2_linear(k, z))
    eps = 0.01
    lp = np.log(c.P_linear(np.array([k * (1 - eps), k * (1 + eps)]), z))
    n_eff = (lp[1] - lp[0]) / (np.log(1 + eps) - np.log(1 - eps))
    return d2, float(n_eff)
