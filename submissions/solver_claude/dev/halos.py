"""Halo machinery for the solver: baryon-based halo assignment and the
Di Cintio et al. 2014 (DC14) feedback-dependent profile."""
import numpy as np
from ubench.physics import G, nfw_concentration, nfw_r200, H0_KPC

LOG10E = np.log10(np.e)


def dc14_shape(X):
    """(alpha, beta, gamma) of the DC14 profile for X = log10(M*/Mhalo).
    Calibrated for -4.1 < X < -1.3; outside it we clamp X (below -4.1 the
    profile tends to NFW: alpha=1, beta=3, gamma=1 approximately)."""
    X = float(np.clip(X, -4.1, -1.3))
    a = 2.94 - np.log10((10 ** (X + 2.33)) ** -1.08 + (10 ** (X + 2.33)) ** 2.29)
    b = 4.23 + 1.34 * X + 0.26 * X ** 2
    g = -0.06 + np.log10((10 ** (X + 2.56)) ** -0.68 + 10 ** (X + 2.56))
    return a, b, max(g, 0.0)


def dc14_conc_boost(X):
    """c_DC14 / c_NFW (Di Cintio+ 2014b): 1 + 3e-5 exp(3.4 (X + 4.5))."""
    X = float(np.clip(X, -4.1, -1.3))
    return 1.0 + 3e-5 * np.exp(3.4 * (X + 4.5))


class Profile:
    """Generalised (alpha, beta, gamma) halo normalised to mass m200 inside r200,
    with scale radius from r_-2 = r200 / c. Enclosed mass tabulated on a
    log grid for speed."""

    def __init__(self, m200, c, a=1.0, b=3.0, g=1.0, h=0.7):
        self.m200, self.r200 = m200, nfw_r200(m200)
        r_m2 = self.r200 / c
        self.rs = r_m2 / ((2 - g) / (b - 2)) ** (1 / a)
        x = np.logspace(-4, np.log10(30 * self.r200 / self.rs), 600)
        rho = x ** -g * (1 + x ** a) ** (-(b - g) / a)
        integrand = rho * x ** 3                       # dM ~ rho x^2 dx = rho x^3 dlnx
        cum = np.concatenate([[rho[0] * x[0] ** 3 / (3 - g)],
                              np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(np.log(x)))])
        cum[1:] += cum[0]
        self.lx, self.lcum = np.log(x), np.log(cum)
        norm = np.exp(np.interp(np.log(self.r200 / self.rs), self.lx, self.lcum))
        self.scale = m200 / norm

    def mass(self, r):
        lr = np.log(np.maximum(np.asarray(r, dtype=float), 1e-12) / self.rs)
        inner = lr < self.lx[0]
        out = self.scale * np.exp(np.interp(lr, self.lx, self.lcum))
        return out

    def vcirc(self, r):
        r = np.asarray(r, dtype=float)
        return np.sqrt(G * self.mass(r) / r)


def halo_from_baryons(mb, p):
    """Inverse of a double power law M_b(M_h) = 2 N M_h / ((M_h/M1)^-beta + (M_h/M1)^gam),
    with global constants p = (log_N, log_M1, beta, gam). Returns M_h."""
    lN, lM1, beta, gam = p
    lmh = np.linspace(6, 16, 2001)
    mh = 10 ** lmh
    mbs = mh * 2 * 10 ** lN / ((mh / 10 ** lM1) ** -beta + (mh / 10 ** lM1) ** gam)
    lmb = np.log10(mbs)
    # monotonic for beta>0, gam<1
    return 10 ** np.interp(np.log10(mb), lmb, lmh)


def half_light_radius(g):
    """Stellar (disc+bulge) half-light radius [kpc] from the SPARC surface
    brightness profile, by cumulative integration over the measured radii."""
    R = g.R
    sb = g.SBdisk + 1.4 * g.SBbul
    ring = 2 * np.pi * R * np.gradient(R) * sb
    cum = np.cumsum(ring)
    if cum[-1] <= 0:
        return float(R[len(R) // 2])
    return float(np.interp(0.5 * cum[-1], cum, R))


class CoreNFW:
    """Read, Agertz & Collins 2016 coreNFW: M(<r) = M_NFW(<r) * tanh(r/r_c)^n."""

    def __init__(self, m200, c, r_c, n=1.0):
        from ubench.physics import nfw_mass_enclosed
        self.m200, self.c, self.r_c, self.n = m200, c, r_c, n
        self._m = nfw_mass_enclosed

    def mass(self, r):
        r = np.asarray(r, dtype=float)
        return self._m(r, self.m200, self.c) * np.tanh(r / self.r_c) ** self.n

    def vcirc(self, r):
        r = np.asarray(r, dtype=float)
        return np.sqrt(G * self.mass(r) / r)
