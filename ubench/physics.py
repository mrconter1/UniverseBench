"""
Shared physics helpers and constants used by tests and baseline models.

Units, unless noted: kpc, km/s, Msun. G in kpc (km/s)^2 / Msun.
"""
import numpy as np
from scipy.special import ellipk, i0, i1, k0, k1

G = 4.30091e-6                      # kpc (km/s)^2 / Msun
G_SI = 6.674e-11                    # m^3 kg^-1 s^-2
C_KMS = 299792.458                  # speed of light, km/s
KPC_M = 3.0857e19                   # metres per kpc
A0_SI_TO_KPC = 1.0 / 3.0857e-14     # (m/s^2) -> (km/s)^2/kpc


def v2_exponential_disc(R, S0, Rd):
    """Exact speed^2 of a razor-thin exponential disc (Freeman 1970).
    S0 in Msun/kpc^2, Rd in kpc."""
    y = R / (2 * Rd)
    return 4 * np.pi * G * S0 * Rd * y**2 * (i0(y) * k0(y) - i1(y) * k1(y))


def v2_thin_disc(R, sigma, a_max=80.0, n=40000):
    """Speed^2 of any razor-thin axisymmetric disc with surface density
    sigma(a) [Msun/kpc^2], by summing thin-ring potentials and
    differentiating. Validated against v2_exponential_disc to 0.02 km/s."""
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


def rar_nu(y):
    """MOND / radial-acceleration-relation boost g_obs/g_N for y = g_N/a0
    (McGaugh, Lelli & Schombert 2016)."""
    y = np.maximum(np.asarray(y, dtype=float), 1e-30)
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))


# ---------------------------------------------------------------- halos
H0_KPC = 0.07                       # 70 km/s/Mpc in km/s/kpc (galaxy-scale work)


def moster_mstar(mh):
    """Moster+2013 stellar-to-halo mass relation at z=0 (Msun)."""
    M1, N, b, gam = 10**11.59, 0.0351, 1.376, 0.608
    return mh * 2 * N / ((mh / M1) ** -b + (mh / M1) ** gam)


_MH_GRID = np.logspace(7, 15, 4000)
_MS_GRID = moster_mstar(_MH_GRID)


def halo_mass_from_stars(mstar):
    """Invert Moster+2013 on its low-mass (monotonic) branch. Below the
    calibrated range (~1e8 Msun of stars) this is an extrapolation."""
    mstar = max(float(mstar), 1e2)
    i = np.searchsorted(_MS_GRID, mstar)
    i = min(max(i, 1), len(_MH_GRID) - 1)
    return float(_MH_GRID[i])


def nfw_concentration(m200, h=0.7):
    """Dutton & Maccio 2014 NFW concentration-mass relation, z=0."""
    return 10 ** (0.905 - 0.101 * np.log10(m200 * h / 1e12))


def nfw_r200(m200):
    return (G * m200 / (100 * H0_KPC**2)) ** (1 / 3)


def nfw_mass_enclosed(r, m200, c):
    """3D mass inside radius r [kpc] for an NFW halo."""
    mu = lambda y: np.log(1 + y) - y / (1 + y)
    x = np.asarray(r) / nfw_r200(m200)
    return m200 * mu(c * x) / mu(c)


def v_nfw(R, m200, c):
    return np.sqrt(G * nfw_mass_enclosed(R, m200, c) / np.asarray(R))


def dispersion_from_mass(M_half, r_half):
    """Line-of-sight velocity dispersion of a dispersion-supported system
    whose mass inside the 3D half-light radius is M_half (Wolf+ 2010:
    M_half = 3 sigma^2 r_half / G). Msun, kpc -> km/s."""
    return np.sqrt(G * M_half / (3.0 * r_half))
