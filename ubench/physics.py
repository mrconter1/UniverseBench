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


MSTAR_SCATTER_DEX = 0.15           # scatter of log M* at fixed halo mass (Moster+ 2013)


def halo_mass_from_stars(mstar):
    """Typical (log-mean) halo mass hosting a central galaxy of stellar mass
    mstar. Naively inverting the mean M*(Mh) relation is biased high where
    that relation is flat (massive galaxies), because small halos vastly
    outnumber big ones. So we average over halo mass with weight
    halo abundance x P(M* | Mh), using dn/dlnMh ~ Mh^-0.9 exp(-Mh/2e14)
    and 0.15 dex lognormal scatter. Below ~1e8 Msun of stars Moster+2013
    is an extrapolation."""
    lm = np.log10(max(float(mstar), 1e2))
    return float(10 ** np.interp(lm, _LMS_TABLE, _LMH_TABLE))


def _mean_log_halo(lm):
    lmh = np.log10(_MH_GRID)
    prior = _MH_GRID ** -0.9 * np.exp(-_MH_GRID / 2e14)
    like = np.exp(-0.5 * ((np.log10(_MS_GRID) - lm) / MSTAR_SCATTER_DEX) ** 2)
    w = prior * like
    if w.sum() <= 0:
        return float(lmh[np.argmin(np.abs(np.log10(_MS_GRID) - lm))])
    return float(np.sum(w * lmh) / np.sum(w))


# precomputed once: log M* -> mean log Mh
_LMS_TABLE = np.arange(2.0, 12.6, 0.005)
_LMH_TABLE = np.array([_mean_log_halo(x) for x in _LMS_TABLE])


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


def nfw_mass_projected(R, m200, c):
    """Mass inside a cylinder of projected radius R [kpc] for an NFW halo
    (Bartelmann 1996; Golse & Kneib 2002)."""
    rs = nfw_r200(m200) / c
    x = np.atleast_1d(np.asarray(R, dtype=float)) / rs
    g = np.empty_like(x)
    lo, hi, eq = x < 1 - 1e-8, x > 1 + 1e-8, np.abs(x - 1) <= 1e-8
    g[lo] = np.log(x[lo] / 2) + 2 / np.sqrt(1 - x[lo]**2) * np.arctanh(np.sqrt((1 - x[lo]) / (1 + x[lo])))
    g[hi] = np.log(x[hi] / 2) + 2 / np.sqrt(x[hi]**2 - 1) * np.arctan(np.sqrt((x[hi] - 1) / (1 + x[hi])))
    g[eq] = 1 + np.log(0.5)
    mu_c = np.log(1 + c) - c / (1 + c)
    out = m200 / mu_c * g
    return out if np.ndim(R) else float(out[0])
