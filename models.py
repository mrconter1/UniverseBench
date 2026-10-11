"""
Baseline models. Every model implements the same interface:

  NAME            str
  GLOBAL_PARAMS   {name: (start, lo, hi)}  shared by ALL galaxies,
                  tuned on public galaxies only
  LOCAL_PARAMS    {name: (default, lo, hi)} per-galaxy knobs; the
                  default is used in test A, they are tuned on inner
                  points in test B
  predict(g, R, gp, lp) -> predicted rotation speed [km/s] at radii R [kpc]

Units: kpc, km/s, Msun. G in kpc (km/s)^2 / Msun.
"""
import numpy as np

G = 4.30091e-6
H0 = 0.07                          # 70 km/s/Mpc in km/s/kpc
A0_SI_TO_KPC = 1.0 / 3.0857e-14    # m/s^2 -> (km/s)^2/kpc


def v_baryon_sq(g, ups_disk, ups_bul):
    """Newtonian speed^2 from visible matter. Signed squares: SPARC's
    gas term can be negative in central holes (net outward pull)."""
    s = lambda v: v * np.abs(v)
    return s(g.Vgas) + ups_disk * s(g.Vdisk) + ups_bul * s(g.Vbul)


UPS = {"ups_disk": (0.5, 0.1, 1.5)}           # global mass-to-light
LOCAL_UPS = {"ups_scale": (1.0, 0.33, 3.0)}   # per-galaxy correction


class Newton:
    """Visible matter only, ordinary gravity. Expected to fail: this is
    the original 'missing mass' problem."""
    NAME = "Newton, visible matter only"
    GLOBAL_PARAMS = dict(UPS)
    LOCAL_PARAMS = dict(LOCAL_UPS)

    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp["ups_scale"]
        return np.sqrt(np.clip(v_baryon_sq(g, u, 1.4 * u), 0, None))

    def predict_cluster_mass(self, c, gp):
        return c.M_bar


class MOND:
    """Modified gravity (MOND), using the McGaugh 2016 radial
    acceleration relation as interpolating function. One universal
    acceleration scale a0 replaces dark matter entirely."""
    NAME = "MOND (one universal a0)"
    GLOBAL_PARAMS = {**UPS, "a0_e10": (1.2, 0.3, 4.0)}   # a0 in 1e-10 m/s^2
    LOCAL_PARAMS = dict(LOCAL_UPS)

    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp["ups_scale"]
        gbar = np.clip(v_baryon_sq(g, u, 1.4 * u), 1e-6, None) / R
        a0 = gp["a0_e10"] * 1e-10 * A0_SI_TO_KPC
        gobs = gbar / (1.0 - np.exp(-np.sqrt(gbar / a0)))
        return np.sqrt(gobs * R)

    def predict_cluster_mass(self, c, gp):
        # same law, same a0: spherical baryon mass -> Newtonian g -> MOND g
        gbar = G * c.M_bar / c.r_kpc**2
        a0 = gp["a0_e10"] * 1e-10 * A0_SI_TO_KPC
        gobs = gbar / (1.0 - np.exp(-np.sqrt(gbar / a0)))
        return gobs * c.r_kpc**2 / G


def _moster_mstar(mh):
    """Moster+2013 stellar-to-halo mass relation at z=0."""
    M1, N, b, gam = 10**11.59, 0.0351, 1.376, 0.608
    return mh * 2 * N / ((mh / M1) ** -b + (mh / M1) ** gam)


_MH_GRID = np.logspace(9, 15, 2000)
_MS_GRID = _moster_mstar(_MH_GRID)


def halo_mass_from_stars(mstar):
    """Invert the relation on its low-mass (monotonic) branch."""
    mstar = max(mstar, 1e5)
    i = np.searchsorted(_MS_GRID, mstar)
    i = min(max(i, 1), len(_MH_GRID) - 1)
    return float(_MH_GRID[i])


def concentration(m200):
    """Dutton & Maccio 2014 NFW concentration-mass relation, z=0."""
    h = 0.7
    return 10 ** (0.905 - 0.101 * np.log10(m200 * h / 1e12))


def v_nfw(R, m200, c):
    r200 = (G * m200 / (100 * H0**2)) ** (1 / 3)
    v200 = np.sqrt(G * m200 / r200)
    x = R / r200
    mu = lambda y: np.log(1 + y) - y / (1 + y)
    return v200 * np.sqrt(mu(c * x) / (x * mu(c)))


class LCDM_NFW:
    """Standard cold dark matter: visible matter plus an NFW halo.
    Test A uses no per-galaxy freedom: halo mass comes from the galaxy's
    stellar mass (abundance matching) and concentration from simulations.
    Test B lets it adjust halo mass and concentration per galaxy."""
    NAME = "Cold dark matter (NFW halo)"
    # Clusters: total mass = baryons / baryon fraction. Cosmology says ~0.16
    # overall; feedback pushes gas out of small groups and cluster cores, so
    # the fraction rises with mass and radius. Two global trend slopes model
    # that known gas physics (fitted, not hand-set).
    GLOBAL_PARAMS = {**UPS,
                     "f_b_cluster": (0.13, 0.03, 0.3),     # at r500, 1e13 Msun of baryons
                     "fb_mass_slope": (0.1, -0.5, 0.5),
                     "fb_radius_slope": (0.1, -0.5, 0.5)}

    def predict_cluster_mass(self, c, gp):
        fb = (gp["f_b_cluster"]
              * (c.M_bar / 1e13) ** gp["fb_mass_slope"]
              * (c.delta / 500) ** (-gp["fb_radius_slope"]))
        return c.M_bar / fb

    LOCAL_PARAMS = {**LOCAL_UPS,
                    # wide bounds: effectively a free halo in test B
                    "dlogM": (0.0, -3.0, 4.0),
                    "dlogc": (0.0, -1.5, 1.0)}

    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp["ups_scale"]
        Ld, Lb = g.luminosity()
        mstar = u * Ld + 1.4 * u * Lb
        m200 = halo_mass_from_stars(mstar) * 10 ** lp["dlogM"]
        c = concentration(m200) * 10 ** lp["dlogc"]
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        return np.sqrt(np.clip(vb2 + v_nfw(R, m200, c) ** 2, 0, None))


BASELINES = [Newton(), MOND(), LCDM_NFW()]
