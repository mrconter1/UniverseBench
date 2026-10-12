"""Baseline: the standard model, cold dark matter (+ dark energy later).

Galaxies: visible matter plus an NFW halo. Test 01 allows no per-galaxy
freedom: halo mass comes from the galaxy's stellar mass (abundance
matching, Moster+ 2013) and concentration from simulations (Dutton & Maccio
2014). Test 02 lets it adjust halo mass and concentration per galaxy.

Clusters: total mass = baryons / baryon fraction. Cosmology says ~0.16
overall; feedback pushes gas out of small groups and cluster cores, so the
fraction rises with mass and radius. Two global trend slopes model that
known gas physics (fitted, not hand-set)."""
import numpy as np

from ..physics import (G, dispersion_from_mass, halo_mass_from_stars, nfw_concentration,
                       nfw_mass_enclosed, nfw_mass_projected, v_nfw)
from .common import COSMO, COSMO_CDM, CosmologyMixin, ELL, GeneralRelativity, LOCAL_UPS, UPS, v_baryon_sq


class LCDM(GeneralRelativity, CosmologyMixin):
    ACTION = {"type": "GR+minimal"}   # GR + cold dark matter + Lambda: test 37
    NAME = "Cold dark matter (ΛCDM)"
    GLOBAL_PARAMS = {**UPS, **ELL, **COSMO, **COSMO_CDM,
                     "f_b_cluster": (0.13, 0.03, 0.3),     # at r500, 1e13 Msun of baryons
                     "fb_mass_slope": (0.1, -0.5, 0.5),
                     "fb_radius_slope": (0.1, -0.5, 0.5)}
    LOCAL_PARAMS = {**LOCAL_UPS,
                    "dlogM": (0.0, -3.0, 4.0),             # wide: a free halo in test 02
                    "dlogc": (0.0, -1.5, 1.0)}

    def halo(self, mstar, lp=None):
        lp = lp or {"dlogM": 0.0, "dlogc": 0.0}
        m200 = halo_mass_from_stars(mstar) * 10 ** lp["dlogM"]
        return m200, nfw_concentration(m200) * 10 ** lp["dlogc"]

    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp["ups_scale"]
        m200, c = self.halo(g.stellar_mass(u), lp)
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        return np.sqrt(np.clip(vb2 + v_nfw(R, m200, c) ** 2, 0, None))

    def predict_cluster_mass(self, c, gp):
        return c.M_bar / self.cluster_baryon_fraction(c.M_bar, gp, c.delta)

    def predict_dwarf_dispersion(self, d, gp):
        """Stars plus an NFW halo from abundance matching (extrapolated to
        these tiny stellar masses), no tidal stripping: the naive standard
        expectation."""
        m200, c = self.halo(d.M_bar)
        return dispersion_from_mass(d.M_bar / 2 + nfw_mass_enclosed(d.r_half, m200, c), d.r_half)

    def predict_wide_binary_boost(self, gN_si, g_ext_si, gp):
        return 1.0          # dark matter density is negligible inside a binary

    def predict_lensing_vcirc(self, M_bar, R_kpc, gp):
        """Baryons plus the galaxy's own NFW halo (abundance matching). The
        data give only baryonic mass, used here as the stellar mass, which
        overestimates halo masses for gas-rich (or hot-gas-rich) galaxies.
        Naive: no contribution from neighbouring halos (two-halo term),
        which matters beyond ~300 kpc."""
        R = np.asarray(R_kpc)
        m200, c = self.halo(M_bar)
        return np.sqrt(G * (M_bar + nfw_mass_enclosed(R, m200, c)) / R)

    def predict_lens_mass(self, L, gp):
        """Stars plus the projected mass of an NFW halo from abundance
        matching (Chabrier stellar masses, as in Moster+ 2013). No
        adiabatic contraction."""
        m200, c = self.halo(L.Mstar_chab)
        return gp["imf_ell"] * L.Mstar_proj_chab + nfw_mass_projected(L.r_Ein, m200, c)

    def cluster_baryon_fraction(self, M_bar, gp, delta=500):
        return (gp["f_b_cluster"] * (M_bar / 1e13) ** gp["fb_mass_slope"]
                * (delta / 500) ** (-gp["fb_radius_slope"]))

    def predict_bullet_offset(self, s, gp):
        """Collisionless dark matter travels with the galaxies. Its amount
        follows from the model's own cluster baryon fraction."""
        M_bar = s.M_gas + s.M_star
        dm = M_bar * (1 / self.cluster_baryon_fraction(M_bar, gp) - 1)
        return s.M_gas / (s.M_gas + s.M_star + dm)

    def dm_particle(self, gp):
        """Generic cold dark matter: effectively infinite thermal-relic mass,
        no specified couplings (so no predicted direct-detection signal)."""
        return {"type": "cold", "wdm_mass_keV": float("inf"), "sigma_SI_cm2": 0.0,
                "sigma_v_cm3s": 0.0, "mass_GeV": None}
