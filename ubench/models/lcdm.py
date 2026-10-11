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

from ..physics import (dispersion_from_mass, halo_mass_from_stars, nfw_concentration,
                       nfw_mass_enclosed, v_nfw)
from .common import LOCAL_UPS, UPS, v_baryon_sq


class LCDM:
    NAME = "Cold dark matter (ΛCDM)"
    GLOBAL_PARAMS = {**UPS,
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
        fb = (gp["f_b_cluster"]
              * (c.M_bar / 1e13) ** gp["fb_mass_slope"]
              * (c.delta / 500) ** (-gp["fb_radius_slope"]))
        return c.M_bar / fb

    def predict_dwarf_dispersion(self, d, gp):
        """Stars plus an NFW halo from abundance matching (extrapolated to
        these tiny stellar masses), no tidal stripping: the naive standard
        expectation."""
        m200, c = self.halo(d.M_bar)
        return dispersion_from_mass(d.M_bar / 2 + nfw_mass_enclosed(d.r_half, m200, c), d.r_half)
