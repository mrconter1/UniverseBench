"""Baseline: standard gravity (Newton / general relativity), visible matter
only. No dark matter, no modified gravity. The original 'missing mass'
problem: expected to fail wherever dark matter is invoked."""
import numpy as np

from ..physics import G, dispersion_from_mass
from .common import LOCAL_UPS, UPS, v_baryon_sq


class VisibleOnly:
    NAME = "Standard gravity, visible matter only"
    GLOBAL_PARAMS = dict(UPS)
    LOCAL_PARAMS = dict(LOCAL_UPS)

    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp["ups_scale"]
        return np.sqrt(np.clip(v_baryon_sq(g, u, 1.4 * u), 0, None))

    def predict_cluster_mass(self, c, gp):
        return c.M_bar

    def predict_dwarf_dispersion(self, d, gp):
        return dispersion_from_mass(d.M_bar / 2, d.r_half)

    def predict_wide_binary_boost(self, gN_si, g_ext_si, gp):
        return 1.0

    def predict_lensing_vcirc(self, M_bar, R_kpc, gp):
        return np.sqrt(G * M_bar / np.asarray(R_kpc))    # all baryons well inside 50 kpc
