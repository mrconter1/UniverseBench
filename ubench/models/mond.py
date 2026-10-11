"""Baseline: MOND (Milgrom 1983), using the McGaugh, Lelli & Schombert 2016
radial acceleration relation as interpolating function. One universal
acceleration scale a0 replaces dark matter entirely."""
import numpy as np

from ..physics import A0_SI_TO_KPC, G, dispersion_from_mass, rar_nu
from .common import LOCAL_UPS, UPS, v_baryon_sq


class MOND:
    NAME = "MOND (a0 fitted to all tests)"
    GLOBAL_PARAMS = {**UPS, "a0_e10": (1.2, 0.3, 4.0)}   # a0 in 1e-10 m/s^2
    LOCAL_PARAMS = dict(LOCAL_UPS)

    def a0(self, gp):
        """a0 in (km/s)^2/kpc."""
        return gp["a0_e10"] * 1e-10 * A0_SI_TO_KPC

    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp["ups_scale"]
        gbar = np.clip(v_baryon_sq(g, u, 1.4 * u), 1e-6, None) / R
        return np.sqrt(gbar * rar_nu(gbar / self.a0(gp)) * R)

    def predict_cluster_mass(self, c, gp):
        # same law, same a0: spherical baryon mass -> Newtonian g -> MOND g
        gbar = G * c.M_bar / c.r_kpc**2
        return c.M_bar * float(rar_nu(gbar / self.a0(gp)))

    def predict_dwarf_dispersion(self, d, gp):
        """Isolated deep-MOND: sigma^4 = (4/81) G M a0 (Milgrom 1994). If the
        host's field dominates, the dwarf is quasi-Newtonian with an
        effective G boosted by nu(g_ext/a0) (external field effect). Regime
        test as in McGaugh & Milgrom 2013: the external field dominates when
        it exceeds the dwarf's internal acceleration g_in = 3 sigma_iso^2 / r."""
        a0 = self.a0(gp)
        s_iso = (4.0 / 81.0 * G * d.M_bar * a0) ** 0.25
        g_in = 3.0 * s_iso**2 / d.r_half
        if d.g_ext <= g_in:
            return s_iso
        return dispersion_from_mass(float(rar_nu(d.g_ext / a0)) * d.M_bar / 2, d.r_half)


class MONDFixed(MOND):
    """MOND as its proponents use it: a0 is a measured constant of nature
    (1.2e-10 m/s^2, McGaugh+ 2016), not tuned to this benchmark."""
    NAME = "MOND (a0 fixed at 1.2, as published)"
    GLOBAL_PARAMS = dict(UPS)

    def a0(self, gp):
        return 1.2e-10 * A0_SI_TO_KPC
