"""Candidate O: MOND plus hot dark matter (~11 eV sterile neutrinos; Angus
2009, Angus & Diaferio 2011).

* Galaxies, dwarfs, galaxy lensing, ellipticals: MOND with the published a0.
  Fermions cannot pack into galaxies (Tremaine-Gunn: the maximum
  phase-space density caps rho_nu <= g (2 pi)^(3/2) m^4 sigma^3 / h^3, tiny
  at galaxy velocity dispersions).
* Clusters (sigma ~ 1000 km/s): the cap is ~1e5 times higher, so neutrinos
  add real mass. Their profile needs a hydrostatic model, so (exactly like
  the LCDM baseline's baryon fraction) neutrino mass inside r is taken as
  eta x M_baryon x (M_baryon/1e13)^s_m x (delta/500)^-s_r with three fitted
  constants; MOND then acts on baryons + neutrinos. (A first version filled
  the volume at the Tremaine-Gunn cap with a self-consistent sigma; it runs
  away and was dropped.)
* Bullet Cluster: the neutrinos are collisionless, so they travel with the
  galaxies.
* Cosmology: a thermal sterile neutrino of mass m has Omega_nu h^2 = m / 94 eV,
  ~0.12 for 11 eV, which mimics cold dark matter in the expansion history
  (non-relativistic by z ~ 2e4). The background tests (CMB distances, BAO,
  supernovae, BBN) use a LCDM engine with omch2 = m/94 eV. Structure growth
  (sigma8, growth rate, cluster counts) is NOT computable here: hot dark
  matter erases small-scale power while MOND speeds up growth, and resolving
  the competition needs MOND N-body simulations. Those tests are reported
  as not computed.
"""
import numpy as np

from ubench.models.common import COSMO, CosmologyMixin
from ubench.models.mond import MONDFixed
from ubench.physics import G, rar_nu

H_PLANCK = 6.62607e-34
EV_KG = 1.78266e-36
MSUN, KPC = 1.98847e30, 3.0857e19


def rho_max_tg(m_ev, sigma_kms, g=2):
    """Tremaine-Gunn maximum mass density [Msun/kpc^3] for fermions of mass m."""
    m = m_ev * EV_KG
    s = sigma_kms * 1e3
    rho = g * (2 * np.pi) ** 1.5 * m ** 4 * s ** 3 / H_PLANCK ** 3        # kg/m^3
    return rho * KPC ** 3 / MSUN


class MONDHotDM(MONDFixed, CosmologyMixin):
    NAME = "MOND + 11 eV sterile neutrinos (hot DM)"
    GLOBAL_PARAMS = {**MONDFixed.GLOBAL_PARAMS, **COSMO,
                     "m_nu_ev": (11.0, 3.0, 40.0),
                     "eta_nu": (3.0, 0.0, 20.0),
                     "snu_mass": (0.0, -0.5, 0.5),
                     "snu_radius": (0.0, -0.5, 0.5)}

    CLASS_DEPENDS = ("H0", "ombh2", "ns", "ln10As", "tau_reio", "m_nu_ev")

    def class_params(self, gp):
        """Full CMB physics in CLASS: no cold dark matter; two massive
        neutrino species, the active one (0.06 eV) and a fully thermalised
        sterile one of mass m_nu_ev (same temperature, so it also adds
        ~1 to N_eff while relativistic). Linear theory with standard gravity:
        valid at recombination, where accelerations are far above a0."""
        return {"h": gp["H0"] / 100.0, "omega_b": gp["ombh2"], "omega_cdm": 1e-6,
                "n_s": gp["ns"], "ln_A_s_1e10": gp["ln10As"], "tau_reio": gp["tau_reio"],
                "N_ur": 2.0328, "N_ncdm": 2, "m_ncdm": f"0.06,{gp['m_nu_ev']:.6f}",
                "T_ncdm": "0.71611,0.71611"}

    def cosmology(self, gp):
        gp = dict(gp)
        gp["omch2"] = gp["m_nu_ev"] / 94.1
        return CosmologyMixin.cosmology(self, gp)

    def _nu_mass(self, M_b, gp, delta=500.0):
        return gp["eta_nu"] * M_b * (M_b / 1e13) ** gp["snu_mass"] * (delta / 500.0) ** (-gp["snu_radius"])

    def predict_cluster_mass(self, c, gp):
        M = c.M_bar + self._nu_mass(c.M_bar, gp, c.delta)
        return M * float(rar_nu(G * M / c.r_kpc ** 2 / self.a0(gp)))

    def predict_bullet_offset(self, s, gp):
        """Neutrinos travel with the galaxies (collisionless)."""
        m_nu = self._nu_mass(s.M_gas + s.M_star, gp)
        return s.M_gas / (s.M_gas + s.M_star + m_nu)

    def dm_particle(self, gp):
        return {"type": "hot (sterile neutrino)", "mass_GeV": gp["m_nu_ev"] * 1e-9,
                "wdm_mass_keV": gp["m_nu_ev"] * 1e-3, "sigma_SI_cm2": 0.0, "sigma_v_cm3s": 0.0}


MODEL = MONDHotDM()
