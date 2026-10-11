"""Baseline: pure Newtonian gravity, visible matter only, and NO relativity.

Identical to the visible-only model for galaxies and clusters, but in the
precision tests it shows why general relativity was needed:
  * light bending and radar echo delay are half of GR's (PPN gamma = 0),
  * no extra precession of Mercury's orbit (2 + 2 gamma - beta = 0, so
    beta = 2 and Nordtvedt eta = 4 beta - gamma - 3 = 5),
  * no gravitational waves: binary pulsar orbits do not decay, and the
    gravitational-wave tests get no prediction,
  * no relativistic cosmology (no prediction for Stage 3).
  * no black hole shadows, no relativistic precession of S2, and no
    gravitational time dilation: clocks tick at the same rate everywhere.
Laboratory tests (equivalence principle, inverse-square law) are passed.
"""
from .visible import VisibleOnly
from .common import ELL, LOCAL_UPS, UPS

NEWTON_VALUES = {
    "gamma_minus_1": -1.0, "beta_minus_1": 1.0, "nordtvedt_eta": 5.0, "gdot_over_g_per_yr": 0.0,
    "pbdot_ratio_B1913+16": 0.0, "quadrupole_ratio_double_pulsar": 0.0,
    "gw_speed_minus_c_over_c": None, "graviton_mass_ev": None,
    "eotvos_eta_Ti_Pt": 0.0, "yukawa_alpha_at_40um": 0.0,
    # no event horizons or shadows; no relativistic precession; clocks tick
    # at the same rate everywhere (redshift deviation alpha = -1)
    "sgrA_shadow_deviation": None, "m87_shadow_to_dynamics_mass_ratio": None,
    "s2_schwarzschild_precession_fSP": 0.0, "s2_gravitational_redshift_f": 0.0,
    "redshift_violation_alpha_galileo": -1.0, "redshift_violation_alpha_skytree": -1.0,
}


class Newtonian(VisibleOnly):
    NAME = "Pure Newtonian gravity, visible matter only"
    GLOBAL_PARAMS = {**UPS, **ELL}
    LOCAL_PARAMS = dict(LOCAL_UPS)

    def predict_gravity(self, quantity, gp):
        return NEWTON_VALUES.get(quantity)

    def cosmology(self, gp):
        return None
