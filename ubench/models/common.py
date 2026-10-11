"""Pieces shared by the baseline models."""
import numpy as np

UPS = {"ups_disk": (0.5, 0.1, 1.5)}           # global mass-to-light (3.6 um)
LOCAL_UPS = {"ups_scale": (1.0, 0.33, 3.0)}   # per-galaxy correction
# stellar-mass factor for giant ellipticals relative to a Chabrier IMF
# (1.0 = Chabrier, ~1.78 = Salpeter)
ELL = {"imf_ell": (1.3, 0.8, 2.0)}


def v_baryon_sq(g, ups_disk, ups_bul):
    """Newtonian speed^2 from visible matter. Signed squares: SPARC's gas
    term can be negative in central holes (net outward pull). Galaxies
    with ups_fixed (the Milky Way) already carry real stellar masses."""
    if getattr(g, "ups_fixed", False):
        ups_disk = ups_bul = 1.0
    s = lambda v: v * np.abs(v)
    return s(g.Vgas) + ups_disk * s(g.Vdisk) + ups_bul * s(g.Vbul)


class GeneralRelativity:
    """Predictions of general relativity for the precision gravity tests
    (11-14): every PPN deviation is 0 and every ratio to GR is 1."""
    GR_VALUES = {"pbdot_ratio_B1913+16": 1.0, "quadrupole_ratio_double_pulsar": 1.0,
                 "m87_shadow_to_dynamics_mass_ratio": 1.0, "s2_schwarzschild_precession_fSP": 1.0,
                 "s2_gravitational_redshift_f": 1.0, "spacetime_dimensions": 4.0}

    def predict_gravity(self, quantity, gp):
        return self.GR_VALUES.get(quantity, 0.0)


# Cosmological parameters (Stage 3). Starting points near Planck 2018.
COSMO = {"H0": (67.4, 50.0, 90.0),
         "ombh2": (0.0224, 0.010, 0.040),
         "ns": (0.965, 0.85, 1.10),
         "ln10As": (3.044, 2.0, 4.0)}
COSMO_CDM = {"omch2": (0.120, 0.01, 0.40)}


class CosmologyMixin:
    """Builds (and caches) a Cosmology from the global parameters. Subclasses
    set COSMO_FIXED for parameters they hold fixed (e.g. omch2 = 0)."""
    COSMO_FIXED = {}

    def cosmology(self, gp):
        from ..cosmology import Cosmology
        keys = ("H0", "ombh2", "omch2", "ns", "ln10As", "w0", "wa", "mnu")
        args = {k: gp[k] for k in keys if k in gp}
        args.update(self.COSMO_FIXED)
        key = tuple(sorted(args.items()))
        cache = self.__dict__.setdefault("_cosmo_cache", {})
        if key not in cache:
            if len(cache) > 8:
                cache.clear()
            cache[key] = Cosmology(**args)
        return cache[key]
