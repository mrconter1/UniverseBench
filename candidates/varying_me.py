"""Candidate 8: a different electron mass at recombination (Hubble tension).

LCDM plus one constant, m_e(early)/m_e(today), switching to today's value
at z = 50 (instantaneous; CLASS 'varying_fundamental_constants').
A heavier early electron makes hydrogen recombine earlier (z* ∝ m_e),
shrinking the sound horizon, r_d ∝ m_e^-0.617 (calibrated with CLASS: the
exponent is -0.61 to -0.63 across cosmologies and for m_e ratios 0.97-1.06).
That lets BAO and the CMB accept a higher H0.
* Test 15 (CMB spectra): exact, through CLASS.
* Fast engine (BAO, compressed CMB fit aid): z* and r_d rescaled by the
  calibrated power laws.
Galaxies, clusters, particles: the LCDM baseline. Action: GR + minimally
coupled matter; the constant's variation is an assumption about the
early-universe value of a Standard Model parameter, with no new long-range
field at late times (an instantaneous switch is a phenomenological
stand-in for a field that settles before z = 50).
"""
from ubench.cosmology import Cosmology
from ubench.models.lcdm import LCDM

RD_EXP = -0.617


class VaryingMeCosmology(Cosmology):
    def __init__(self, me_ratio=1.0, **kw):
        self.me_ratio = me_ratio
        super().__init__(**kw)

    def z_star(self):
        return super().z_star() * self.me_ratio

    def r_d(self):
        return super().r_d() * self.me_ratio ** RD_EXP


class VaryingMe(LCDM):
    NAME = "Varying electron mass at recombination"
    GLOBAL_PARAMS = {**LCDM.GLOBAL_PARAMS, "me_ratio": (1.0, 0.90, 1.10)}
    CLASS_DEPENDS = LCDM.CLASS_DEPENDS + ("me_ratio",)
    ACTION = {"type": "GR+minimal"}

    def class_params(self, gp):
        p = super().class_params(gp)
        if p is not None:
            p.update({"varying_fundamental_constants": "instantaneous",
                      "varying_me": gp["me_ratio"], "varying_transition_redshift": 50})
        return p

    def cosmology(self, gp):
        keys = ("H0", "ombh2", "omch2", "ns", "ln10As", "mnu")
        args = {k: gp[k] for k in keys if k in gp}
        args["me_ratio"] = gp["me_ratio"]
        key = tuple(sorted(args.items()))
        cache = self.__dict__.setdefault("_cosmo_cache", {})
        if key not in cache:
            if len(cache) > 8:
                cache.clear()
            cache[key] = VaryingMeCosmology(**args)
        return cache[key]


MODEL = VaryingMe()
