"""Candidate C, dipolar dark matter / gravitational polarization (Blanchet &
Le Tiec 2008-09; Blanchet & Heisenberg): a medium that behaves like cold dark
matter on cosmological scales but reproduces MOND inside galaxies and
clusters. Benchmark proxy: the published-a0 MOND baseline everywhere below
cosmological scales, plus a cold-dark-matter cosmology (the medium's
linear-order behaviour). The medium carries no particle: no direct or
indirect detection signal."""
from ubench.models.common import COSMO, COSMO_CDM, CosmologyMixin
from ubench.models.mond import MONDFixed


class DipolarDM(MONDFixed, CosmologyMixin):
    NAME = "Dipolar dark matter (MOND galaxies + CDM cosmology)"
    GLOBAL_PARAMS = {**MONDFixed.GLOBAL_PARAMS, **COSMO, **COSMO_CDM}

    def dm_particle(self, gp):
        return {"type": "polarizable medium", "wdm_mass_keV": float("inf"),
                "sigma_SI_cm2": 0.0, "sigma_v_cm3s": 0.0, "mass_GeV": None}


MODEL = DipolarDM()
