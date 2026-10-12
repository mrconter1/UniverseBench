"""Idea 6: dark matter heated (or cooled) by baryons, Famaey, Khoury & Penco
2017 (arXiv:1712.01316) and Famaey, Khoury, Penco & Sharma 2019
(arXiv:1912.07626, "BIDM").

The papers' content, as far as the benchmark can use it:
* Where baryons and dark matter have reached thermal equilibrium (relaxation
  time ~ dynamical time; heating rate dE/dt / m = C a0 v rho_b/rho, their
  eq. 28), the dark matter settles so that the total field obeys the mass
  discrepancy-acceleration relation (their eqs. 36, 46-47, 59-61). a0 is an
  INPUT ("fixed empirically"), not derived; the scatter is not predicted.
* Systems that have not equilibrated keep ordinary collisionless-like halos:
  galaxy clusters explicitly (rho ~ 10 rho_b), and outer halos where the
  density is low.
* No particle mass or cross-section is specified; cosmology is CDM-like.

Benchmark proxy (a BEST CASE for the idea): MOND's relation in the
equilibrated, baryon-dominated regions (rotation curves 01/02/04 and the
inner lensing of ellipticals 08); cold dark matter everywhere else
(clusters, Bullet Cluster, outer-halo lensing 07, dwarf spheroidals,
cosmology). The dwarf choice is a judgement call: rho_b/rho is small in
dwarf spheroidals, so heating is weak; a variant with MOND dwarfs is
`BIDMMondDwarfs`.
"""
from ubench.models.lcdm import LCDM
from ubench.models.mond import MOND


class BIDM(LCDM):
    NAME = "Baryon-heated dark matter (FKP proxy)"
    ACTION = {"type": "other", "note": "no Lagrangian specified (FKP leave the microphysics open)"}
    GLOBAL_PARAMS = {**LCDM.GLOBAL_PARAMS, "a0_e10": MOND.GLOBAL_PARAMS["a0_e10"]}

    LOCAL_PARAMS = dict(MOND.LOCAL_PARAMS)   # no per-galaxy halo freedom needed
    a0 = MOND.a0
    predict = MOND.predict                      # equilibrated disks: the MDAR
    predict_lens_mass = MOND.predict_lens_mass  # inner ellipticals: baryon-dominated


class BIDMMondDwarfs(BIDM):
    NAME = "Baryon-heated dark matter (FKP proxy, MOND dwarfs)"
    predict_dwarf_dispersion = MOND.predict_dwarf_dispersion


MODEL = BIDMMondDwarfs()
