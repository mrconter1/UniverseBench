"""
Test 07: weak gravitational lensing around isolated galaxies.

Data: Mistele, McGaugh, Lelli, Schombert & Li 2024 (ApJL 969, L3), Table 2:
the average circular velocity of isolated KiDS-1000 lens galaxies in four
baryonic-mass bins, inferred from weak lensing over 50-300 kpc and 50-1000
kpc, about ten times farther out than rotation curves reach. The curves
stay roughly flat to ~1 Mpc. A published comment (arXiv:2407.18263) argues
that neighbouring halos (the "two-halo term") contaminate the outer parts.

Uncertainty: statistical error plus the paper's 0.1 dex stellar
mass-to-light systematic, propagated through each model's own prediction.
Too few numbers to hide, so fitted and scored in-sample.

Model method: predict_lensing_vcirc(M_bar, R_kpc, gp) -> circular velocity
[km/s] at radii R_kpc (array) around an isolated galaxy of baryonic mass
M_bar [Msun]. The test averages it over the paper's radial range.
"""
import csv

import numpy as np

from ..core import Test, data_path

ML_SYS_DEX = 0.1


class GalaxyLensing(Test):
    ID = "07_galaxy_lensing"
    LABEL = "Lensing around galaxies"
    HIDDEN_SPLIT = False
    DATA_DIRS = ("galaxy_lensing",)

    def load(self):
        path = data_path("galaxy_lensing", "mistele2024_btfr.csv")
        self.rows = [{k: float(v) for k, v in r.items()}
                     for r in csv.DictReader(l for l in open(path) if not l.startswith("#"))]

    def _vmean(self, model, M, r, gp):
        R = np.logspace(np.log10(r["R_min_kpc"]), np.log10(r["R_max_kpc"]), 25)
        return float(np.mean(model.predict_lensing_vcirc(M, R, gp)))

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_lensing_vcirc"):
            return None
        z = []
        for r in self.rows:
            M = 10 ** r["log10_Mb"]
            v = self._vmean(model, M, r, gp)
            v_up = self._vmean(model, M * 10**ML_SYS_DEX, r, gp)
            z.append((v - r["Vflat"]) / np.hypot(r["Vflat_err"], v_up - v))
        return np.array(z)
