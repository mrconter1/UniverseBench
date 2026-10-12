"""
Test 15: the cosmic microwave background, full power spectra.

Data: Planck 2018 plik-lite binned TT (ℓ 30–2508), TE and EE (ℓ 30–1996)
spectra with their full 613×613 covariance, two low-ℓ TT bins, and the
low-ℓ EE optical depth τ = 0.0544 ± 0.0073 (data/cmb_planck_lite/README.md).

Prediction: the model's cosmology is passed to the CLASS Boltzmann code
(model.class_params(gp) -> CLASS input dict; CosmologyMixin provides a
standard one), which computes the lensed spectra. Unlike the compressed
summary it replaces, this does not presuppose LCDM perturbations: any physics
CLASS can express (massive or sterile neutrinos, evolving dark energy,
curvature, dark radiation, ...) is judged on the spectra themselves.

Score: % of the 616 decorrelated (whitened) points within 2 sigma. Planck's
best-fit LCDM gives 96.4%. In-sample (‡).

Fitting: one CLASS run takes ~4 s, so this test is EXPENSIVE: the framework
fits constants on the fast tests first, then polishes the constants this
test depends on (model.CLASS_DEPENDS) with it included (core.fit_global).
"""
import numpy as np

from ..core import Test

_LIK = None


class CMBSpectra(Test):
    ID = "15_cmb"
    LABEL = "Microwave background (Planck spectra)"
    STAGE = 3
    HIDDEN_SPLIT = False
    DATA_DIRS = ("cmb_planck_lite",)
    EXPENSIVE = True

    def load(self):
        global _LIK
        from ..cmb_lite import PlanckLite
        if _LIK is None:
            _LIK = PlanckLite()
        self.lik = _LIK

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "class_params"):
            return None
        from ..class_engine import run_class
        params = model.class_params(gp)
        if params is None:
            return None
        out = run_class(params)
        if out is None:
            return np.full(616, 1e6)
        return self.lik.residuals(out["Dl_tt"], out["Dl_te"], out["Dl_ee"], params.get("tau_reio", 0.0544))
