"""
UniverseBench core: hidden splits, the Test base class, joint fitting of
global constants, and running a model through every test.

Every test is a table of measured numbers with error bars. A model must
supply a prediction for each; how it computes it is its own business.

Rules enforced here:
  * Hidden data: tests with many data points are split by a secret salt
    (BENCH_SALT) into public (for fitting / development) and hidden (for
    scoring). Tests with only a handful of numbers cannot be hidden; they
    are fitted and scored on everything and marked "in-sample".
  * One set of constants: a model's GLOBAL_PARAMS are fitted once, jointly
    on the public data of every test, then used everywhere.
  * No prediction = 0: if a model has no method for a test, it scores 0
    there (shown as "—").
"""
import hashlib
import os

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SPLIT_SALT = os.environ.get("BENCH_SALT", "demo-salt")


def data_path(*parts):
    return os.path.join(DATA, *parts)


def is_hidden(key, salt):
    return int(hashlib.sha256((salt + str(key)).encode()).hexdigest(), 16) % 2 == 1


def summarize(z):
    """Summary of standardized residuals z = (prediction - data) / sigma."""
    z = np.asarray(np.concatenate(z) if isinstance(z, list) else z, dtype=float)
    if len(z) == 0:
        return {"within_2sigma_pct": 0.0, "n_points": 0}
    return {
        "within_2sigma_pct": round(100 * float(np.mean(np.abs(z) < 2)), 1),
        "reduced_chi2_median": round(float(np.median(z**2)), 3),
        "n_points": int(len(z)),
    }


NO_PREDICTION = {"within_2sigma_pct": 0.0, "no_prediction": True, "n_points": 0}


def default_locals(model):
    return {k: v[0] for k, v in getattr(model, "LOCAL_PARAMS", {}).items()}


class Test:
    """Base class. Subclasses set the class attributes and implement
    residuals(model, gp, part, salt) returning standardized residuals for
    part in {"public", "hidden", "all"}, or None if the model makes no
    prediction."""
    ID = ""            # e.g. "05_dwarfs"
    LABEL = ""         # short column title
    STAGE = 1
    STATUS = "settled"  # or "contested"
    HIDDEN_SPLIT = True  # False: too few points to hide; scored in-sample
    USED_IN_FIT = True   # False: pure out-of-sample test
    DATA_DIRS = ()       # folders under data/ this test reads

    def __init__(self):
        self._loaded = False

    def ensure_loaded(self):
        if not self._loaded:
            self.load()
            self._loaded = True

    def load(self):
        pass

    def residuals(self, model, gp, part, salt):
        raise NotImplementedError

    def export_public(self, dest_data, salt):
        """Copy the data an agent may see. Default: everything (tests whose
        data cannot be hidden). Tests with hidden splits override this."""
        import shutil
        for d in self.DATA_DIRS:
            shutil.copytree(data_path(d), os.path.join(dest_data, d), dirs_exist_ok=True)

    # ---- used by the framework
    def fit_term(self, model, gp, salt):
        if not self.USED_IN_FIT:
            return None
        part = "public" if self.HIDDEN_SPLIT else "all"
        z = self.residuals(model, gp, part, salt)
        if z is None or len(z) == 0:
            return None
        z = np.asarray(z, dtype=float)
        return float(np.mean(z**2)) if np.all(np.isfinite(z)) else 1e12

    def score(self, model, gp, salt):
        part = "hidden" if self.HIDDEN_SPLIT else "all"
        z = self.residuals(model, gp, part, salt)
        if z is None:
            return dict(NO_PREDICTION)
        z = np.asarray(z, dtype=float)
        z[~np.isfinite(z)] = 1e6
        out = summarize(z)
        out["in_sample"] = (not self.HIDDEN_SPLIT) and self.USED_IN_FIT
        return out


def _vec_to_dict(x, spec):
    return {k: float(np.clip(v, lo, hi)) for v, (k, (_, lo, hi)) in zip(x, spec.items())}


def fit_global(model, tests, salt):
    """Fit GLOBAL_PARAMS jointly on the public data of all tests. Each test
    contributes its mean chi^2 per point, so tests weigh equally."""
    spec = getattr(model, "GLOBAL_PARAMS", {})
    if not spec:
        return {}

    def f(x):
        gp = _vec_to_dict(x, spec)
        total = 0.0
        for t in tests:
            term = t.fit_term(model, gp, salt)
            if term is not None:
                total += term
        return total

    x0 = np.array([v[0] for v in spec.values()], dtype=float)
    bounds = [(lo, hi) for (_, lo, hi) in spec.values()]
    best = minimize(f, x0, method="Powell", bounds=bounds,
                    options={"maxiter": 20000, "xtol": 1e-5, "ftol": 1e-7})
    # polish with Nelder-Mead from Powell's answer
    nm = minimize(f, best.x, method="Nelder-Mead",
                  options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-8})
    x = nm.x if nm.fun < best.fun else best.x
    return _vec_to_dict(x, spec)


def run(model, tests, salt=SPLIT_SALT, verbose=False):
    for t in tests:
        t.ensure_loaded()
    gp = fit_global(model, tests, salt)
    results = {t.ID: t.score(model, gp, salt) for t in tests}
    vals = [results[t.ID]["within_2sigma_pct"] for t in tests]
    out = {
        "model": model.NAME,
        "global_params": {k: round(v, 5) for k, v in gp.items()},
        "n_global_params": len(getattr(model, "GLOBAL_PARAMS", {})),
        "n_local_params": len(getattr(model, "LOCAL_PARAMS", {})),
        "tests": results,
        "LOWEST": min(vals), "MEDIAN": float(np.median(vals)), "HIGHEST": max(vals),
    }
    if verbose:
        s = "  ".join(f"{t.ID[:2]}={results[t.ID]['within_2sigma_pct']:5.1f}" for t in tests)
        print(f"{model.NAME:40s} {s}")
    return out
