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
    GUARD_RAIL = False   # True: a gate every viable theory must pass; counts
                         # toward a model's lowest score but not its median
    SIMPLIFIED = ""      # non-empty: a known simplification, flagged in the scorecard
    N_BOUND = 0          # points that are one-sided limits (a perfect theory
                         # always passes them); used for the noise ceiling

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

    def n_bound(self):
        return self.N_BOUND

    def scored(self, z):
        """Which residuals count toward the score (default: all). Tests may
        keep some points for fitting only."""
        return z

    def score(self, model, gp, salt):
        part = "hidden" if self.HIDDEN_SPLIT else "all"
        z = self.residuals(model, gp, part, salt)
        if z is None:
            return dict(NO_PREDICTION)
        z = np.asarray(self.scored(np.asarray(z, dtype=float)), dtype=float)
        z[~np.isfinite(z)] = 1e6
        out = summarize(z)
        out["in_sample"] = (not self.HIDDEN_SPLIT) and self.USED_IN_FIT
        return out

    def noise_ceiling(self, n_points, sims=20000, seed=0):
        """What a perfect theory would score on this test: Gaussian points
        land within 2 sigma 95.4% of the time, one-sided limits always pass.
        Returns the expected score and the 5th-95th percentile range, which
        is wide for tests with only a few points."""
        nb = min(self.n_bound(), n_points)
        ng = n_points - nb
        if n_points == 0:
            return {"expected": 0.0, "p5": 0.0, "p95": 0.0}
        rng = np.random.default_rng(seed)
        inside = (np.abs(rng.standard_normal((sims, ng))) < 2).sum(axis=1) + nb
        pct = 100.0 * inside / n_points
        return {"expected": round(float(pct.mean()), 1),
                "p5": round(float(np.percentile(pct, 5)), 1),
                "p95": round(float(np.percentile(pct, 95)), 1)}


def _vec_to_dict(x, spec):
    return {k: float(np.clip(v, lo, hi)) for v, (k, (_, lo, hi)) in zip(x, spec.items())}


def fit_global(model, tests, salt):
    """Fit GLOBAL_PARAMS jointly on the public data of all tests. Each test
    contributes its mean chi^2 per point, so tests weigh equally.
    Independent groups of constants (see _param_groups) are fitted
    separately, which is exact and much faster."""
    spec = getattr(model, "GLOBAL_PARAMS", {})
    if not spec:
        return {}
    tests = [t for t in tests if t.USED_IN_FIT]
    x = np.array([v[0] for v in spec.values()], dtype=float)
    bounds_all = [(lo, hi) for (_, lo, hi) in spec.values()]

    for idx, group_tests in _param_groups(model, tests, salt, spec, x):
        if not group_tests:
            continue

        def f(sub, idx=idx, group_tests=group_tests):
            xx = x.copy()
            xx[idx] = sub
            gp = _vec_to_dict(xx, spec)
            total = 0.0
            for t in group_tests:
                term = t.fit_term(model, gp, salt)
                if term is not None:
                    total += term
            return total

        bounds = [bounds_all[i] for i in idx]
        best = minimize(f, x[idx], method="Powell", bounds=bounds,
                        options={"maxfev": 6000, "xtol": 1e-4, "ftol": 1e-7})
        # polish with Nelder-Mead from Powell's answer
        nm = minimize(f, best.x, method="Nelder-Mead",
                      options={"maxfev": 3000, "xatol": 1e-5, "fatol": 1e-8})
        x[idx] = nm.x if nm.fun < best.fun else best.x
    return _vec_to_dict(x, spec)


def _param_groups(model, tests, salt, spec, x0):
    """Split the global constants into independent groups: two constants
    belong together if some test responds to both. Detected by nudging each
    constant and recording which tests' fit terms change."""
    base_gp = _vec_to_dict(x0, spec)
    base = {t.ID: t.fit_term(model, base_gp, salt) for t in tests}
    touches = []
    for i, (_, lo, hi) in enumerate(spec.values()):
        step = 0.013 * (hi - lo)
        x = np.array(x0, dtype=float)
        x[i] = x0[i] + step if x0[i] + step <= hi else x0[i] - step
        gp = _vec_to_dict(x, spec)
        hit = set()
        for t in tests:
            b, v = base[t.ID], t.fit_term(model, gp, salt)
            if (b is None) != (v is None) or (b is not None and abs(v - b) > 1e-12 * max(1.0, abs(b))):
                hit.add(t.ID)
        touches.append(hit)
    parent = list(range(len(touches)))

    def find(a):
        while parent[a] != a:
            a = parent[a]
        return a

    for i in range(len(touches)):
        for j in range(i + 1, len(touches)):
            if touches[i] & touches[j]:
                parent[find(i)] = find(j)
    groups = {}
    for i in range(len(touches)):
        groups.setdefault(find(i), []).append(i)
    return [(idx, [t for t in tests if t.ID in set().union(*(touches[i] for i in idx))])
            for idx in groups.values()]


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
