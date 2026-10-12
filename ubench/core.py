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
    def fit_chi2(self, model, gp, salt):
        """(sum of chi^2, number of points) on the fitting part, or None."""
        if not self.USED_IN_FIT:
            return None
        part = "public" if self.HIDDEN_SPLIT else "all"
        z = self.residuals(model, gp, part, salt)
        if z is None or len(z) == 0:
            return None
        z = np.asarray(z, dtype=float)
        return (float(np.sum(z**2)) if np.all(np.isfinite(z)) else 1e12 * len(z), len(z))

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


# Cosmology data sets are one joint likelihood: inside this block every data
# point weighs the same (summed chi^2), so Planck's 616 spectral points are
# not outvoted by 13 BAO numbers when constants are inferred. The block as a
# whole keeps the weight it would have under equal per-test weighting.
# (Lesson from idea 8, analysis/RESULTS.md: equal per-test weighting made a
# worse summed fit look preferred.)
JOINT_IDS = {"15_cmb", "16_bbn", "17_bao", "18_supernovae", "21_cluster_counts", "28_inflation",
             "29_neutrinos", "30_curvature", "31_growth_rate", "34_standard_sirens"}


def _objective(terms):
    """terms: list of (test ID, (chi2_sum, n) or None)."""
    total, js, jn, jk = 0.0, 0.0, 0, 0
    for tid, v in terms:
        if v is None:
            continue
        if tid in JOINT_IDS:
            js, jn, jk = js + v[0], jn + v[1], jk + 1
        else:
            total += v[0] / v[1]
    if jn:
        total += jk * js / jn
    return total


def fit_global(model, tests, salt):
    """Fit GLOBAL_PARAMS jointly on the public data of all tests. Each test
    contributes its mean chi^2 per point, so tests weigh equally, except the
    cosmology data sets (JOINT_IDS), which form one summed likelihood.

    Independent groups of constants (see _param_groups) are fitted
    separately, which is exact and much faster. Within a group, each test's
    term is cached by the values of the constants it responds to, so moving
    one constant only re-runs the tests that depend on it.

    The search is robust to step-like, multi-modal objectives (e.g. a
    threshold radius): every constant is scanned across its range, Powell
    runs from the model's start and from the scan's best points, and the
    answer is the best point evaluated anywhere. (scipy's bounded Powell
    line search never re-evaluates its current point, so on its own it can
    end WORSE than it started.)"""
    spec = getattr(model, "GLOBAL_PARAMS", {})
    if not spec:
        return {}
    tests = [t for t in tests if t.USED_IN_FIT]
    expensive = [t for t in tests if getattr(t, "EXPENSIVE", False)]
    tests = [t for t in tests if not getattr(t, "EXPENSIVE", False)]
    x = np.array([v[0] for v in spec.values()], dtype=float)
    bounds_all = [(lo, hi) for (_, lo, hi) in spec.values()]
    groups, touches = _param_groups(model, tests, salt, spec, x, return_touches=True)

    for idx, group_tests in groups:
        if not group_tests:
            continue
        bounds = [bounds_all[i] for i in idx]
        lo_b = np.array([b[0] for b in bounds])
        hi_b = np.array([b[1] for b in bounds])
        # positions (within idx) of the constants each test responds to
        deps = {t.ID: [j for j, i in enumerate(idx) if t.ID in touches[i]] for t in group_tests}
        cache = {t.ID: {} for t in group_tests}
        best = {"f": np.inf, "x": None}

        def f(sub, idx=idx, group_tests=group_tests, deps=deps, cache=cache, best=best,
              lo_b=lo_b, hi_b=hi_b):
            sub = np.clip(np.asarray(sub, dtype=float), lo_b, hi_b)
            xx = x.copy()
            xx[idx] = sub
            gp = None
            terms = []
            for t in group_tests:
                key = tuple(np.round(sub[deps[t.ID]], 12))
                c = cache[t.ID]
                if key not in c:
                    if gp is None:
                        gp = _vec_to_dict(xx, spec)
                    if len(c) > 20000:
                        c.clear()
                    c[key] = t.fit_chi2(model, gp, salt)
                terms.append((t.ID, c[key]))
            total = _objective(terms)
            if np.isfinite(total) and total < best["f"]:
                best["f"], best["x"] = total, sub.copy()
            return total

        x0 = np.clip(x[idx].copy(), lo_b, hi_b)
        f0 = f(x0)
        # scan each constant across its range (others at their starts)
        scan_best, combined = [], x0.copy()
        for j, (lo, hi) in enumerate(bounds):
            grid = lo + (hi - lo) * np.linspace(0.05, 0.95, 7)
            vals = []
            for g in grid:
                trial = x0.copy()
                trial[j] = g
                vals.append(f(trial))
            vals = np.array(vals, dtype=float)
            if np.isfinite(vals).any() and np.nanmin(vals) < f0 - 1e-6:
                k = int(np.nanargmin(vals))
                combined[j] = grid[k]
                scan_best.append((vals[k], j, grid[k]))
        starts = [x0]
        if scan_best:
            v, j, g = min(scan_best)
            single = x0.copy()
            single[j] = g
            starts.append(single)
            if len(scan_best) > 1:
                starts.append(combined)
        for x_start in starts:
            minimize(f, x_start, method="Powell", bounds=bounds,
                     options={"maxfev": 4000, "xtol": 1e-4, "ftol": 1e-7})
        # coordinate sweeps: scan each constant finely across its whole range
        # with the others at the best point so far, then re-run Powell. This
        # finds narrow basins of step-like constants (a threshold radius
        # crossing individual objects) that line searches step over. Cheap,
        # because only the tests depending on that constant are re-run.
        for _ in range(2):
            before = best["f"]
            for j, (lo, hi) in enumerate(bounds):
                for g in np.linspace(lo, hi, 41):
                    trial = best["x"].copy()
                    trial[j] = g
                    f(trial)
            if best["f"] < before - 1e-6:
                minimize(f, best["x"], method="Powell", bounds=bounds,
                         options={"maxfev": 4000, "xtol": 1e-4, "ftol": 1e-7})
            else:
                break
        # polish with Nelder-Mead from the best point found
        minimize(f, best["x"], method="Nelder-Mead",
                 options={"maxfev": 3000, "xatol": 1e-5, "fatol": 1e-8})
        x[idx] = best["x"]
    if expensive:
        x = _polish_expensive(model, tests, expensive, salt, spec, x, touches)
    return _vec_to_dict(x, spec)


POLISH_MAXFEV = int(os.environ.get("UBENCH_POLISH_MAXFEV", 120))
_POLISH_DIR = os.path.join(ROOT, ".polish_cache")


def _polish_expensive(model, tests, expensive, salt, spec, x, touches):
    """Second fitting stage for EXPENSIVE tests (e.g. the CLASS-based CMB
    spectra, ~4 s per evaluation). Only the constants such a test depends on
    (model.CLASS_DEPENDS) are moved, by a short Nelder-Mead starting from the
    stage-one fit, with the expensive test added to every fast test those
    constants affect. Accepted only if it improves the joint objective.
    Results are cached on disk: the cosmology fit is the same for every
    split unless a hidden-split test shares its constants."""
    import hashlib, inspect, json, pickle
    keys = list(spec)
    names = [k for k in keys if k in getattr(model, "CLASS_DEPENDS", ())]
    if not names:
        return x
    idx = [keys.index(k) for k in names]
    touched = set().union(*(touches[i] for i in idx))
    obj_tests = [t for t in tests if t.ID in touched] + list(expensive)
    lo_b = np.array([spec[k][1] for k in names])
    hi_b = np.array([spec[k][2] for k in names])

    def f(sub):
        xx = x.copy()
        xx[idx] = np.clip(sub, lo_b, hi_b)
        gp = _vec_to_dict(xx, spec)
        return _objective([(t.ID, t.fit_chi2(model, gp, salt)) for t in obj_tests])

    try:
        src = open(getattr(model, "_source_file", None) or inspect.getfile(type(model))).read()
    except (TypeError, OSError):
        src = model.NAME
    salt_dep = salt if any(t.HIDDEN_SPLIT for t in obj_tests) else ""
    key = hashlib.sha256(json.dumps([model.NAME, src, {k: list(v) for k, v in spec.items()}, names,
                                     list(np.round(x, 8)), [t.ID for t in obj_tests], salt_dep,
                                     POLISH_MAXFEV, sorted(JOINT_IDS)]).encode()).hexdigest()[:24]
    path = os.path.join(_POLISH_DIR, key + ".pkl")
    if os.path.exists(path):
        return pickle.load(open(path, "rb"))
    x0 = x[idx].copy()
    f0 = f(x0)
    steps = 0.015 * (hi_b - lo_b)
    simplex = [x0] + [np.clip(x0 + np.eye(len(x0))[i] * steps[i], lo_b, hi_b) for i in range(len(x0))]
    res = minimize(f, x0, method="Nelder-Mead",
                   options={"maxfev": POLISH_MAXFEV, "xatol": 1e-5, "fatol": 1e-3,
                            "initial_simplex": np.array(simplex)})
    out = x.copy()
    if res.fun < f0:
        out[idx] = np.clip(res.x, lo_b, hi_b)
    os.makedirs(_POLISH_DIR, exist_ok=True)
    pickle.dump(out, open(path, "wb"))
    return out


def _param_groups(model, tests, salt, spec, x0, return_touches=False):
    """Split the global constants into independent groups: two constants
    belong together if some test responds to both. Detected by nudging each
    constant and recording which tests' fit terms change."""

    base_gp = _vec_to_dict(x0, spec)
    base = {t.ID: t.fit_term(model, base_gp, salt) for t in tests}
    touches = []
    for i, (_, lo, hi) in enumerate(spec.values()):
        # probe a small nudge AND points spread across the whole range, so
        # threshold-like constants (which only matter when they cross an
        # object) are still detected
        probes = [x0[i] + 0.013 * (hi - lo) if x0[i] + 0.013 * (hi - lo) <= hi else x0[i] - 0.013 * (hi - lo)]
        probes += list(lo + (hi - lo) * np.array([0.1, 0.35, 0.65, 0.9]))
        hit = set()
        for val in probes:
            x = np.array(x0, dtype=float)
            x[i] = val
            gp = _vec_to_dict(x, spec)
            for t in tests:
                if t.ID in hit:
                    continue
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
    out = [(idx, [t for t in tests if t.ID in set().union(*(touches[i] for i in idx))])
           for idx in groups.values()]
    return (out, touches) if return_touches else out


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
