import sys
sys.path[:0] = ['.', 'dev']
import numpy as np, hashlib
from rc import *
from fit_rc import fit, KEYS, BOUNDS
from scipy.optimize import minimize

def folds(k=5):
    f = [int(hashlib.sha256(("cv" + g.name).encode()).hexdigest(), 16) % k for g in GALS]
    return [([g for g, ff in zip(GALS, f) if ff != i], [g for g, ff in zip(GALS, f) if ff == i]) for i in range(k)]

def cv(vfun, x0, bounds, keys, label):
    zt = []
    for tr, te in folds():
        p, _ = fit(vfun, tr, x0, bounds, keys)
        zt += zs(lambda g: vfun(g, p), te)
    return report(zt, label + " (5-fold CV)")

def mond_fun(g, p):
    return v_mond(g, p["ups"], p["a0"] * 1e-10 * 3.0857e13)

if __name__ == "__main__":
    which = sys.argv[1:]
    if "mond" in which:
        cv(mond_fun, [0.5, 1.2], [(0.2, 1.0), (0.6, 2.0)], ["ups", "a0"], "MOND, ups and a0 fitted")
    if "nfw" in which:
        cv(lambda g, p: v_halo_model(g, p, "nfw"), [0.5, -1.0, 11.6, 1.4, 0.6, 0.0], BOUNDS, KEYS, "NFW + baryonic halo connection")
    if "dc14" in which:
        cv(lambda g, p: v_halo_model(g, p, "dc14"), [0.5, -1.0, 11.6, 1.4, 0.6, 0.0], BOUNDS, KEYS, "DC14 feedback profile")
    if "core" in which:
        cv(v_corenfw, [0.5, -1.0, 11.6, 1.4, 0.6, -0.3, 0.2], BOUNDS + [(-1.0, 1.0)], KEYS + ["leta"], "coreNFW, r_c = eta R_half")
