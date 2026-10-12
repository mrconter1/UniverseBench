"""Idea 9: cosmologically coupled black holes (k = 3) vs LCDM on the
cosmology tests, public split 0. Same constants in both models (the CCBH
dark-energy shape is fixed by the star-formation history), so the comparison
is direct. Reports scores and the summed chi^2 of CMB + BAO + SN."""
import sys, pickle, time
import numpy as np
sys.path[:0] = [".", "candidates"]
from ubench import core
from ubench.tests import all_tests
from ubench.models.lcdm import LCDM
from coupled_bh import MODEL as CCBH, W0, WA

keep = {"15_cmb", "15a_cmb_compressed", "16_bbn", "17_bao", "18_supernovae", "19_hubble_tension",
        "20_s8", "21_cluster_counts", "28_inflation", "29_neutrinos", "30_curvature",
        "31_growth_rate", "34_standard_sirens"}
tests = [t for t in all_tests() if t.ID in keep]
out = {}
for m in (LCDM(), CCBH):
    t0 = time.time()
    r = core.run(m, tests, salt=core.SPLIT_SALT)
    gp = {k: v for k, v in r["global_params"].items()}
    chi = {}
    for t in tests:
        if t.ID in ("15_cmb", "17_bao", "18_supernovae"):
            z = t.residuals(m, gp, "all", core.SPLIT_SALT)          # all points
            chi[t.ID] = float((np.asarray(z) ** 2).sum())
    out[m.NAME] = (r, chi)
    print(m.NAME, f"{time.time()-t0:.0f}s", flush=True)
    print("  params", {k: gp[k] for k in ("H0", "omch2", "ombh2") if k in gp})
    print("  scores", {t.ID: round(r["tests"][t.ID]["within_2sigma_pct"], 1) for t in tests})
    print("  chi2", {k: round(v, 1) for k, v in chi.items()}, "sum", round(sum(chi.values()), 1), flush=True)
print(f"CCBH CPL used for CLASS: w0 = {W0:.3f}, wa = {WA:.3f}")
pickle.dump(out, open("/tmp/ccbh_cosmo.pkl", "wb"))
