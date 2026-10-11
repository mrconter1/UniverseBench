"""
Score models on the HIDDEN galaxies (run this yourself, never the agent).

    BENCH_SALT='my-private-secret' python score.py              # baselines
    BENCH_SALT='my-private-secret' python score.py my_model.py  # + a candidate

A candidate file must define MODEL = <instance with the models.py interface>.
"""
import importlib.util
import json
import sys

import benchmark as B
import models as M

candidates = list(M.BASELINES)
for path in sys.argv[1:]:
    spec = importlib.util.spec_from_file_location("candidate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    candidates.append(mod.MODEL)

gals = B.load_all()
pub, hid = B.split(gals)
print(f"{len(gals)} galaxies: {len(pub)} public, {len(hid)} hidden\n")
results = [B.run(m, gals) for m in candidates]
results.sort(key=lambda r: -r["TOTAL"])

print("\nLEADERBOARD (hidden data, % of points within 2 sigma; ~95 = perfect)")
for i, r in enumerate(results, 1):
    gap = r["A_predict_public"]["within_2sigma_pct"] - r["A_predict_hidden"]["within_2sigma_pct"]
    print(f"{i}. {r['model']:30s} TOTAL {r['TOTAL']:5.1f}   "
          f"predict {r['A_predict_hidden']['within_2sigma_pct']:5.1f}   "
          f"extrapolate {r['B_extrapolate_hidden']['within_2sigma_pct']:5.1f}   "
          f"clusters {r['C_clusters_hidden']['within_2sigma_pct']:5.1f}   "
          f"overfit gap {gap:+.1f}   knobs: {r['n_global_params']} global, "
          f"{r['n_local_params']} per galaxy")

with open("results.json", "w") as f:
    json.dump(results, f, indent=1)
