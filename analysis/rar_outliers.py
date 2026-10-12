"""Idea 5: which quality-1 SPARC galaxies sit off the RAR, and why?
Reuses the per-galaxy offset and error budget of rar_age_v2.py, adds
distance method (1 = Hubble flow, least reliable) and inclination."""
import sys, runpy, io, contextlib
sys.path.insert(0, '.')
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path("analysis/rar_age_v2.py")
rows, tau_ml, T1 = ns["rows"], ns["tau_ml"], ns["T1"]
fD = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    if len(f) == 19 and f[1].isdigit():
        fD[f[0]] = int(f[4])
q1 = [r for r in rows if r["Q"] == 1]
for r in q1:
    r["z"] = r["off"] / r["err"]; r["fD"] = fD[r["name"]]; r["inc"] = T1[r["name"]]["inc"]
print(f"quality-1 galaxies: {len(q1)}, intrinsic tau {tau_ml(q1):.3f} dex\n")
print("largest outliers (offset in units of its own error):")
for r in sorted(q1, key=lambda r: -abs(r["z"]))[:10]:
    meth = {1: "Hubble flow", 2: "TRGB", 3: "Cepheid", 4: "UMa cluster", 5: "SN"}[r["fD"]]
    print(f"  {r['name']:12s} offset {r['off']:+.2f} dex ({r['z']:+.1f} sigma)  type T={r['T']:2d}  inc {r['inc']:.0f}  distance: {meth:11s}  gas frac {r['fgas']:.2f}  logMb {r['lm']:.1f}")
good = [r for r in q1 if r["fD"] != 1]
hub = [r for r in q1 if r["fD"] == 1]
print(f"\nby distance method: Hubble-flow distances tau {tau_ml(hub):.3f} (n={len(hub)}) | "
      f"direct distances (TRGB, Cepheid, cluster, SN) tau {tau_ml(good):.3f} (n={len(good)})")
z = np.array([r["z"] for r in q1])
print(f"fraction of quality-1 galaxies beyond 3 sigma: {np.mean(abs(z)>3):.2f} (Gaussian with no intrinsic scatter: 0.003)")
