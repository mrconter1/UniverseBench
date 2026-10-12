"""Idea 10 (user's): a slowly rolling light field drives both evolving dark
energy and the dark-matter self-interaction, so sigma/m drifts with time and
a0 should drift with it. How sensitive is the emergent a0 to sigma/m?
Frozen cSIDM, abundance-matched halo link (analysis/a0_derived.py)."""
import sys, io, contextlib, runpy
sys.path.insert(0, '.')
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    ns = runpy.run_path("analysis/a0_derived.py", run_name="not_main")
m, evaluate, mh_of_mb = ns["m"], ns["evaluate"], ns["mh_of_mb"]
base = m.SIGMA0
out = []
for f in (0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 10.0):
    m.SIGMA0 = base * f
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        evaluate(mh_of_mb, "x")
    a0 = float(buf.getvalue().split("emergent a0 = ")[1].split()[0])
    out.append((f, a0))
    print(f"sigma/m x {f:5.1f}: emergent a0 = {a0:.3e} m/s^2")
m.SIGMA0 = base
f, a = np.array(out).T
s = np.gradient(np.log(a), np.log(f))
print("local slope d ln a0 / d ln sigma:", ", ".join(f"{x:.2f}" for x in s))
