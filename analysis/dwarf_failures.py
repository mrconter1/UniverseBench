"""Which dwarf galaxies do the two dark-matter solver models miss (test 05,
all dwarfs, models' hidden-run constants)? Binned by stellar mass and by
host (satellite of the Milky Way / M31 vs isolated)."""
import sys, importlib.util, copy
sys.path.insert(0, ".")
import numpy as np
from ubench.tests import all_tests


def load(p):
    s = importlib.util.spec_from_file_location(p, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


gf, cs = load("submissions/solver_galform/my_model.py"), load("submissions/solver_claude/my_model.py")
GP = {"halo_A": 11.55989, "halo_B": 0.8545, "halo_lfloor": 9.8397, "dlogc": -0.11012}
t = [x for x in all_tests() if x.ID == "05_dwarfs"][0]; t.ensure_loaded()
rows = []
for d in t.dwarfs:
    if d.is_limit:
        continue
    out = []
    for m, gp in ((gf.MODEL, GP), (cs.MODEL, {})):
        p = m.predict_dwarf_dispersion(d, gp)
        du = copy.copy(d); du.M_bar = d.M_bar * 10 ** d.log_Mb_err
        pu = m.predict_dwarf_dispersion(du, gp)
        e = d.sigma_ep if p > d.sigma_obs else d.sigma_em
        out.append((p - d.sigma_obs) / np.hypot(e, pu - p))
    rows.append((np.log10(d.M_bar), d.host in ("MW", "M31"), out[0], out[1], d.name))
r = rows
print(f"{len(r)} dwarfs with measured dispersions. z = (pred - obs)/err; '>0' = model too high")
for lab, sel in [("log Mb < 4", lambda x: x[0] < 4), ("4-5", lambda x: 4 <= x[0] < 5), ("5-6", lambda x: 5 <= x[0] < 6),
                 ("6-7", lambda x: 6 <= x[0] < 7), (">= 7", lambda x: x[0] >= 7),
                 ("satellites", lambda x: x[1]), ("isolated/field", lambda x: not x[1])]:
    s = [x for x in r if sel(x)]
    if not s:
        continue
    zg, zc = np.array([x[2] for x in s]), np.array([x[3] for x in s])
    print(f"  {lab:15s} n={len(s):3d} | galform: median z {np.median(zg):+.2f}, within 2s {100*np.mean(abs(zg)<2):3.0f}% "
          f"| cSIDM: median z {np.median(zc):+.2f}, within 2s {100*np.mean(abs(zc)<2):3.0f}%")
worst = sorted(r, key=lambda x: -abs(x[2]))[:8]
print("largest galform misses:", ", ".join(f"{x[4]} ({x[2]:+.1f})" for x in worst))
