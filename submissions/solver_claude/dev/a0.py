"""Which acceleration scale does the SIDM model produce? Fit the RAR form
g = g_b / (1 - exp(-sqrt(g_b/a0))) to the MODEL's own predicted accelerations
(not to the data), over all public galaxies."""
import sys, json
sys.path[:0] = ['.', 'dev']
import numpy as np
from scipy.optimize import minimize_scalar
from fitpar import MODELS, expand, GALS
import rc
from ubench.models.common import v_baryon_sq
KPC = 3.0857e19
x = json.load(open("dev/ck_sidm_all.json"))["x"]
f, keys, _, _ = MODELS["sidm"]
p = expand(dict(zip(keys, x)))
gb, gm, gd = [], [], []
for g in GALS:
    v = f(g, p)
    vb2 = np.clip(v_baryon_sq(g, p["ups"], 1.4 * p["ups"]), 1e-6, None)
    gb += list(vb2 / g.R); gm += list(v ** 2 / g.R); gd += list(g.Vobs ** 2 / g.R)
to_si = 1e6 / KPC
gb, gm, gd = np.array(gb) * to_si, np.array(gm) * to_si, np.array(gd) * to_si
ok = gb > 1e-13
def rar(gb, a0): return gb / (-np.expm1(-np.sqrt(gb / a0)))
for lab, y in [("model", gm), ("data", gd)]:
    r = minimize_scalar(lambda la: np.mean((np.log10(y[ok]) - np.log10(rar(gb[ok], 10 ** la))) ** 2), bounds=(-11, -9), method="bounded")
    res = np.log10(y[ok]) - np.log10(rar(gb[ok], 10 ** r.x))
    print(f"{lab:5s}: effective a0 = {10**r.x:.3e} m/s^2, scatter about RAR {np.std(res):.3f} dex")
H0 = 67.4e3 / 3.0857e22
print(f"cH0/2pi = {2.998e8*H0/2/np.pi:.3e} m/s^2;  MOND a0 = 1.2e-10")
print("\nlog g_bar bin | median log(g_model/g_obs) | n")
lb = np.log10(gb)
for lo in np.arange(-12.5, -8.9, 0.5):
    k = (lb >= lo) & (lb < lo + 0.5)
    if k.sum() > 5:
        print(f"  {lo:+.1f}..{lo+0.5:+.1f}   {np.median(np.log10(gm[k]/gd[k])):+.3f}   {k.sum()}")
