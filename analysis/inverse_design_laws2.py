"""Inverse design, step 4: the required halo is nearly isothermal with
sigma^4 ~ M_b (analysis/inverse_eos.py). The remaining question is what
sets its central density. Physically motivated closures, one or two global
constants each, no per-galaxy freedom:
  E: core set by the MOND radius:  rho0 = k sigma0^2 / (G r_M^2), r_M = sqrt(G M_b / a0)
  F: core set by the baryons' size: r0 = lam * R_disk  (r0^2 = 9 sigma0^2 / (4 pi G rho0))
  S: singular isothermal sphere, V_dm = sqrt(2) sigma0 at all radii (no core)
sigma0^4 = G M_b a0 / 4 throughout (the deep-MOND limit is then automatic)."""
import sys, json
sys.path[:0] = [".", "analysis"]
import numpy as np
from scipy.optimize import minimize
from inverse_design_laws import GL, VB2, LMB, LRD, KPC, model_v, rar_scatter, chi2, law_sigma
from ubench.physics import G


def lrho_E(p, g):
    a0 = 10 ** p[0] / KPC
    s2 = 10 ** (2 * law_sigma(p[0], g))
    rM2 = G * 10 ** LMB[g.name] / a0
    return np.log10(10 ** p[1] * s2 / (G * rM2))


def lrho_F(p, g):
    s2 = 10 ** (2 * law_sigma(p[0], g))
    r0 = 10 ** p[1] * 10 ** LRD[g.name]
    return np.log10(9 * s2 / (4 * np.pi * G * r0 ** 2))


def v_S(p, g):
    s0 = 10 ** law_sigma(p[0], g)
    return np.sqrt(VB2[g.name] + 2 * s0 ** 2)


LAWS = {
    "E: core at the MOND radius": ([np.log10(1.2e-10), 0.0], lambda p, g: model_v(g, lrho_E(p, g), law_sigma(p[0], g))),
    "F: core = lam * R_disk": ([np.log10(1.2e-10), 0.0], lambda p, g: model_v(g, lrho_F(p, g), law_sigma(p[0], g))),
    "S: singular isothermal (no core)": ([np.log10(1.2e-10)], v_S),
}
out = {}
for name, (p0, fn) in LAWS.items():
    obj = lambda p: chi2([fn(p, g) for g in GL])
    res = minimize(obj, p0, method="Nelder-Mead", options={"maxfev": 150, "xatol": 1e-3, "fatol": 1})
    vs = [fn(res.x, g) for g in GL]
    sc, _ = rar_scatter(vs)
    out[name] = dict(params=res.x.tolist(), chi2=res.fun, scatter=sc)
    print(f"{name}: chi2 {res.fun:.0f}, model RAR scatter {sc:.3f} dex, a0 = {10**res.x[0]:.2e}, params {np.round(res.x, 3)}",
          flush=True)
json.dump(out, open("analysis/inverse_design_laws2.json", "w"), indent=1)
