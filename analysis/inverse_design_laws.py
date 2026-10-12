"""Inverse design, step 2: impose candidate LAWS for the halo's thermal state
(no per-galaxy freedom) and ask (a) how well they fit SPARC and (b) how much
scatter they produce around the radial acceleration relation (RAR), versus
the observed intrinsic 0.058 dex. Halo: isothermal dark matter in equilibrium
in the galaxy's baryonic potential (analysis/inverse_design.py)."""
import sys, json
sys.path.insert(0, "."); sys.path.insert(0, "analysis")
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from inverse_design import GALS, U, cs, iso_vdm2, T1
from ubench.models.common import v_baryon_sq
from ubench.physics import G
from ubench.physics import rar_nu

fits = {r["name"]: r for r in json.load(open("analysis/inverse_design_fits.json"))}
GL = [g for g in GALS if fits[g.name]["chi2"] / fits[g.name]["n"] < 3]
VB2 = {g.name: np.clip(v_baryon_sq(g, U, 1.4 * U), 0, None) for g in GL}
LMB = {g.name: fits[g.name]["logMb"] for g in GL}
LRD = {g.name: np.log10(fits[g.name]["Rdisk"]) for g in GL}
KPC = 1e6 / 3.0857e19


def model_v(g, lrho, lsig):
    return np.sqrt(VB2[g.name] + iso_vdm2(g, np.array([lrho]), np.array([lsig]))[0])


def rar_scatter(vs):
    gb, gm, gi = [], [], []
    for j, g in enumerate(GL):
        ok = VB2[g.name] > 0
        gb += list(VB2[g.name][ok] / g.R[ok]); gm += list(vs[j][ok] ** 2 / g.R[ok]); gi += [j] * int(ok.sum())
    gb, gm, gi = map(np.array, (gb, gm, gi))
    rar = lambda a0: gb / (-np.expm1(-np.sqrt(gb / a0)))
    la = minimize_scalar(lambda la: np.mean((np.log10(gm / rar(10 ** la))) ** 2), bounds=(-1, 4), method="bounded").x
    res = np.log10(gm / rar(10 ** la))
    offs = np.array([np.median(res[gi == j]) for j in range(len(GL))])
    return float(offs.std()), 10 ** la * KPC


def chi2(vs):
    return float(sum((((v - g.Vobs) / g.sigma) ** 2).sum() for v, g in zip(vs, GL)))


def law_sigma(la0, g):
    # MOND law: v_flat^4 = G Mb a0 with v_flat = sqrt(2) sigma0  ->  sigma0^4 = G Mb a0 / 4
    a0 = 10 ** la0 / KPC
    return 0.25 * np.log10(G * 10 ** LMB[g.name] * a0 / 4)


LAWS = {
    # name: (initial params, function params,g -> (lrho, lsig))
    "A: sigma0 MOND law, rho0 = const": ([np.log10(1.2e-10), 7.8], lambda p, g: (p[1], law_sigma(p[0], g))),
    "B: sigma0 MOND law, rho0 ~ Mb^a Rd^b": ([np.log10(1.2e-10), 7.8, 0.3, -0.9],
        lambda p, g: (p[1] + p[2] * (LMB[g.name] - 10) + p[3] * LRD[g.name], law_sigma(p[0], g))),
    "C: free power laws sigma0 ~ Mb^s Rd^t, rho0 ~ Mb^a Rd^b": ([1.6, 0.3, 0.0, 7.8, 0.3, -0.9],
        lambda p, g: (p[3] + p[4] * (LMB[g.name] - 10) + p[5] * LRD[g.name],
                      p[0] + p[1] * (LMB[g.name] - 10) + p[2] * LRD[g.name])),
}

if __name__ == "__main__":
    vs_fit = [model_v(g, fits[g.name]["lrho0"], fits[g.name]["lsig0"]) for g in GL]
    print(f"{len(GL)} galaxies, {sum(len(g.R) for g in GL)} points")
    print(f"per-galaxy fits (2 free per galaxy): chi2 {chi2(vs_fit):.0f}, RAR between-galaxy scatter {rar_scatter(vs_fit)[0]:.3f}")
    vs_mond = []
    for g in GL:
        gb = VB2[g.name] / g.R
        vs_mond.append(np.sqrt(np.clip(gb * rar_nu(gb / (1.2e-10 / KPC)), 0, None) * g.R))
    print(f"MOND RAR (a0 = 1.2e-10, 0 free):     chi2 {chi2(vs_mond):.0f}")
    out = {}
    for name, (p0, fn) in LAWS.items():
        obj = lambda p: chi2([model_v(g, *fn(p, g)) for g in GL])
        res = minimize(obj, p0, method="Nelder-Mead", options={"maxfev": 600, "xatol": 1e-3, "fatol": 0.5})
        vs = [model_v(g, *fn(res.x, g)) for g in GL]
        sc, a0 = rar_scatter(vs)
        out[name] = dict(params=res.x.tolist(), chi2=res.fun, scatter=sc)
        print(f"{name}: chi2 {res.fun:.0f}, model RAR scatter {sc:.3f} dex (observed 0.058), params {np.round(res.x, 3)}"
              + (f", a0 = {10**res.x[0]:.2e}" if name[0] in "AB" else ""), flush=True)
    json.dump(out, open("analysis/inverse_design_laws.json", "w"), indent=1)
