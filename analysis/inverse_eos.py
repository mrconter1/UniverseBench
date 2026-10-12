"""Inverse design, step 3: if galaxies obey the radial acceleration relation
EXACTLY, what dark-matter density and pressure does each one need, and is
there one equation of state behind them?

For each SPARC quality-1 galaxy: g_tot = RAR(g_bar) (a0 = 1.2e-10 m/s^2,
McGaugh+ 2016 form) gives M_dm(<r) = (g_tot - g_bar) r^2 / G and the
required density rho(r) (spherical). Hydrostatic equilibrium dP/dr = -rho g_tot,
integrated inward from the last point (boundary: isothermal tail,
P = rho v^2 / 2), gives the required dark-matter 'temperature'
sigma^2(r) = P/rho. Then we ask which variable sigma^2 is a function of:
  * the local density rho        -> a dark fluid with a universal equation of state
  * the local baryonic field g_bar -> a local coupling to baryons
  * the galaxy's baryonic mass M_b -> a global thermostat
The smallest scatter (across 87 galaxies, all radii) points to the mechanism."""
import sys, json
sys.path[:0] = [".", "analysis"]
import numpy as np
from inverse_design_laws import GL, VB2, LMB, KPC
from ubench.physics import G, rar_nu

a0 = 1.2e-10 / KPC
rows = []
for g in GL:
    R = g.R
    gb = VB2[g.name] / R
    gt = gb * rar_nu(np.clip(gb, 1e-9, None) / a0)
    Mdm = np.maximum.accumulate(np.clip((gt - gb) * R ** 2 / G, 0, None))
    r = np.linspace(R[0], R[-1], 300)
    M = np.interp(r, R, Mdm)
    lnM = np.log(np.maximum(M, 1e-3))
    rho = np.gradient(M, r) / (4 * np.pi * r ** 2)
    gtr = np.interp(r, R, gt)
    gbr = np.interp(r, R, gb)
    good = (rho > 0) & (gbr > 0)
    if good.sum() < 50:
        continue
    r, rho, gtr, gbr = r[good], rho[good], gtr[good], gbr[good]
    P = np.empty_like(r)
    P[-1] = rho[-1] * gtr[-1] * r[-1] / 2
    for i in range(len(r) - 2, -1, -1):
        P[i] = P[i + 1] + 0.5 * (rho[i] * gtr[i] + rho[i + 1] * gtr[i + 1]) * (r[i + 1] - r[i])
    s2 = P / rho
    sel = slice(None, None, 10)
    for x in zip(np.log10(rho[sel]), np.log10(s2[sel]), np.log10(gbr[sel] * KPC), r[sel] / r[-1]):
        rows.append((g.name, LMB[g.name]) + x)

names = sorted(set(r[0] for r in rows))
gi = np.array([names.index(r[0]) for r in rows])
lMb, lrho, ls2, lgb, rr = (np.array([r[i] for r in rows]) for i in range(1, 6))
print(f"{len(names)} galaxies, {len(rows)} radial samples")


def scatter_vs(x_cols, label):
    X = np.column_stack([np.ones_like(ls2)] + [c for c in x_cols] + [c ** 2 for c in x_cols])
    b, *_ = np.linalg.lstsq(X, ls2, rcond=None)
    res = ls2 - X @ b
    per_gal = np.array([res[gi == j].mean() for j in range(len(names))])
    print(f"  sigma^2 as a function of {label:32s}: total rms {res.std():.3f} dex, between-galaxy {per_gal.std():.3f}")


print("Required dark-matter temperature sigma^2 = P/rho, scatter around the best quadratic law in:")
scatter_vs([lrho], "local density (equation of state)")
scatter_vs([lgb], "local baryonic field g_bar")
scatter_vs([lMb], "baryonic mass M_b (global)")
scatter_vs([lrho, lMb], "density and M_b")
scatter_vs([lgb, lMb], "g_bar and M_b")
# Isothermality: how much does sigma^2 vary within a galaxy?
within = np.array([ls2[gi == j].std() for j in range(len(names))])
print(f"  within-galaxy variation of sigma^2: median {np.median(within):.3f} dex (0 = exactly isothermal)")
b = np.polyfit(lMb, ls2, 1)
print(f"  global law: log sigma^2 = {b[1]:.2f} + {b[0]:.3f} log M_b  (MOND scaling predicts slope 0.5)")
