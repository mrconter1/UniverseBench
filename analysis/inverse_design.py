"""Inverse design: what dark-matter equilibrium does each galaxy require?

For every SPARC quality-1 galaxy (inclination >= 30 deg), fit a self-gravitating
ISOTHERMAL dark-matter halo in equilibrium inside the galaxy's own baryonic
potential (the fully thermalized limit of any collisional dark matter:
self-interacting, baryon-heated, superfluid at finite temperature ...):
    sigma0^2 d ln rho/dr = -G (M_b(<r) + M_dm(<r)) / r^2,   rho(0) = rho0.
Two numbers per galaxy, (rho0, sigma0); stellar M/L fixed at 0.463 (disc),
1.4x for bulges, as in the benchmark. Then ask the data which law links
(rho0, sigma0) to the baryons, and how tight it is. The scaling argument:
MOND's v^4 = G M_b a0 needs sigma0^4 proportional to M_b, galaxy by galaxy.
"""
import sys, importlib.util, json
sys.path.insert(0, ".")
import numpy as np
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G

spec = importlib.util.spec_from_file_location("csidm", "submissions/solver_claude/my_model.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
U = cs.UPS_DISK
T1 = {}
for line in open("data/sparc_table1/SPARC_Lelli2016c.mrt"):
    f = line.split()
    if len(f) == 19 and f[1].isdigit():
        T1[f[0]] = dict(inc=float(f[5]), Q=int(f[17]), Reff=float(f[9]), Rdisk=float(f[11]),
                        MHI=float(f[13]), L36=float(f[7]))
GALS = [g for g in load_sparc() if g.name in T1 and T1[g.name]["Q"] == 1 and T1[g.name]["inc"] >= 30]


def iso_vdm2(g, lrho, lsig, rmax_factor=1.0):
    """Vectorized RK4 for arrays of (log10 rho0 [Msun/kpc^3], log10 sigma0 [km/s]).
    Returns V_dm^2 at the galaxy's radii, shape (n_params, n_R)."""
    Mb = cs._mb_profile(g, U)
    r = np.concatenate([np.geomspace(1e-3, g.R[0], 40, endpoint=False), np.linspace(g.R[0], g.R[-1], 360)])
    rho0, s2 = 10.0 ** lrho, (10.0 ** lsig) ** 2
    y = np.zeros_like(rho0)                      # ln(rho/rho0)
    m = 4 / 3 * np.pi * r[0] ** 3 * rho0
    out = np.empty((len(r), len(rho0)))
    out[0] = m
    Mbr = Mb(r)

    def f(ri, Mbi, y, m):
        return -G * (Mbi + m) / (s2 * ri ** 2), 4 * np.pi * ri ** 2 * rho0 * np.exp(y)
    for i in range(len(r) - 1):
        h, rm = r[i + 1] - r[i], 0.5 * (r[i] + r[i + 1])
        Mbm = 0.5 * (Mbr[i] + Mbr[i + 1])
        k1 = f(r[i], Mbr[i], y, m)
        k2 = f(rm, Mbm, y + h / 2 * k1[0], m + h / 2 * k1[1])
        k3 = f(rm, Mbm, y + h / 2 * k2[0], m + h / 2 * k2[1])
        k4 = f(r[i + 1], Mbr[i + 1], y + h * k3[0], m + h * k3[1])
        y = y + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        m = m + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        out[i + 1] = m
    Mdm = np.array([np.interp(g.R, r, out[:, j]) for j in range(len(rho0))])
    return G * Mdm / g.R


def fit(g):
    vb2 = np.clip(v_baryon_sq(g, U, 1.4 * U), 0, None)
    best = None
    lr_c, ls_c, wr, ws = 7.5, 1.6, 2.5, 0.8
    for it in range(4):
        LR, LS = np.meshgrid(np.linspace(lr_c - wr, lr_c + wr, 41), np.linspace(ls_c - ws, ls_c + ws, 41))
        LR, LS = LR.ravel(), LS.ravel()
        v = np.sqrt(vb2 + iso_vdm2(g, LR, LS))
        chi = (((v - g.Vobs) / g.sigma) ** 2).sum(1)
        j = int(np.nanargmin(chi))
        best = (LR[j], LS[j], chi[j])
        lr_c, ls_c, wr, ws = LR[j], LS[j], wr / 4, ws / 4
    return best


if __name__ == "__main__":
    rows = []
    for g in GALS:
        lr, ls, chi = fit(g)
        Ld, Lb = cs._lum_total(g)
        Mb = U * Ld + 1.4 * U * Lb + cs._gas_mass(g)
        t = T1[g.name]
        rows.append(dict(name=g.name, lrho0=float(lr), lsig0=float(ls), chi2=float(chi), n=len(g.R),
                         logMb=float(np.log10(Mb)), Rdisk=t["Rdisk"], Reff=t["Reff"],
                         fgas=float(cs._gas_mass(g) / Mb), vflat=float(np.median(g.Vobs[-3:]))))
        print(f"{g.name:12s} log rho0 {lr:6.2f} log sigma0 {ls:5.2f}  chi2/n {chi/len(g.R):5.2f}", flush=True)
    json.dump(rows, open("analysis/inverse_design_fits.json", "w"), indent=1)
