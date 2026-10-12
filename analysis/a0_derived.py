"""Idea 6: derive the galaxy-halo connection, and hence the emergent a0,
from cosmology instead of fitting it.
Abundance matching (centrals, no scatter): n_gal(>M_star) = n_halo(>M_h), with
  * stellar mass function: GAMA, Baldry et al. 2012 (MNRAS 421, 621), double
    Schechter, M* = 10^10.66, phi1 = 3.96e-3, a1 = -0.35, phi2 = 0.79e-3,
    a2 = -1.47 [Mpc^-3, h = 0.7] (data/abundance_matching/, verified against the paper)
  * halo mass function: Sheth-Tormen from the benchmark's Planck-like cosmology
  * gas: Papastergis et al. 2012, log(M_HI/M*) = -0.43 log M* + 3.75,
    M_b = M* + 1.4 M_HI.
The frozen cSIDM model is then evaluated on all SPARC galaxies with this
connection replacing its fitted one (A, B, floor), so nothing about galaxies
is fitted except the cSIDM cross-section, M/L and concentration offset."""
import sys, importlib.util
sys.path.insert(0, '.')
import numpy as np
from scipy.optimize import minimize_scalar
from ubench.cosmology import Cosmology
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G

import csv
def _row(path):
    return {k: float(v) for k, v in next(csv.DictReader(l for l in open(path) if not l.startswith("#"))).items()}
SMF = _row("data/abundance_matching/baldry2012_gsmf.csv")
GAS = _row("data/abundance_matching/papastergis2012_gas.csv")
lms = np.linspace(5.0, 12.5, 600)
Ms = 10 ** lms
x = Ms / 10 ** SMF["log10_Mstar_char"]
phi = np.log(10) * np.exp(-x) * (SMF["phi1"] * x ** (SMF["alpha1"] + 1) + SMF["phi2"] * x ** (SMF["alpha2"] + 1))   # per dex
n_gal = np.cumsum((phi * (lms[1] - lms[0]))[::-1])[::-1]                                     # n(>M*)

c = Cosmology()
lmh = np.linspace(7.0, 15.5, 700)
dndlnM = c.halo_mass_function(10 ** lmh, 0.0)
n_halo = np.cumsum((dndlnM * np.log(10) * (lmh[1] - lmh[0]))[::-1])[::-1]

def mstar_of_mh(mh):
    n = np.interp(np.log10(mh), lmh, np.log10(n_halo))
    return 10 ** np.interp(-n, -np.log10(n_gal), lms)

def mb_of_mstar(ms):
    fhi = 10 ** (GAS["slope"] * np.log10(ms) + GAS["intercept"])
    return ms * (1 + GAS["helium_factor"] * fhi)

grid_mh = np.logspace(9, 15, 400)
grid_mb = mb_of_mstar(mstar_of_mh(grid_mh))
def mh_of_mb(mb):
    return float(10 ** np.interp(np.log10(mb), np.log10(grid_mb), np.log10(grid_mh)))

print("Abundance-matched connection (z = 0, centrals):")
for lmb in (8, 9, 10, 10.5, 11, 11.5):
    print(f"   M_b = 1e{lmb:<4}  ->  M_h = {mh_of_mb(10**lmb):.2e}   (fitted cSIDM connection: {10**11.60*(10**lmb/1e10)**0.884 + 10**9.431:.2e})")

spec = importlib.util.spec_from_file_location("csidm", "submissions/solver_claude/my_model.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
GALS = load_sparc()

def evaluate(halo_fn, label):
    m.halo_mass = halo_fn
    M = m.CollapsingSIDM()
    gb, gm, go, zs = [], [], [], []
    for g in GALS:
        v = M.predict(g, g.R, {}, {"ups_scale": 1.0})
        vb2 = v_baryon_sq(g, m.UPS_DISK, 1.4 * m.UPS_DISK)
        ok = vb2 > 0
        gb += list(vb2[ok] / g.R[ok]); gm += list(v[ok] ** 2 / g.R[ok]); go += list(g.Vobs[ok] ** 2 / g.R[ok])
        zs += list((v - g.Vobs) / g.sigma)
    conv = 1e6 / 3.0857e19
    gb, gm, go = np.array(gb) * conv, np.array(gm) * conv, np.array(go) * conv
    k = gb > 1e-13
    def fit(y):
        rar = lambda a0: gb[k] / (-np.expm1(-np.sqrt(gb[k] / a0)))
        return 10 ** minimize_scalar(lambda la: np.mean((np.log10(y[k]) - np.log10(rar(10 ** la))) ** 2),
                                     bounds=(-11.5, -8.5), method="bounded").x
    zs = np.array(zs)
    print(f"{label:38s} emergent a0 = {fit(gm):.2e} m/s^2 | rotation curves within 2 sigma {100*np.mean(abs(zs)<2):.1f}% (all 165 galaxies)")
    return fit(go)

orig = m.halo_mass
a0_data = evaluate(orig, "fitted connection (frozen model)")
evaluate(mh_of_mb, "abundance-matched link (concentration offset still fitted)")
print(f"\ndata (same galaxies, same M/L): a0 = {a0_data:.2e};  MOND standard 1.2e-10;  cH0/2pi = {2.998e8*67.4e3/3.0857e22/2/np.pi:.2e}")

print("\nSensitivity of the emergent a0 to the fitted cSIDM constants (abundance-matched connection):")
base_s, base_c = m.SIGMA0, m.DLOGC
for lab, s0, dc in [("sigma0 / 10", base_s / 10, base_c), ("sigma0 x 3", base_s * 3, base_c),
                    ("concentration -0.2 dex", base_s, base_c - 0.2), ("concentration +0.2 dex", base_s, base_c + 0.2),
                    ("standard LCDM concentration (offset 0)", base_s, 0.0)]:
    m.SIGMA0, m.DLOGC = s0, dc
    evaluate(mh_of_mb, f"  {lab}")
m.SIGMA0, m.DLOGC = base_s, base_c
