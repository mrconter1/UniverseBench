"""Idea 7: a core of degenerate keV fermions in place of Sgr A*'s black hole
(Argüelles, Ruffini et al.; Becerra-Vergara+ 2020; Crespi+ 2025 compare
56 and 300 keV). For each fermion mass, solve the TOV equation for a fully
degenerate (T = 0) spin-1/2 Fermi gas holding Sgr A*'s mass. T = 0 is the
MOST compact configuration on the stable branch (thermal pressure only
puffs the core up), so it is the best case for the idea.
Tests: (1) S2 must see a point mass: core radius << S2 pericentre;
(2) GRAVITY's flares orbit at 6-10 r_g with a 45 +/- 15 min period, so most
of the mass must lie inside ~8 r_g; (3) a photon sphere (R < 3 r_g) is
needed for a black-hole-like shadow."""
import csv, pathlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = {r["quantity"]: float(r["value"]) for r in
     csv.DictReader(l for l in open(ROOT / "data/sgr_a/sgra_constraints.csv") if not l.startswith("#"))}
G, c, hbar, Msun, AU, keV = 6.674e-11, 2.998e8, 1.0546e-34, 1.989e30, 1.496e11, 1.602e-16
M_SGR = D["mass_Msun"] * Msun
rg = G * M_SGR / c**2


def _eps(x):
    return np.where(x < 1e-3, 8 / 3 * x**3 * (1 + 0.3 * x**2),
                    x * (2 * x**2 + 1) * np.sqrt(1 + x**2) - np.arcsinh(x))


def _P(x):
    return np.where(x < 1e-3, 8 / 15 * x**5 * (1 - 5 / 14 * x**2),
                    (x * (2 * x**2 - 3) * np.sqrt(1 + x**2) + 3 * np.arcsinh(x)) / 3)


def star_dimless(xc):
    """TOV in units where eps, P are in k = m^4 c^5/(8 pi^2 hbar^3),
    r in L = c^2/sqrt(4 pi G k), m in L c^2/G. Variable: Fermi momentum x = p_F/(m c)
    (dP/dx = (8/3) x^4/sqrt(1+x^2))."""
    def rhs(r, y):
        x, mm = y
        x = max(x, 0.0)
        e, p = float(_eps(x)), float(_P(x))
        dPdr = -(e + p) * (mm + r**3 * p) / (r**2 * (1 - 2 * mm / r))
        dPdx = 8 / 3 * x**4 / np.sqrt(1 + x**2) if x > 0 else 1e-300
        return [dPdr / dPdx, r**2 * e]
    surf = lambda r, y: y[0] - 1e-7 * xc
    surf.terminal = True
    r0 = 1e-6 / max(xc, 1e-3)
    return solve_ivp(rhs, [r0, 1e12], [xc, r0**3 * float(_eps(xc)) / 3], events=surf, rtol=1e-9, atol=1e-30,
                     dense_output=True)


def units(m_keV):
    m = m_keV * keV / c**2
    k = m**4 * c**5 / (8 * np.pi**2 * hbar**3)
    L = c**2 / np.sqrt(4 * np.pi * G * k)
    return L, L * c**2 / G


LX = np.linspace(-3, 1.0, 81)
SOLS = [star_dimless(10**l) for l in LX]
MD = np.array([s.y[1, -1] for s in SOLS])        # dimensionless total masses
I_MAX = int(np.argmax(MD))


def config(m_keV):
    L, Mu = units(m_keV)
    Mmax = MD[I_MAX] * Mu
    if Mmax < M_SGR:
        return None, Mmax
    target = M_SGR / Mu
    lx = brentq(lambda l: star_dimless(10**l).y[1, -1] - target, LX[0], LX[I_MAX], xtol=1e-7)
    s = star_dimless(10**lx)
    Rd = s.t[-1]
    return (Rd * L, lambda r: M_SGR if r >= Rd * L else float(s.sol(r / L)[1]) * Mu), Mmax


print(f"Sgr A*: M = {D['mass_Msun']:.3e} Msun, r_g = {rg/AU:.4f} AU, S2 pericentre = {D['s2_pericenter_AU']:.0f} AU "
      f"= {D['s2_pericenter_AU']*AU/rg:.0f} r_g")
r_fl = 0.5 * (D["flare_radius_rg_lo"] + D["flare_radius_rg_hi"]) * rg
P_bh = 2 * np.pi * np.sqrt(r_fl**3 / (G * M_SGR)) / 60
print(f"Black hole: Keplerian period at {r_fl/rg:.0f} r_g = {P_bh:.0f} min (observed {D['flare_period_min']:.0f} "
      f"+/- {D['flare_period_err_min']:.0f})")
for m_keV in (10, 56, 100, 200, 300, 345, 360, 400):
    cfg, Mmax = config(m_keV)
    if cfg is None:
        print(f"m = {m_keV:4d} keV: maximum degenerate mass {Mmax/Msun:.2e} Msun < Sgr A*: must collapse to a black hole")
        continue
    R, Menc = cfg
    Pfl = 2 * np.pi * np.sqrt(r_fl**3 / (G * Menc(r_fl))) / 60
    z = (Pfl - D["flare_period_min"]) / D["flare_period_err_min"]
    print(f"m = {m_keV:4d} keV: core radius {R/AU:9.3f} AU = {R/rg:8.1f} r_g "
          f"({'inside' if R < D['s2_pericenter_AU']*AU else 'OUTSIDE'} S2 pericentre); "
          f"M(<8 r_g)/M = {float(Menc(r_fl))/M_SGR:.2e}; flare period {Pfl:9.0f} min ({z:+.1f} sigma); "
          f"photon sphere: {'yes' if R < 3 * rg else 'no'} (M_OV = {Mmax/Msun:.2e})")
