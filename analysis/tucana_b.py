"""Pre-registered predictions for the isolated ultra-faint dwarf Tucana B
(Sand et al. 2022, ApJL 935, L17: L_V = 5(+4/-2)e4 Lsun, R_half = 80 +- 40 pc,
D = 1.4 Mpc, nearest galaxy ~500 kpc, quenched). No velocity dispersion is
published (Oct 2026). Stellar M/L_V = 1.5-2.5 for an old, metal-poor population."""
import sys, importlib.util
sys.path.insert(0, '.')
import numpy as np
from ubench.physics import G, nfw_concentration, nfw_mass_enclosed

spec = importlib.util.spec_from_file_location("csidm", "submissions/solver_claude/my_model.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
A0 = 1.2e-10 * 3.0857e13          # (km/s)^2/kpc

def mond_iso(M):
    return (4 / 81 * G * M * A0) ** 0.25

def cdm_halo(M, r3, mh):
    return np.sqrt(G * (M / 2 + nfw_mass_enclosed(r3, mh, nfw_concentration(mh))) / (3 * r3))

def csidm(M, r3):
    mh = m.halo_mass(M)
    h = m.SIDMHalo(mh, m.concentration(mh), m._plummer(M, r3))
    return float(np.sqrt(G * (M / 2 + float(h.mass(np.array([r3]))[0])) / (3 * r3)))

rows = []
for L in (3e4, 5e4, 9e4):
    for ml in (1.5, 2.5):
        for Rh in (40, 80, 120):
            M, r3 = L * ml, 4 / 3 * Rh / 1000
            rows.append((mond_iso(M), cdm_halo(M, r3, 1e9), csidm(M, r3), np.sqrt(G * M / 2 / (3 * r3))))
r = np.array(rows)
M, r3 = 5e4 * 2.0, 4 / 3 * 0.080
print(f"Central values (L_V 5e4, M/L 2, R_half 80 pc):")
print(f"  MOND (isolated, deep regime)     sigma = {mond_iso(M):.1f} km/s")
print(f"  stars alone (Newton)             sigma = {np.sqrt(G*M/2/(3*r3)):.2f} km/s")
print(f"  CDM halo 1e9 Msun (NFW)          sigma = {cdm_halo(M, r3, 1e9):.1f} km/s")
print(f"  cSIDM (frozen model, floor halo) sigma = {csidm(M, r3):.1f} km/s")
print(f"Ranges over L_V 3-9e4, M/L 1.5-2.5, R_half 40-120 pc:")
for i, lab in enumerate(["MOND", "CDM 1e9 halo", "cSIDM", "stars alone"]):
    print(f"  {lab:14s} {r[:,i].min():.1f} - {r[:,i].max():.1f} km/s")
