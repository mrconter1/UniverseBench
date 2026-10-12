"""Candidate I: does the cSIDM 'age clock' make the effective a0 rise with
redshift? The frozen cSIDM model (submissions/solver_claude) is evaluated on
all SPARC baryon distributions with halos placed at z = 0 or z = 1:
  clock   : halo age 10 Gyr -> 4.3 Gyr (universe age ratio 5.9/13.8)
  density : r200 at rho_crit(z) = rho_crit0 E(z)^2, concentration from
            Dutton & Maccio 2014 at z (a = 0.520 + 0.385 exp(-0.617 z^1.21),
            b = -0.101 + 0.026 z)
The galaxy-halo connection is kept fixed (an assumption).
The effective a0 is fitted to the MODEL's accelerations with the RAR form."""
import sys, importlib.util
sys.path.insert(0, '.')
import numpy as np
from scipy.optimize import minimize_scalar
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench import physics

spec = importlib.util.spec_from_file_location("csidm", "submissions/solver_claude/my_model.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
GALS = load_sparc()
OM = 0.31
E = lambda z: np.sqrt(OM * (1 + z) ** 3 + 1 - OM)
H0K = physics.H0_KPC

def r200_z(M, z):
    return (physics.G * M / (100 * (H0K * E(z)) ** 2)) ** (1 / 3)

def c_z(M, z, h=0.7):
    a = 0.520 + (0.905 - 0.520) * np.exp(-0.617 * z ** 1.21)
    b = -0.101 + 0.026 * z
    return 10 ** (a + b * np.log10(M * h / 1e12)) * 10 ** m.DLOGC

def setup(z_density, t_age):
    m.T_AGE = t_age
    m.nfw_r200 = lambda M: r200_z(M, z_density)
    m.concentration = lambda M: c_z(M, z_density)

def effective_a0():
    M = m.CollapsingSIDM()
    gb, gm = [], []
    for g in GALS:
        v = M.predict(g, g.R, {}, {"ups_scale": 1.0})
        vb2 = v_baryon_sq(g, m.UPS_DISK, 1.4 * m.UPS_DISK)
        ok = vb2 > 0
        gb += list(vb2[ok] / g.R[ok]); gm += list(v[ok] ** 2 / g.R[ok])
    gb, gm = np.array(gb) * 1e6 / 3.0857e19, np.array(gm) * 1e6 / 3.0857e19
    k = gb > 1e-13
    rar = lambda a0: gb[k] / (-np.expm1(-np.sqrt(gb[k] / a0)))
    r = minimize_scalar(lambda la: np.mean((np.log10(gm[k]) - np.log10(rar(10 ** la))) ** 2),
                        bounds=(-11.5, -8.5), method="bounded")
    lo = gb[k] < 1e-11
    return 10 ** r.x, np.median(np.log10(gm[k][lo] / gb[k][lo]))

out = {}
for lab, zd, t in [("z=0 (reference)", 0.0, 10.0), ("clock only (age 4.3 Gyr)", 0.0, 4.3),
                   ("density only (halos at z=1)", 1.0, 10.0), ("clock + density (z=1)", 1.0, 4.3)]:
    setup(zd, t)
    a0, boost = effective_a0()
    out[lab] = a0
    print(f"{lab:30s} effective a0 = {a0:.2e} m/s^2   ratio to z=0: {a0/out['z=0 (reference)']:.2f}   (median log g_obs/g_bar at g_bar<1e-11: {boost:.2f})")
print(f"\nfor comparison: H(z=1)/H0 = {E(1.0):.2f};  MUSE-DARK III: ~2.4-2.6 at z~1;  MOND: 1.00")
