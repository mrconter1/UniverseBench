import sys
sys.path[:0] = ['.', 'dev']
import numpy as np
from ubench.tests.dwarfs import load_dwarfs
from ubench.physics import G, nfw_concentration, nfw_r200, nfw_mass_enclosed
from sidm import SIDMHalo
D = load_dwarfs()

def plummer_m(M, rh3):
    a = rh3 / 1.305           # Plummer: 3D half-mass radius = 1.305 a
    return lambda r: M * r ** 3 / (r ** 2 + a ** 2) ** 1.5

def sig_wolf(Mhalf, rh):
    return np.sqrt(G * Mhalf / (3 * rh))

def zscore(pred, d):
    if d.is_limit:
        return 0.0 if pred <= d.sigma_ul else 2 * pred / d.sigma_ul
    e = d.sigma_ep if pred > d.sigma_obs else d.sigma_em
    return (pred - d.sigma_obs) / np.hypot(e, 0.1 * pred)

def run(lmh_fixed=None, lsig=None, w=1e5, cfac=1.0, label=""):
    zz = []
    for d in D:
        mh = 10 ** lmh_fixed
        c = nfw_concentration(mh) * cfac
        if lsig is None:
            mdm = nfw_mass_enclosed(d.r_half, mh, c)
        else:
            h = SIDMHalo(mh, c, nfw_r200(mh), plummer_m(d.M_bar, d.r_half), lsig, w)
            mdm = h.mass(d.r_half)
        zz.append(zscore(sig_wolf(d.M_bar / 2 + mdm, d.r_half), d))
    zz = np.array(zz)
    print(f"{label:45s} within2s {100*np.mean(abs(zz)<2):5.1f}%  median z {np.median(zz):+.2f}  n={len(zz)}")

if __name__ == "__main__":
    for l in (8.5, 9.0, 9.5):
        run(l, None, label=f"NFW, all halos 10^{l}")
    for l in (9.0, 9.5, 10.0):
        for ls in (0.0, 0.88):
            run(l, ls, label=f"SIDM {10**ls:.1f} cm2/g, halos 10^{l}")
