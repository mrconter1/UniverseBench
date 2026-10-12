"""Joint fit: rotation curves (public SPARC) + dwarfs (public LVDB), each
weighted as a mean z^2, as the benchmark does.
usage: python dev/joint.py TAG [eval]"""
import sys, json, os, time
sys.path[:0] = ['.', 'dev']
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import rc
from rc import GALS
from ubench.physics import G, nfw_concentration, nfw_r200
from ubench.tests.dwarfs import load_dwarfs
from sidm import SIDMHalo, collapse_tau, nfw_m, sigma_eff, nfw_vmax

DW = load_dwarfs()
import csv as _csv
LENS = [{k: float(v) for k, v in r.items()} for r in _csv.DictReader(
    l for l in open("data/galaxy_lensing/mistele2024_btfr.csv") if not l.startswith("#"))]


def v_lens(M, R, p):
    mh = halo_mass(M, p)
    c = nfw_concentration(mh) * 10 ** p["dlogc"]
    return np.sqrt(G * (M + nfw_m(np.asarray(R, dtype=float), mh, c, nfw_r200(mh))) / np.asarray(R))


def z_lens(p):
    z = []
    for r in LENS:
        R = np.logspace(np.log10(r["R_min_kpc"]), np.log10(r["R_max_kpc"]), 25)
        M = 10 ** r["log10_Mb"]
        v = np.mean(v_lens(M, R, p)); vu = np.mean(v_lens(M * 10 ** 0.1, R, p))
        z.append((v - r["Vflat"]) / np.hypot(r["Vflat_err"], vu - v))
    return np.array(z)
KEYS = ["ups", "A", "B", "dlogc", "lsig", "lw", "lfloor", "lgx"]
BOUNDS = [(0.2, 1.0), (10.0, 13.0), (0.3, 1.5), (-0.8, 0.5), (-1.0, 3.5), (1.0, 3.0), (7.0, 10.5), (0.0, 0.6)]
X0 = [0.416, 11.547, 0.889, -0.074, 2.223, 2.028, 9.151, 0.2]


def lum_total(g):
    """Disc + 1.4 x bulge luminosity [Lsun] including an exponential tail beyond
    the last measured radius (scale length from the outer half of the profile)."""
    R, sb = g.R, g.SBdisk
    Ld, Lb = g.luminosity()
    tail = 0.0
    k = (sb > 0) & (R >= R[len(R) // 2])
    if k.sum() >= 3:
        slope = np.polyfit(R[k], np.log(sb[k]), 1)[0]
        if slope < -1e-3:
            h = -1 / slope * 1000.0                         # pc
            Rl = R[-1] * 1000.0
            S_last = np.exp(np.polyval(np.polyfit(R[k], np.log(sb[k]), 1), R[-1]))
            tail = 2 * np.pi * S_last * h * (Rl + h)
    return Ld + tail, Lb


_LT = {}


def baryon_mass(g, p):
    if g.name not in _LT:
        _LT[g.name] = lum_total(g)
    Ld, Lb = _LT[g.name]
    u = p["ups"]
    return u * Ld + 1.4 * u * Lb + 10 ** p["lgx"] * rc.gas_mass(g)


def halo_mass(mb, p):
    return 10 ** p["A"] * (mb / 1e10) ** p["B"] + 10 ** p["lfloor"]


rc.halo_from_baryons = None
WARM_D = {}


def v_rc(g, p):
    u = p["ups"]
    vb2 = rc.v_baryon_sq(g, u, 1.4 * u)
    mh = halo_mass(baryon_mass(g, p), p)
    c = nfw_concentration(mh) * 10 ** p["dlogc"]
    h, ls = make_halo(mh, c, rc.mb_profile(g, u), p, rc._WARM, g.name)
    return np.sqrt(np.clip(vb2, 0, None) + rc.dm_mass_blend(h, g.R, ls, 1e5) * G / g.R)


def make_halo(mh, c, mb_of_r, p, warm, key):
    """One effective cross-section per halo (Yang & Yu 2022 weighting at the
    halo's own Vmax), used for both the thermalisation radius and the
    collapse clock."""
    r200 = nfw_r200(mh)
    ls = np.log10(sigma_eff(p["lsig"], 10 ** p["lw"], nfw_vmax(mh, c, r200)))
    h = SIDMHalo(mh, c, r200, mb_of_r, ls, 1e5, guess=warm.get(key))
    if h.ok:
        warm[key] = (h.sigma0, np.log(h.rho0))
    return h, ls


def plummer(M, rh3):
    a = rh3 / 1.305
    return lambda r: M * r ** 3 / (r ** 2 + a ** 2) ** 1.5


def sigma_dwarf(d, p, mstar=None):
    M = d.M_bar if mstar is None else mstar
    mh = halo_mass(M, p)
    c = nfw_concentration(mh) * 10 ** p["dlogc"]
    key = (d.name, mstar is None)
    h, ls = make_halo(mh, c, plummer(M, d.r_half), p, WARM_D, key)
    mdm = float(rc.dm_mass_blend(h, np.array([d.r_half]), ls, 1e5)[0])
    return np.sqrt(G * (M / 2 + mdm) / (3 * d.r_half))


def z_dwarf(d, p):
    s = sigma_dwarf(d, p)
    if d.is_limit:
        return 0.0 if s <= d.sigma_ul else 2 * s / d.sigma_ul
    s_up = sigma_dwarf(d, p, d.M_bar * 10 ** d.log_Mb_err)
    e = d.sigma_ep if s > d.sigma_obs else d.sigma_em
    return (s - d.sigma_obs) / np.hypot(e, s_up - s)


def job(args):
    kind, i, x = args
    p = dict(zip(KEYS, np.clip(x, [b[0] for b in BOUNDS], [b[1] for b in BOUNDS])))
    try:
        if kind == "g":
            g = GALS[i]
            z = (v_rc(g, p) - g.Vobs) / g.sigma
            return kind, z if np.all(np.isfinite(z)) else np.full(len(g.R), 30.0)
        z = z_dwarf(DW[i], p)
        return kind, np.array([z if np.isfinite(z) else 30.0])
    except Exception:
        return kind, (np.full(len(GALS[i].R), 30.0) if kind == "g" else np.array([30.0]))


def evaluate(pool, x, gal_idx=None, dw_idx=None):
    gal_idx = range(len(GALS)) if gal_idx is None else gal_idx
    dw_idx = range(len(DW)) if dw_idx is None else dw_idx
    tasks = [("g", i, x) for i in gal_idx] + [("d", i, x) for i in dw_idx]
    out = pool.map(job, tasks, chunksize=4)
    zg = np.concatenate([z for k, z in out if k == "g"])
    zd = np.concatenate([z for k, z in out if k == "d"])
    return zg, zd


if __name__ == "__main__":
    tag = sys.argv[1]
    ck = f"dev/ck_{tag}.json"
    state = json.load(open(ck)) if os.path.exists(ck) else {"x": X0, "f": None, "done": False, "nfev": 0}
    t0 = time.time()
    with Pool(2) as pool:
        if "eval" not in sys.argv and not state["done"]:
            def obj(x):
                zg, zd = evaluate(pool, x)
                zl = z_lens(dict(zip(KEYS, x)))
                v = float(np.mean(zg ** 2) + np.mean(zd ** 2) + np.mean(zl ** 2))
                state["nfev"] += 1
                if state["f"] is None or v < state["f"]:
                    state["f"], state["x"] = v, list(map(float, x))
                    json.dump(state, open(ck, "w"))
                if time.time() - t0 > float(os.environ.get("TLIM", 540)):
                    raise TimeoutError
                return v
            x0 = np.array(state["x"])
            simplex = [x0] + [x0 + np.eye(len(x0))[i] * 0.08 * (BOUNDS[i][1] - BOUNDS[i][0])
                              * (1 if x0[i] + 0.08 * (BOUNDS[i][1] - BOUNDS[i][0]) <= BOUNDS[i][1] else -1)
                              for i in range(len(x0))]
            try:
                minimize(obj, x0, method="Nelder-Mead",
                         options={"maxfev": 500, "xatol": 1e-3, "fatol": 1e-3, "initial_simplex": np.array(simplex)})
                state["done"] = True
                json.dump(state, open(ck, "w"))
            except TimeoutError:
                print("checkpointed")
        zg, zd = evaluate(pool, state["x"])
        zl = z_lens(dict(zip(KEYS, state["x"])))
        print(f"{tag}: RC within2s {100*np.mean(abs(zg)<2):.1f}% (z2 {np.mean(zg**2):.2f}) | dwarfs within2s {100*np.mean(abs(zd)<2):.1f}% (z2 {np.mean(zd**2):.2f}) | lensing {100*np.mean(abs(zl)<2):.0f}% z={np.round(zl,1)}")
        print({k: round(v, 3) for k, v in zip(KEYS, state["x"])}, "nfev", state["nfev"], "done", state["done"])
