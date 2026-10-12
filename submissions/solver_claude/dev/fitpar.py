"""Checkpointed, 2-core parallel fit of rotation-curve models on public SPARC.
usage: python dev/fitpar.py MODEL TAG [fold k|all] """
import sys, json, os, time, hashlib
sys.path[:0] = ['.', 'dev']
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import rc
from rc import GALS, zs

MODELS = {
  # name: (function, keys, x0, bounds)
  "sidm": (lambda g, p: rc.v_sidm(g, p), ["ups", "A", "B", "dlogc", "lsig"],
           [0.464, 11.356, 0.874, -0.051, 0.883], [(0.2, 1.0), (10.0, 13.0), (0.3, 1.5), (-0.8, 0.5), (-1.0, 1.5)]),
  "vsidm": (lambda g, p: rc.v_sidm(g, p), ["ups", "A", "B", "dlogc", "lsig", "lw"],
           [0.464, 11.356, 0.874, -0.051, 2.0, 1.9], [(0.2, 1.0), (10.0, 13.0), (0.3, 1.5), (-0.8, 0.5), (-1.0, 3.5), (1.0, 3.0)]),
  "cnfw": (lambda g, p: rc.v_corenfw(g, p), ["ups", "A", "B", "dlogc", "leta"],
           [0.5, 11.4, 0.75, 0.0, 0.2], [(0.2, 1.0), (10.0, 13.0), (0.3, 1.5), (-0.8, 0.5), (-1.0, 1.0)]),
  "mond": (lambda g, p: rc.v_mond(g, p["ups"], p["a0"] * 1e-10 * 3.0857e13), ["ups", "a0"],
           [0.5, 1.2], [(0.2, 1.0), (0.6, 2.0)]),
}


def power_law_connection(p):
    # map (A, B) -> the double power law used by rc: lN, lM1 such that Mh = 10^A (Mb/1e10)^B
    # implemented by overriding rc.halo_from_baryons
    return p


def halo_pl(mb, p4):
    A, B = p4[0], p4[1]
    return 10 ** A * (mb / 1e10) ** B


rc.halo_from_baryons = lambda mb, p4: halo_pl(mb, p4)


def expand(p):
    q = dict(p)
    if "A" in q:
        q["lN"], q["lM1"], q["beta"], q["gam"] = q["A"], q["B"], 0, 0
    return q


def per_gal(args):
    name, x, model = args
    f, keys, _, bounds = MODELS[model]
    x = np.clip(x, [b[0] for b in bounds], [b[1] for b in bounds])
    p = expand(dict(zip(keys, x)))
    g = GAL[name]
    try:
        v = f(g, p)
        z = (np.asarray(v) - g.Vobs) / g.sigma
        return z if np.all(np.isfinite(z)) else np.full(len(g.R), 30.0)
    except Exception:
        return np.full(len(g.R), 30.0)


GAL = {g.name: g for g in GALS}


def fold_of(name, k=5):
    return int(hashlib.sha256(("cv" + name).encode()).hexdigest(), 16) % k


if __name__ == "__main__":
    model, tag, fold = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else "all")
    f, keys, x0, bounds = MODELS[model]
    train = [g.name for g in GALS if fold == "all" or fold_of(g.name) != int(fold)]
    test = [g.name for g in GALS if fold != "all" and fold_of(g.name) == int(fold)]
    ck = f"dev/ck_{tag}.json"
    state = json.load(open(ck)) if os.path.exists(ck) else {"x": x0, "f": None, "done": False, "nfev": 0}
    if state["done"]:
        print("already done", state); sys.exit()
    t0 = time.time()
    with Pool(2) as pool:
        def obj(x):
            zz = pool.map(per_gal, [(n, x, model) for n in train])
            v = float(np.mean(np.concatenate(zz) ** 2))
            state["nfev"] += 1
            if state["f"] is None or v < state["f"]:
                state["f"], state["x"] = v, list(map(float, x))
                json.dump(state, open(ck, "w"))
            if time.time() - t0 > float(os.environ.get("TLIM", 540)):
                raise TimeoutError
            return v
        try:
            x0 = np.array(state["x"], dtype=float)
            simplex = [x0] + [x0 + np.eye(len(x0))[i] * 0.1 * (bounds[i][1] - bounds[i][0])
                              * (1 if x0[i] + 0.1 * (bounds[i][1] - bounds[i][0]) <= bounds[i][1] else -1)
                              for i in range(len(x0))]
            r = minimize(obj, x0, method="Nelder-Mead",
                         options={"maxfev": int(os.environ.get("MAXFEV", 400)), "xatol": 1e-3, "fatol": 1e-3,
                                  "initial_simplex": np.array(simplex)})
            state["done"] = True
            json.dump(state, open(ck, "w"))
        except TimeoutError:
            print("checkpointed; rerun to continue")
        x = state["x"]
        for label, names in [("train", train), ("test", test)]:
            if names:
                zz = np.concatenate(pool.map(per_gal, [(n, x, model) for n in names]))
                print(f"{model} {tag} {label}: within2s {100*np.mean(np.abs(zz)<2):.1f}%  mean z2 {np.mean(zz**2):.2f}  n={len(names)}")
    print({k: round(v, 3) for k, v in zip(keys, x)}, "nfev", state["nfev"], "done", state["done"])
