import sys, json
sys.path[:0] = ['.', 'dev']
import numpy as np
from fitpar import MODELS, per_gal, fold_of, GALS
from multiprocessing import Pool
model = sys.argv[1]
zz = []
with Pool(2) as pool:
    for k in range(5):
        x = json.load(open(f"dev/ck_{model}_f{k}.json"))["x"]
        names = [g.name for g in GALS if fold_of(g.name) == k]
        zz += pool.map(per_gal, [(n, x, model) for n in names])
z = np.concatenate(zz)
print(f"{model} pooled 5-fold CV: within2s {100*np.mean(np.abs(z)<2):.1f}%  mean z2 {np.mean(z**2):.2f}  median|z| {np.median(np.abs(z)):.2f}")
np.save(f"dev/cvz_{model}.npy", np.array([zg.tolist() for zg in zz], dtype=object), allow_pickle=True)
