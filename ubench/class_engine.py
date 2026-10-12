"""Thin, cached wrapper around the CLASS Boltzmann code (classy).
A model describes its cosmology as a dict of CLASS input parameters
(CosmologyMixin.class_params builds one from the standard constants);
results are cached by those parameters, in memory and on disk."""
import hashlib
import json
import os
import pickle

import numpy as np

T_CMB_MUK = 2.7255e6
LMAX = 2508
_MEM = {}
CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".class_cache")


def _key(params):
    return hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()[:24]


def run_class(params):
    """Returns dict with Dl_tt, Dl_te, Dl_ee (l = 0..LMAX, muK^2), sigma8,
    or None if CLASS fails for these parameters."""
    params = {k: (round(v, 8) if isinstance(v, float) else v) for k, v in params.items()}
    k = _key(params)
    if k in _MEM:
        return _MEM[k]
    path = os.path.join(CACHE_DIR, k + ".pkl")
    if os.path.exists(path):
        _MEM[k] = pickle.load(open(path, "rb"))
        return _MEM[k]
    from classy import Class, CosmoSevereError, CosmoComputationError
    c = Class()
    full = {"output": "tCl,pCl,lCl,mPk", "lensing": "yes", "l_max_scalars": 2600}
    full.update(params)
    try:
        c.set(full)
        c.compute()
        cl = c.lensed_cl(LMAX)
        l = np.arange(LMAX + 1)
        f = l * (l + 1) / (2 * np.pi) * T_CMB_MUK ** 2
        out = {"Dl_tt": f * cl["tt"], "Dl_te": f * cl["te"], "Dl_ee": f * cl["ee"], "sigma8": c.sigma8()}
    except (CosmoSevereError, CosmoComputationError):
        out = None
    finally:
        c.struct_cleanup()
        c.empty()
    os.makedirs(CACHE_DIR, exist_ok=True)
    pickle.dump(out, open(path, "wb"))
    _MEM[k] = out
    return out
