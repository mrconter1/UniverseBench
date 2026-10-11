"""
Stage 3 cosmology tests, part 1: the expansion history and the early universe.

Models supply a cosmology via model.cosmology(gp), returning an object with
the interface of ubench.cosmology.Cosmology (or None: no prediction).
The Big Bang nucleosynthesis test also accepts model.predict_bbn(gp).

15 CMB: compressed Planck 2018 constraints: the acoustic angle theta*, the
   shift parameter R, the baryon density, and the primordial spectrum (n_s,
   A_s). Fitted and scored in-sample.
16 BBN: primordial helium, deuterium and lithium (PDG 2025). Lithium is the
   long-standing "lithium problem": standard physics overpredicts it ~3x.
17 BAO: DESI DR2 baryon acoustic oscillation distances (13 numbers with full
   covariance; residuals are decorrelated before scoring).
18 Supernovae: 1580 Pantheon+ type Ia supernovae (z > 0.01). Only relative
   distances matter: the supernova absolute magnitude is a nuisance fitted on
   public supernovae, then applied to hidden ones.
"""
import csv

import numpy as np

from ..core import Test, data_path, is_hidden


def _rows(*path):
    return list(csv.DictReader(l for l in open(data_path(*path)) if not l.startswith("#")))


def _cosmo(model, gp):
    if not hasattr(model, "cosmology"):
        return None
    try:
        return model.cosmology(gp)
    except (ValueError, ZeroDivisionError, FloatingPointError):
        return "bad"


class _Cosmo(Test):
    STAGE = 3
    HIDDEN_SPLIT = False


class CMB(_Cosmo):
    ID, LABEL = "15_cmb", "Microwave background"
    DATA_DIRS = ("cmb",)

    def load(self):
        self.rows = _rows("cmb", "planck2018.csv")

    def residuals(self, model, gp, part, salt):
        c = _cosmo(model, gp)
        if c is None:
            return None
        if c == "bad":
            return np.full(len(self.rows), 1e6)
        p = c.cmb_distance_priors()
        pred = {"theta100": 100 * c.theta_star(), "R": p["R"], "omega_b": c.ombh2,
                "n_s": c.ns, "ln10As": c.ln10As}
        return np.array([(pred[r["quantity"]] - float(r["value"]))
                         / np.hypot(float(r["sigma"]), float(r["engine_floor"])) for r in self.rows])


class BBN(_Cosmo):
    ID, LABEL = "16_bbn", "Big Bang nucleosynthesis"
    DATA_DIRS = ("bbn",)

    def load(self):
        self.rows = _rows("bbn", "pdg2025.csv")

    def residuals(self, model, gp, part, salt):
        if hasattr(model, "predict_bbn"):
            y = model.predict_bbn(gp)
        else:
            c = _cosmo(model, gp)
            if c is None:
                return None
            if c == "bad":
                return np.full(len(self.rows), 1e6)
            y = c.bbn()
        return np.array([(y[r["quantity"]] - float(r["value"]))
                         / np.hypot(float(r["sigma"]), float(r["theory_floor"])) for r in self.rows])


class BAO(_Cosmo):
    ID, LABEL = "17_bao", "Baryon acoustic oscillations"
    DATA_DIRS = ("bao",)

    def load(self):
        rows = [l.split() for l in open(data_path("bao", "desi_dr2_mean.txt")) if not l.startswith("#")]
        self.z = np.array([float(r[0]) for r in rows])
        self.val = np.array([float(r[1]) for r in rows])
        self.kind = [r[2] for r in rows]
        cov = np.loadtxt(data_path("bao", "desi_dr2_cov.txt"))
        self.Linv = np.linalg.inv(np.linalg.cholesky(cov))

    def residuals(self, model, gp, part, salt):
        c = _cosmo(model, gp)
        if c is None:
            return None
        if c == "bad":
            return np.full(len(self.val), 1e6)
        rd = c.r_d()
        f = {"DV_over_rs": c.D_V, "DM_over_rs": c.D_M, "DH_over_rs": c.D_H}
        pred = np.array([float(f[k](z)) / rd for z, k in zip(self.z, self.kind)])
        return self.Linv @ (pred - self.val)          # decorrelated residuals


class Supernovae(_Cosmo):
    ID, LABEL = "18_supernovae", "Supernova distances"
    HIDDEN_SPLIT = True

    def load(self):
        self.rows = _rows("supernovae", "pantheon_plus.csv")
        self.cid = np.array([r["CID"] for r in self.rows])
        self.zhd = np.array([float(r["zHD"]) for r in self.rows])
        self.zhel = np.array([float(r["zHEL"]) for r in self.rows])
        self.m = np.array([float(r["m_b_corr"]) for r in self.rows])
        self.err = np.array([float(r["err"]) for r in self.rows])
        self._hidden_cache = {}

    def hidden_mask(self, salt):
        if salt not in self._hidden_cache:
            self._hidden_cache[salt] = np.array([is_hidden("SN" + c, salt) for c in self.cid])
        return self._hidden_cache[salt]

    def export_public(self, dest_data, salt):
        import os
        os.makedirs(os.path.join(dest_data, "supernovae"), exist_ok=True)
        with open(data_path("supernovae", "pantheon_plus.csv")) as fi, \
                open(os.path.join(dest_data, "supernovae", "pantheon_plus.csv"), "w") as fo:
            for line in fi:
                if line.startswith(("#", "CID,")) or not is_hidden("SN" + line.split(",")[0], salt):
                    fo.write(line)

    def residuals(self, model, gp, part, salt):
        c = _cosmo(model, gp)
        if c is None:
            return None
        hid = self.hidden_mask(salt)
        sel = {"all": np.ones_like(hid), "public": ~hid, "hidden": hid}[part]
        if c == "bad":
            return np.full(int(sel.sum()), 1e6)
        mu = c.distance_modulus(self.zhd, self.zhel)
        # absolute magnitude: weighted mean offset on the public supernovae
        w = 1 / self.err**2
        ref = ~hid if part != "all" else np.ones_like(hid)
        M = np.sum(w[ref] * (self.m - mu)[ref]) / np.sum(w[ref])
        return ((mu + M - self.m) / self.err)[sel]
