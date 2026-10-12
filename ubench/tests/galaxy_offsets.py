"""
Test 01b: whole-galaxy offsets with each galaxy's own error budget (SPARC).

Why: in test 01 every point is independent with a 5% floor, but the dominant
errors (distance, inclination, stellar mass-to-light) move a galaxy's whole
rotation curve together. So test 01 barely discriminates (every reasonable
model scores 55-58). Here each galaxy is one data point:

    offset = median over radii of log10(V_obs^2 / V_model^2)
    error^2 = stat^2 + (e_D/D)^2 + (2 cot(i) e_i)^2 [/ln 10] + (0.1 dex M/L x stellar share x 0.75)^2

with distances, inclinations and their errors from SPARC Table 1
(data/sparc_errors/). z = offset / error, scored within 2 sigma. Quality-1
galaxies with inclination >= 30 deg (Lelli+ 2016 flag; quality-2 galaxies
show ~0.17 dex of extra scatter from non-circular motions and are left out).
No per-galaxy freedom (default local parameters). A model with the right
mean relation and no intrinsic scatter scores ~93 here on SPARC (MOND with
a0 = 1.2e-10: 92.5%); halo models that shift galaxies systematically by
mass, or scatter them, score lower.

Model method: predict(g, R, gp, lp), as in test 01.
"""
import csv
import os

import numpy as np

from ..core import Test, data_path, default_locals, is_hidden
from ..models.common import v_baryon_sq
from .rotation import load_sparc

U_ERR = 0.5      # M/L used only to estimate each galaxy's stellar share for the error budget


class GalaxyOffsets(Test):
    ID = "01b_galaxy_offsets"
    LABEL = "Rotation curves: whole-galaxy offsets"
    DATA_DIRS = ("sparc", "sparc_errors")

    def load(self):
        T = {r["name"]: r for r in csv.DictReader(l for l in open(data_path("sparc_errors", "sparc_errors.csv"))
                                                   if not l.startswith("#"))}
        self.rows = []
        for g in load_sparc():
            t = T.get(g.name)
            if t is None or int(t["Q"]) != 1 or float(t["inc"]) < 30:
                continue
            vb2 = v_baryon_sq(g, U_ERR, 1.4 * U_ERR)
            ok = (vb2 > 0) & (g.Vobs > 0)
            if ok.sum() < 5:
                continue
            vstar2 = U_ERR * g.Vdisk * np.abs(g.Vdisk) + 1.4 * U_ERR * g.Vbul * np.abs(g.Vbul)
            fstar = np.median(np.clip(vstar2[ok] / vb2[ok], 0, 1))
            D, eD, inc, einc = (float(t[k]) for k in ("D", "e_D", "inc", "e_inc"))
            stat = np.median(2 * g.sigma[ok] / g.Vobs[ok] / np.log(10)) / np.sqrt(ok.sum())
            err = np.sqrt(stat ** 2 + (eD / D / np.log(10)) ** 2
                          + (2 / np.tan(np.radians(inc)) * np.radians(einc) / np.log(10)) ** 2
                          + (0.1 * fstar * 0.75) ** 2)
            self.rows.append((g, ok, err))

    def export_public(self, dest_data, salt):
        import shutil
        os.makedirs(os.path.join(dest_data, "sparc_errors"), exist_ok=True)
        names = {g.name for g, _, _ in self.rows if not is_hidden(g.name, salt)}
        src = data_path("sparc_errors", "sparc_errors.csv")
        with open(src) as f, open(os.path.join(dest_data, "sparc_errors", "sparc_errors.csv"), "w") as o:
            for line in f:
                if line.startswith("#") or line.startswith("name") or line.split(",")[0] in names:
                    o.write(line)
        os.makedirs(os.path.join(dest_data, "sparc"), exist_ok=True)
        for n in names:
            shutil.copy(data_path("sparc", f"{n}_rotmod.dat"), os.path.join(dest_data, "sparc"))

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict"):
            return None
        lp = default_locals(model)
        z = []
        for g, ok, err in self.rows:
            if part != "all" and is_hidden(g.name, salt) != (part == "hidden"):
                continue
            v = np.asarray(model.predict(g, g.R, gp, lp), dtype=float)
            if np.any(v[ok] <= 0) or not np.all(np.isfinite(v[ok])):
                z.append(1e3)
                continue
            z.append(np.median(np.log10(g.Vobs[ok] ** 2 / v[ok] ** 2)) / err)
        return np.array(z)
