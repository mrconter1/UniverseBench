"""
Test 03: galaxy cluster masses from X-ray gas.

Data: 13 relaxed galaxy clusters observed with Chandra (Vikhlinin et al.
2006, ApJ 640, 691; z from Table 1, the rest from Table 3). X-rays measure
the hot gas, and its temperature profile gives the total mass holding it in
place (hydrostatic equilibrium). Up to two radii per cluster, r2500 and r500
(mean density 2500x / 500x the critical density): 23 measurements.

Input: visible (baryonic) mass inside the radius = hot gas + stars; stars
are estimated as 15% of the gas (typical for clusters, Gonzalez+ 2013).
Uncertainty: quoted error plus a 15% floor for hydrostatic bias and the
stellar estimate.

Model method: predict_cluster_mass(c, gp) -> total mass [Msun] inside
c.r_kpc, given c.M_bar [Msun], c.z, c.delta (2500 or 500).
"""
import csv
import math

import numpy as np

from ..core import Test, data_path, is_hidden

STAR_TO_GAS = 0.15
FLOOR = 0.15
_H, _OM = 0.71, 0.3                      # cosmology used by the paper
_RHO_C0 = 2.775e11 * _H**2 / 1e9         # Msun / kpc^3


class ClusterPoint:
    def __init__(self, name, z, delta, r_kpc, M_tot, M_tot_err, f_gas):
        self.name, self.z, self.delta = name, z, delta
        self.r_kpc = r_kpc
        self.M_tot, self.sigma = M_tot, math.hypot(M_tot_err, FLOOR * M_tot)
        self.M_gas = f_gas * M_tot
        self.M_star = STAR_TO_GAS * self.M_gas
        self.M_bar = self.M_gas + self.M_star


def _radius(M, z, delta):
    rho = _RHO_C0 * (_OM * (1 + z) ** 3 + 1 - _OM)
    return (3 * M / (4 * math.pi * delta * rho)) ** (1 / 3)


def load_clusters(path=None):
    path = path or data_path("clusters", "vikhlinin2006.csv")
    rows = csv.DictReader(l for l in open(path) if not l.startswith("#"))
    pts = []
    for r in rows:
        z = float(r["z"])
        if r["M2500"]:
            M = float(r["M2500"]) * 1e14
            pts.append(ClusterPoint(r["name"], z, 2500, _radius(M, z, 2500), M,
                                    float(r["M2500_err"]) * 1e14, float(r["fg2500"])))
        if r["M500"]:
            M = float(r["M500"]) * 1e14
            pts.append(ClusterPoint(r["name"], z, 500, float(r["r500"]), M,
                                    float(r["M500_err"]) * 1e14, float(r["fg500"])))
    return pts


class ClusterMasses(Test):
    ID = "03_clusters"
    LABEL = "Galaxy clusters"

    def load(self):
        self.points = load_clusters()

    def export_public(self, dest_data, salt):
        import os
        os.makedirs(os.path.join(dest_data, "clusters"), exist_ok=True)
        with open(data_path("clusters", "vikhlinin2006.csv")) as fi, \
                open(os.path.join(dest_data, "clusters", "vikhlinin2006.csv"), "w") as fo:
            for line in fi:
                name = line.split(",")[0]
                if line.startswith(("#", "name,")) or not is_hidden(name, salt):
                    fo.write(line)

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_cluster_mass"):
            return None
        pts = self.points if part == "all" else \
            [c for c in self.points if is_hidden(c.name, salt) == (part == "hidden")]
        return np.array([(model.predict_cluster_mass(c, gp) - c.M_tot) / c.sigma for c in pts])
