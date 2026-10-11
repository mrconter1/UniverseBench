"""
Test C: galaxy cluster masses.

Data: 13 relaxed galaxy clusters observed with Chandra (Vikhlinin et al.
2006). X-rays measure the hot gas, and the gas temperature profile gives
the total mass holding it in place (hydrostatic equilibrium). For each
cluster we use up to two radii, r2500 and r500 (where the mean density is
2500x and 500x the cosmic critical density).

A model gets the visible (baryonic) mass inside the radius and must predict
the total mass there. Baryons = hot gas + stars; stars are estimated as 15%
of the gas mass (typical for clusters; Gonzalez et al. 2013).

Uncertainty per point: the quoted error plus a 15% floor for hydrostatic
mass bias and the stellar estimate.

Models implement:  predict_cluster_mass(c, gp) -> total mass [Msun]
where c has c.r_kpc (radius), c.M_bar (baryons inside it, Msun), c.z.
"""
import csv
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
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
        self.M_bar = self.M_gas * (1 + STAR_TO_GAS)


def _radius(M, z, delta):
    rho = _RHO_C0 * (_OM * (1 + z) ** 3 + 1 - _OM)
    return (3 * M / (4 * math.pi * delta * rho)) ** (1 / 3)


def load_clusters(path=None):
    path = path or os.path.join(HERE, "data", "clusters", "vikhlinin2006.csv")
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


def residuals(model, pts, gp):
    if not hasattr(model, "predict_cluster_mass"):
        return None
    return np.array([(model.predict_cluster_mass(c, gp) - c.M_tot) / c.sigma for c in pts])
