"""Planck 2018 plik-lite likelihood (binned TT/TE/EE, foregrounds marginalised)
plus two low-ell TT bins and a Gaussian tau prior. An independent
implementation of the public plik-lite logic (Planck 2018 V, arXiv:1907.12875;
Prince & Dunkley 2019). Returns WHITENED residuals so the benchmark can score
the fraction of decorrelated data points within 2 sigma, like the BAO test."""
import os
import numpy as np
from scipy.io import FortranFile
from scipy.linalg import cholesky, solve_triangular

from .core import data_path

PLMIN, NTT, NTE, NEE = 30, 215, 199, 199
TAU_OBS, TAU_ERR = 0.0544, 0.0073


class PlanckLite:
    def __init__(self):
        d = data_path("cmb_planck_lite", "planck2018_plik_lite")
        lo = data_path("cmb_planck_lite", "planck2018_low_ell")
        _, X, _ = np.genfromtxt(os.path.join(d, "cl_cmb_plik_v22.dat"), unpack=True)
        self.blmin = np.loadtxt(os.path.join(d, "blmin.dat")).astype(int)
        self.blmax = np.loadtxt(os.path.join(d, "blmax.dat")).astype(int)
        self.w = np.loadtxt(os.path.join(d, "bweight.dat"))
        n = NTT + NTE + NEE
        cov = FortranFile(os.path.join(d, "c_matrix_plik_v22.dat"), "r").read_reals(dtype=float).reshape(n, n)
        cov = np.triu(cov) + np.triu(cov, 1).T if np.allclose(np.tril(cov, -1), 0) else np.tril(cov) + np.tril(cov, -1).T
        self.X_hi = X
        self.L = cholesky(cov, lower=True)
        _, self.X_lo, self.sig_lo = np.genfromtxt(os.path.join(lo, "CTT_bin_low_ell_2018.dat"), unpack=True)
        self.blmin_lo = np.loadtxt(os.path.join(lo, "blmin_low_ell.dat")).astype(int)
        self.blmax_lo = np.loadtxt(os.path.join(lo, "blmax_low_ell.dat")).astype(int)
        self.w_lo = np.loadtxt(os.path.join(lo, "bweight_low_ell.dat"))
        self.plmin_lo = 2

    def _bin(self, cl, blmin, blmax, w, plmin, nb):
        return np.array([np.sum(cl[blmin[i] + plmin:blmax[i] + plmin + 1] * w[blmin[i]:blmax[i] + 1])
                         for i in range(nb)])

    def residuals(self, Dltt, Dlte, Dlee, tau):
        """D_l arrays indexed from l = 0 (muK^2). Returns whitened residuals:
        2 low-l TT bins, 613 decorrelated high-l bins, 1 tau point."""
        l = np.arange(len(Dltt))
        fac = np.where(l > 0, l * (l + 1) / (2 * np.pi), 1.0)
        tt, te, ee = Dltt / fac, Dlte / fac, Dlee / fac
        hi = np.concatenate([self._bin(tt, self.blmin, self.blmax, self.w, PLMIN, NTT),
                             self._bin(te, self.blmin, self.blmax, self.w, PLMIN, NTE),
                             self._bin(ee, self.blmin, self.blmax, self.w, PLMIN, NEE)])
        lo = self._bin(tt, self.blmin_lo, self.blmax_lo, self.w_lo, self.plmin_lo, len(self.X_lo))
        z_hi = solve_triangular(self.L, hi - self.X_hi, lower=True)
        z_lo = (lo - self.X_lo) / self.sig_lo
        return np.concatenate([z_lo, z_hi, [(tau - TAU_OBS) / TAU_ERR]])
