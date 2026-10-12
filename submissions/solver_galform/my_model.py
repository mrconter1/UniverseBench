"""Solver model 2: cold dark matter + galaxy formation ("galform").

The leading hypothesis after 12 Oct (analysis/RESULTS.md): ordinary cold,
collisionless dark matter, a cosmological constant, general relativity. The
galaxy-scale regularities come from how galaxy formation reshapes halos,
taken from hydrodynamic simulations rather than invented:

* Halo connection: M_200 = 10^A (M_b / 1e10)^B + floor (fitted, as in the
  cSIDM solver model), concentration Dutton & Maccio 2014 x 10^dlogc.
* Halo response to galaxy formation: the DC14 profile (Di Cintio et al.
  2014b, MaGICC simulations; data/halo_response/dc14.csv). Its inner slope,
  outer slope and sharpness depend only on X = log10(M_star / M_halo):
  supernova feedback carves cores at X ~ -2.7 (dwarfs, gamma ~ 0.2), halos
  stay cuspy at X < -4.1 (ultra-faints), and contract at X -> -1.3 (Milky
  Way mass, concentration x2.6). Above X = -1.3 the paper reverts to NFW;
  here X is clamped at -1.3 instead (contraction keeps increasing with
  stellar fraction; reverting would be a discontinuity).
* Everything else (clusters, Bullet Cluster, cosmology, particle tests):
  the LCDM baseline. Dark matter particle: generic cold, no couplings.
No per-galaxy halo freedom: only a per-galaxy mass-to-light correction,
like every other model.
"""
import pathlib

import numpy as np

from ubench.models.common import LOCAL_UPS, UPS, v_baryon_sq
from ubench.models.lcdm import LCDM
from ubench.physics import G, nfw_concentration, nfw_r200, dispersion_from_mass

_D = {}
for line in open(pathlib.Path(__file__).resolve().parents[2] / "data/halo_response/dc14.csv"):
    if line.startswith("#") or line.startswith("name"):
        continue
    k, v = line.strip().split(",")
    _D[k] = float(v)


def dc14_shape(X):
    X = float(np.clip(X, _D["X_min"], _D["X_max"]))
    y = 10 ** (X + _D["a1"])
    a = _D["a0"] - np.log10(y ** _D["a2"] + y ** _D["a3"])
    b = _D["b0"] + _D["b1"] * X + _D["b2"] * X ** 2
    z = 10 ** (X + _D["g1"])
    g = _D["g0"] + np.log10(z ** _D["g2"] + z)
    cfac = 1 + _D["k0"] * np.exp(_D["k1"] * (X + 4.5))
    return a, b, g, cfac


_XS = np.geomspace(1e-5, 50.0, 1500)          # r / r_s grid


class DC14Halo:
    """M(<r) and projected M(<R) for a DC14 halo normalised to M_200."""

    def __init__(self, m200, c_dm, X):
        a, b, g, cfac = dc14_shape(X)
        self.r200 = nfw_r200(m200)
        rm2 = self.r200 / (c_dm * cfac)                         # r_-2
        self.rs = rm2 / ((2 - g) / (b - 2)) ** (1 / a)
        x = _XS
        rho = x ** -g * (1 + x ** a) ** (-(b - g) / a)
        integrand = x ** 3 * rho                                 # d M / d ln x  (per 4 pi rho_s rs^3)
        inner = x[0] ** (3 - g) / (3 - g)
        cum = inner + np.concatenate([[0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(np.log(x)))])
        norm = m200 / np.interp(self.r200 / self.rs, x, cum)
        self.lx, self.lM = np.log(x), np.log(cum * norm)
        self.rho_n = rho * norm / (4 * np.pi * self.rs ** 3)    # physical density on the grid

    def mass(self, r):
        lr = np.log(np.maximum(np.asarray(r, dtype=float) / self.rs, _XS[0]))
        return np.exp(np.interp(lr, self.lx, self.lM))

    def mass_projected(self, R):
        """Mass inside a cylinder of radius R (kpc)."""
        x = _XS
        r = x * self.rs
        out = float(self.mass(R))
        sel = r > R
        rr, rho = r[sel], self.rho_n[sel]
        f = 4 * np.pi * rr ** 3 * rho * (1 - np.sqrt(1 - (R / rr) ** 2))        # per d ln r
        return out + float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(np.log(rr))))


def _gas_mass(g):
    return max(float(g.Vgas[-1] * abs(g.Vgas[-1]) * g.R[-1] / G), 0.0)


def _lum_total(g):
    """Disc and bulge luminosity with an exponential tail beyond the last radius
    (same as the cSIDM solver model)."""
    R, sb = g.R, g.SBdisk
    Ld, Lb = g.luminosity()
    k = (sb > 0) & (R >= R[len(R) // 2])
    if k.sum() >= 3:
        fit = np.polyfit(R[k], np.log(sb[k]), 1)
        if fit[0] < -1e-3:
            h, Rl = -1 / fit[0] * 1000.0, R[-1] * 1000.0
            Ld += 2 * np.pi * np.exp(np.polyval(fit, R[-1])) * h * (Rl + h)
    return Ld, Lb


class GalForm(LCDM):
    NAME = "Cold dark matter + galaxy formation (DC14 halos)"
    GLOBAL_PARAMS = {**{k: v for k, v in LCDM.GLOBAL_PARAMS.items()},
                     "halo_A": (11.6, 10.5, 13.0),        # log10 M_200 at M_b = 1e10 Msun
                     "halo_B": (0.9, 0.3, 1.8),
                     "halo_lfloor": (9.4, 8.0, 11.0),
                     "dlogc": (0.0, -0.4, 0.3)}
    LOCAL_PARAMS = dict(LOCAL_UPS)

    def __init__(self):
        self._cache, self._lum = {}, {}

    def halo_mass(self, mb, gp):
        return 10 ** gp["halo_A"] * (mb / 1e10) ** gp["halo_B"] + 10 ** gp["halo_lfloor"]

    def dc14(self, mb, mstar, gp):
        mh = self.halo_mass(mb, gp)
        c = nfw_concentration(mh) * 10 ** gp["dlogc"]
        X = np.log10(max(mstar, 1.0) / mh)
        key = (round(np.log10(mh), 4), round(np.log10(c), 4), round(X, 3))
        h = self._cache.get(key)
        if h is None:
            if len(self._cache) > 50000:
                self._cache.clear()
            h = self._cache[key] = DC14Halo(mh, c, X)
        return h

    # rotation curves (01, 02, 04)
    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp.get("ups_scale", 1.0)
        R = np.asarray(R, dtype=float)
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        if getattr(g, "ups_fixed", False):
            mstar = g.stellar_mass()
        else:
            if g.name not in self._lum:
                self._lum[g.name] = _lum_total(g)
            Ld, Lb = self._lum[g.name]
            mstar = u * Ld + 1.4 * u * Lb
        h = self.dc14(mstar + _gas_mass(g), mstar, gp)
        return np.sqrt(np.clip(vb2 + G * h.mass(R) / R, 0, None))

    # dwarfs (05, 10): stars + DC14 halo (cuspy for ultra-faints, cored near X ~ -2.7)
    def predict_dwarf_dispersion(self, d, gp):
        mstar = max(d.M_bar - getattr(d, "M_gas", 0.0), 1.0)
        h = self.dc14(d.M_bar, mstar, gp)
        return dispersion_from_mass(d.M_bar / 2 + float(h.mass(d.r_half)), d.r_half)

    # weak lensing around isolated galaxies (07)
    def predict_lensing_vcirc(self, M_bar, R_kpc, gp):
        R = np.asarray(R_kpc, dtype=float)
        h = self.dc14(M_bar, M_bar, gp)
        return np.sqrt(G * (M_bar + h.mass(R)) / R)

    # strong-lensing ellipticals (08): contracted halos (X near -1.3 or above)
    def predict_lens_mass(self, L, gp):
        mstar = gp["imf_ell"] * L.Mstar_chab
        h = self.dc14(mstar, mstar, gp)
        return gp["imf_ell"] * L.Mstar_proj_chab + h.mass_projected(L.r_Ein)


MODEL = GalForm()
