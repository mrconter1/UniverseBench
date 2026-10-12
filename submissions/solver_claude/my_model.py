"""
Collapsing self-interacting dark matter (cSIDM).

Dark matter is a cold particle that scatters off itself with a cross-section
that falls with speed (light-mediator, Rutherford-like form):

    sigma/m (v) = sigma0 / (1 + (v / w)^2)^2.

One curve gives three regimes:
  * galaxies (V ~ 50-300 km/s): a few cm^2/g. Scattering thermalises the
    inner halo inside the total (dark + baryonic) potential, so the dark
    matter follows where the ordinary matter sits (Kaplinghat+ 2014). This is
    how the galaxy-scale regularity emerges without a new force.
  * dwarfs (V ~ 10-30 km/s): ~100 cm^2/g. Their halos have already gone
    through core formation into gravothermal collapse and are dense again
    (Balberg+ 2002; Essig+ 2019; Nishikawa+ 2020), so dwarf spheroidals look
    like cuspy cold-dark-matter halos.
  * clusters (V ~ 1000 km/s): ~0.1 cm^2/g or less, effectively collisionless
    (Bullet Cluster, cluster masses, cosmology are those of cold dark matter).

Galaxy-scale constants were fitted by the author on the public half of the
benchmark's data (rotation curves, dwarfs and galaxy lensing jointly) and
are frozen here (see THEORY.md); the benchmark fits the rest.
"""
import numpy as np
from scipy.optimize import brentq

from ubench.models.common import (COSMO, COSMO_CDM, CosmologyMixin, ELL, GeneralRelativity,
                                  LOCAL_UPS, v_baryon_sq)
from ubench.physics import G, nfw_concentration, nfw_r200, nfw_mass_projected

try:
    from numba import njit
except ImportError:                                   # pure-python fallback (slow)
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap(a[0]) if a and callable(a[0]) else wrap

# ---------------------------------------------------------------- frozen constants
UPS_DISK = 0.463           # 3.6 um stellar mass-to-light (disc); bulge 1.4x
HALO_A = 11.60             # log10 halo mass of a galaxy with 1e10 Msun of baryons
HALO_B = 0.884             # slope d log M_halo / d log M_baryon
HALO_FLOOR = 10 ** 9.431   # minimum halo mass hosting a galaxy (reionisation / cooling)
DLOGC = -0.123             # halo concentration offset from Dutton & Maccio 2014 (dex)
SIGMA0 = 10 ** 2.106       # cm^2/g, low-speed cross-section
W_TURN = 10 ** 2.178       # km/s, turnover speed
T_AGE = 10.0               # Gyr, typical halo age
TAU_A, TAU_B = 0.2, 1.0    # collapse progress over which the core turns back into a cusp

KMS_KPC_GYR = 1.0227
CM2G_TO_KPC2_MSUN = 2.0886e-10
LRB = np.linspace(np.log(1e-5), np.log(3e3), 400)


# ---------------------------------------------------------------- halo helpers
def nfw_rs_rhos(m200, c, r200):
    rs = r200 / c
    mu = np.log(1 + c) - c / (1 + c)
    return rs, m200 / (4 * np.pi * rs ** 3 * mu)


def nfw_rho(r, m200, c, r200):
    rs, rhos = nfw_rs_rhos(m200, c, r200)
    x = r / rs
    return rhos / (x * (1 + x) ** 2)


def nfw_m(r, m200, c, r200):
    x = r * c / r200
    mu = np.log(1 + c) - c / (1 + c)
    return m200 * (np.log(1 + x) - x / (1 + x)) / mu


def nfw_vmax(m200, c, r200):
    rs, _ = nfw_rs_rhos(m200, c, r200)
    r = 2.163 * rs
    return np.sqrt(G * nfw_m(r, m200, c, r200) / r)


def sigma_m(v):
    return SIGMA0 / (1 + (v / W_TURN) ** 2) ** 2


_VG = np.linspace(1e-3, 12, 400)


def sigma_eff(vmax):
    """Effective cross-section of a halo (Yang & Yu 2022): sigma(v) averaged
    with weight v^5 over a Maxwellian of relative speeds, nu = 0.64 Vmax."""
    nu = 0.64 * vmax
    v = _VG * nu
    wt = v ** 7 * np.exp(-v ** 2 / (4 * nu ** 2))
    return float(np.sum(sigma_m(v) * wt) / np.sum(wt))


def collapse_tau(m200, c, r200, se):
    """t_age / t_c, t_c = (150/C) / (sigma rho_s r_s) / sqrt(4 pi G rho_s), C=0.75."""
    rs, rhos = nfw_rs_rhos(m200, c, r200)
    t_c = (150 / 0.75) / (se * CM2G_TO_KPC2_MSUN * rhos * rs) / (np.sqrt(4 * np.pi * G * rhos) * KMS_KPC_GYR)
    return T_AGE / t_c


def halo_mass(mb):
    return 10 ** HALO_A * (mb / 1e10) ** HALO_B + HALO_FLOOR


def concentration(mh):
    return nfw_concentration(mh) * 10 ** DLOGC


def _cubic_root(K):
    p, q = -1.0 / 3.0, -2.0 / 27.0 - K
    sd = np.sqrt(max((q / 2) ** 2 + (p / 3) ** 3, 0.0))
    return np.cbrt(-q / 2 + sd) + np.cbrt(-q / 2 - sd) - 2.0 / 3.0


@njit(cache=False)
def _shoot(s0s, r1s, l0s, lrb, rhob, n, Gc):
    k = s0s.shape[0]
    hout = np.empty(k)
    mout = np.empty(k)
    for j in range(k):
        s0, r1, h = s0s[j], r1s[j], l0s[j]
        K = 4.0 * np.pi * Gc / (s0 * s0)
        lr0 = np.log(r1) - 9.0
        dl = 9.0 / (n - 1)
        u = 0.0
        m = 4.0 / 3.0 * np.pi * np.exp(lr0) ** 3 * np.exp(h)
        for i in range(n - 1):
            lra = lr0 + i * dl
            lrm = lra + 0.5 * dl
            lrb_ = lra + dl
            ra, rm, rb = np.exp(lra), np.exp(lrm), np.exp(lrb_)
            pa = np.interp(lra, lrb, rhob)
            pm = np.interp(lrm, lrb, rhob)
            pb = np.interp(lrb_, lrb, rhob)
            k1h = u / ra
            k1u = -K * ra ** 3 * (np.exp(h) + pa)
            h2, u2 = h + 0.5 * dl * k1h, u + 0.5 * dl * k1u
            k2h = u2 / rm
            k2u = -K * rm ** 3 * (np.exp(h2) + pm)
            h3, u3 = h + 0.5 * dl * k2h, u + 0.5 * dl * k2u
            k3h = u3 / rm
            k3u = -K * rm ** 3 * (np.exp(h3) + pm)
            h4, u4 = h + dl * k3h, u + dl * k3u
            k4h = u4 / rb
            k4u = -K * rb ** 3 * (np.exp(h4) + pb)
            hn = h + dl / 6.0 * (k1h + 2 * k2h + 2 * k3h + k4h)
            u = u + dl / 6.0 * (k1u + 2 * k2u + 2 * k3u + k4u)
            m += 4.0 * np.pi * 0.5 * (ra ** 3 * np.exp(h) + rb ** 3 * np.exp(hn)) * dl
            h = hn
        hout[j] = h
        mout[j] = m
    return hout, mout


@njit(cache=False)
def _profile(s0, r1, l0, lrb, rhob, n, Gc):
    K = 4.0 * np.pi * Gc / (s0 * s0)
    lr0 = np.log(r1) - 9.0
    dl = 9.0 / (n - 1)
    rr = np.empty(n)
    mm = np.empty(n)
    h, u = l0, 0.0
    r = np.exp(lr0)
    m = 4.0 / 3.0 * np.pi * r ** 3 * np.exp(h)
    rr[0], mm[0] = r, m
    for i in range(n - 1):
        lra = lr0 + i * dl
        lrm = lra + 0.5 * dl
        lrb_ = lra + dl
        ra, rm, rb = np.exp(lra), np.exp(lrm), np.exp(lrb_)
        pa = np.interp(lra, lrb, rhob)
        pm = np.interp(lrm, lrb, rhob)
        pb = np.interp(lrb_, lrb, rhob)
        k1h = u / ra
        k1u = -K * ra ** 3 * (np.exp(h) + pa)
        h2, u2 = h + 0.5 * dl * k1h, u + 0.5 * dl * k1u
        k2h = u2 / rm
        k2u = -K * rm ** 3 * (np.exp(h2) + pm)
        h3, u3 = h + 0.5 * dl * k2h, u + 0.5 * dl * k2u
        k3h = u3 / rm
        k3u = -K * rm ** 3 * (np.exp(h3) + pm)
        h4, u4 = h + dl * k3h, u + dl * k3u
        k4h = u4 / rb
        k4u = -K * rb ** 3 * (np.exp(h4) + pb)
        hn = h + dl / 6.0 * (k1h + 2 * k2h + 2 * k3h + k4h)
        u = u + dl / 6.0 * (k1u + 2 * k2u + 2 * k3u + k4u)
        m += 4.0 * np.pi * 0.5 * (ra ** 3 * np.exp(h) + rb ** 3 * np.exp(hn)) * dl
        h = hn
        rr[i + 1], mm[i + 1] = rb, m
    return rr, mm


class SIDMHalo:
    """NFW outside r1; isothermal (dispersion sigma0) inside r1 in the total
    potential, matched in density and mass at r1, where each particle has
    scattered once in T_AGE. If no cored solution exists (scattering region
    far beyond r_s), the halo is past its core-forming phase and is treated
    as collapsing (NFW-like). The core is blended back toward the cusp as the
    gravothermal clock tau advances."""
    N = 160

    def __init__(self, m200, c, mb_of_r, guess=None):
        self.m200, self.c, self.r200 = m200, c, nfw_r200(m200)
        r = np.exp(LRB)
        mb = np.maximum.accumulate(np.clip(mb_of_r(r), 0, None))
        self.rhob = np.clip(np.gradient(mb, r) / (4 * np.pi * r ** 2), 0, None)
        self.se = sigma_eff(nfw_vmax(m200, c, self.r200))
        self.tau = collapse_tau(m200, c, self.r200, self.se)
        self.ok = False
        if guess is None or not self._polish(guess):
            self._solve()

    def r1_of(self, s0):
        v = 4 / np.sqrt(np.pi) * s0
        k = self.se * CM2G_TO_KPC2_MSUN * v * KMS_KPC_GYR * T_AGE
        rs, rhos = nfw_rs_rhos(self.m200, self.c, self.r200)
        r1 = _cubic_root(rhos * k) * rs
        if r1 < 1e-4 * self.r200:
            return None
        return min(r1, 3 * self.r200)

    def _F(self, x):
        s0, l0 = x
        if s0 <= 0:
            return None
        r1 = self.r1_of(s0)
        if r1 is None:
            return None
        hh, mm = _shoot(np.array([s0]), np.array([r1]), np.array([l0]), LRB, self.rhob, self.N, G)
        return np.array([hh[0] - np.log(nfw_rho(r1, self.m200, self.c, self.r200)),
                         np.log(mm[0] / nfw_m(r1, self.m200, self.c, self.r200))])

    def _finish(self, s0, l0):
        r1 = self.r1_of(s0)
        rr, mm = _profile(s0, r1, l0, LRB, self.rhob, self.N, G)
        self.sigma0, self.r1, self.rho0 = s0, r1, np.exp(l0)
        self.lr_in, self.lm_in = np.log(rr), np.log(mm)
        self.ok = True

    def _polish(self, guess):
        x = np.array(guess, dtype=float)
        for _ in range(12):
            F = self._F(x)
            if F is None or not np.all(np.isfinite(F)):
                return False
            if np.max(np.abs(F)) < 1e-4:
                break
            J = np.empty((2, 2))
            for j, hstep in enumerate((1e-4 * x[0], 1e-4)):
                xx = x.copy()
                xx[j] += hstep
                Fj = self._F(xx)
                if Fj is None:
                    return False
                J[:, j] = (Fj - F) / hstep
            try:
                dx = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                return False
            sc = min(1.0, 0.2 * x[0] / max(abs(dx[0]), 1e-12), 0.5 / max(abs(dx[1]), 1e-12))
            x = x + sc * dx
        else:
            return False
        F = self._F(x)
        if F is None or np.max(np.abs(F)) > 1e-3:
            return False
        self._finish(*x)
        return True

    def _scan(self, ns0, noff, n):
        m200, c, r200 = self.m200, self.c, self.r200
        vmax = 1.2 * np.sqrt(G * m200 / r200)
        offs = np.linspace(0, 9, noff)
        pts = []
        for s0 in np.linspace(0.15, 1.6, ns0) * vmax:
            r1 = self.r1_of(s0)
            if r1 is None:
                continue
            rho1, m1 = nfw_rho(r1, m200, c, r200), nfw_m(r1, m200, c, r200)
            L = np.log(rho1) + offs
            h, m = _shoot(np.full(noff, s0), np.full(noff, r1), L, LRB, self.rhob, n, G)
            d = h - np.log(rho1)
            mm = np.log(m / m1)
            for j in np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]:
                t = d[j] / (d[j] - d[j + 1])
                pts.append((L[j] + t * (L[j + 1] - L[j]), s0, mm[j] + t * (mm[j + 1] - mm[j])))
        return pts

    def _solve(self):
        pts = self._scan(16, 50, self.N)
        if len(pts) < 2:
            return
        pts = np.array(sorted(pts))
        cross = np.where(np.sign(pts[:-1, 2]) != np.sign(pts[1:, 2]))[0]
        if len(cross) == 0:
            return
        i = cross[0]
        t = pts[i, 2] / (pts[i, 2] - pts[i + 1, 2])
        x0 = [pts[i, 1] + t * (pts[i + 1, 1] - pts[i, 1]), pts[i, 0] + t * (pts[i + 1, 0] - pts[i, 0])]
        if not self._polish(x0):
            self._finish(*x0)

    def mass(self, r):
        r = np.asarray(r, dtype=float)
        mn = nfw_m(r, self.m200, self.c, self.r200)
        if not self.ok:
            return mn
        mc = np.where(r < self.r1, np.exp(np.interp(np.log(r), self.lr_in, self.lm_in)), mn)
        s = float(np.clip((self.tau - TAU_A) / (TAU_B - TAU_A), 0.0, 1.0))
        return np.exp((1 - s) * np.log(mc) + s * np.log(mn))


# ---------------------------------------------------------------- baryon bookkeeping
def _gas_mass(g):
    return max(float(g.Vgas[-1] * abs(g.Vgas[-1]) * g.R[-1] / G), 0.0)


def _lum_total(g):
    """Disc and bulge luminosity with an exponential tail beyond the last radius."""
    R, sb = g.R, g.SBdisk
    Ld, Lb = g.luminosity()
    k = (sb > 0) & (R >= R[len(R) // 2])
    if k.sum() >= 3:
        fit = np.polyfit(R[k], np.log(sb[k]), 1)
        if fit[0] < -1e-3:
            h, Rl = -1 / fit[0] * 1000.0, R[-1] * 1000.0
            Ld += 2 * np.pi * np.exp(np.polyval(fit, R[-1])) * h * (Rl + h)
    return Ld, Lb


def _mb_profile(g, u):
    vb2 = np.clip(v_baryon_sq(g, u, 1.4 * u), 0, None)
    M = np.maximum(np.maximum.accumulate(g.R * vb2 / G), 1e-3)
    R = g.R

    def f(r):
        r = np.asarray(r, dtype=float)
        out = np.exp(np.interp(np.log(r), np.log(R), np.log(M)))
        return np.where(r < R[0], M[0] * (r / R[0]) ** 3, out)
    return f


def _plummer(M, rh3):
    a = rh3 / 1.305
    return lambda r: M * r ** 3 / (r ** 2 + a ** 2) ** 1.5


class CollapsingSIDM(GeneralRelativity, CosmologyMixin):
    NAME = "Collapsing self-interacting dark matter (cSIDM)"
    GLOBAL_PARAMS = {**ELL, **COSMO, **COSMO_CDM,
                     "f_b_cluster": (0.13, 0.03, 0.3),
                     "fb_mass_slope": (0.1, -0.5, 0.5),
                     "fb_radius_slope": (0.1, -0.5, 0.5)}
    LOCAL_PARAMS = {**LOCAL_UPS}

    def __init__(self):
        self._warm, self._lum = {}, {}

    def _halo(self, mh, mb_of_r, key):
        h = SIDMHalo(mh, concentration(mh), mb_of_r, guess=self._warm.get(key))
        if h.ok:
            self._warm[key] = (h.sigma0, np.log(h.rho0))
        return h

    # galaxies (01, 02, 04)
    def predict(self, g, R, gp, lp):
        u = UPS_DISK * lp.get("ups_scale", 1.0)
        R = np.asarray(R, dtype=float)
        vb2 = v_baryon_sq(g, u, 1.4 * u)
        if getattr(g, "ups_fixed", False):
            mb = g.stellar_mass() + _gas_mass(g)
        else:
            if g.name not in self._lum:
                self._lum[g.name] = _lum_total(g)
            Ld, Lb = self._lum[g.name]
            mb = u * Ld + 1.4 * u * Lb + _gas_mass(g)
        key = (g.name, len(g.R), round(u, 4))
        h = self._halo(halo_mass(mb), _mb_profile(g, 1.0 if getattr(g, "ups_fixed", False) else u),
                       (g.name, len(g.R)))
        return np.sqrt(np.clip(vb2 + G * h.mass(R) / R, 0, None))

    # dwarfs (05, 10)
    def predict_dwarf_dispersion(self, d, gp):
        h = self._halo(halo_mass(d.M_bar), _plummer(d.M_bar, d.r_half), ("dwarf", d.name))
        return float(np.sqrt(G * (d.M_bar / 2 + float(h.mass(np.array([d.r_half]))[0])) / (3 * d.r_half)))

    # weak lensing around isolated galaxies (07): outer halo, unaffected by scattering
    def predict_lensing_vcirc(self, M_bar, R_kpc, gp):
        R = np.asarray(R_kpc, dtype=float)
        mh = halo_mass(M_bar)
        c = concentration(mh)
        return np.sqrt(G * (M_bar + nfw_m(R, mh, c, nfw_r200(mh))) / R)

    # strong-lensing ellipticals (08): ~300 km/s halos, sigma/m < 1: effectively NFW
    def predict_lens_mass(self, L, gp):
        mstar = gp["imf_ell"] * L.Mstar_chab
        mh = halo_mass(mstar)
        return gp["imf_ell"] * L.Mstar_proj_chab + nfw_mass_projected(L.r_Ein, mh, concentration(mh))

    # clusters (03, 09): sigma/m ~ 0.05 cm^2/g at 1000 km/s: cold dark matter
    def cluster_baryon_fraction(self, M_bar, gp, delta=500):
        return (gp["f_b_cluster"] * (M_bar / 1e13) ** gp["fb_mass_slope"]
                * (delta / 500) ** (-gp["fb_radius_slope"]))

    def predict_cluster_mass(self, c, gp):
        return c.M_bar / self.cluster_baryon_fraction(c.M_bar, gp, c.delta)

    def predict_bullet_offset(self, s, gp):
        M_bar = s.M_gas + s.M_star
        dm = M_bar * (1 / self.cluster_baryon_fraction(M_bar, gp) - 1)
        return s.M_gas / (s.M_gas + s.M_star + dm)

    def predict_wide_binary_boost(self, gN_si, g_ext_si, gp):
        return 1.0

    def dm_particle(self, gp):
        """Cold, self-interacting through a light dark mediator; no coupling to
        Standard Model particles is required (none assumed)."""
        return {"type": "cold, self-interacting", "wdm_mass_keV": float("inf"),
                "sigma_SI_cm2": 0.0, "sigma_v_cm3s": 0.0, "mass_GeV": None}


MODEL = CollapsingSIDM()
