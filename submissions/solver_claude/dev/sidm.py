"""Semi-analytic SIDM halos (Kaplinghat, Keeley, Linden & Yu 2014, PRL 113,
021302; Kamada+ 2017; Ren+ 2019).

Outside r1 the halo is the cold-dark-matter NFW halo. Inside r1, where each
particle has scattered about once over the halo's age t,
    rho_NFW(r1) * (sigma/m)(v) * <v> * t = 1,   <v> = (4/sqrt(pi)) sigma0,
the dark matter is isothermal with 1D dispersion sigma0 in the TOTAL
potential (dark matter + baryons):
    (1/r^2) d/dr (r^2 d ln rho / dr) = -4 pi G (rho + rho_b) / sigma0^2.
(rho0, sigma0) are fixed by continuity of density and enclosed mass at r1.
Of the solutions, the least dense (cored) one is taken.

Units: kpc, km/s, Msun, Gyr.
"""
import numpy as np
from scipy.optimize import brentq, root
from ubench.physics import G

try:
    from numba import njit
except ImportError:                       # pure-python fallback (slow)
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap(a[0]) if a and callable(a[0]) else wrap

KMS_KPC_GYR = 1.0227
CM2G_TO_KPC2_MSUN = 2.0886e-10
LRB = np.linspace(np.log(1e-5), np.log(3e3), 400)     # log-r grid for baryon density


def nfw_rho(r, m200, c, r200):
    rs = r200 / c
    mu = np.log(1 + c) - c / (1 + c)
    x = r / rs
    return m200 / (4 * np.pi * rs ** 3 * mu) / (x * (1 + x) ** 2)


def nfw_m(r, m200, c, r200):
    x = r * c / r200
    mu = np.log(1 + c) - c / (1 + c)
    return m200 * (np.log(1 + x) - x / (1 + x)) / mu


def _cubic_root(K):
    """Positive real root of x^3 + 2x^2 + x - K = 0 (Cardano; one real root for K > 0)."""
    # substitute x = y - 2/3: y^3 + p y + q = 0
    p = -1.0 / 3.0
    q = -2.0 / 27.0 - K
    disc = (q / 2) ** 2 + (p / 3) ** 3
    sd = np.sqrt(disc) if disc > 0 else 0.0
    y = np.cbrt(-q / 2 + sd) + np.cbrt(-q / 2 - sd)
    return y - 2.0 / 3.0


def sigma_m(v, lsig, w):
    """sigma/m [cm^2/g] at relative speed v [km/s]: 10^lsig at low speed,
    falling as v^-4 above the turnover speed w."""
    return 10 ** lsig / (1 + (v / w) ** 2) ** 2


def baryon_density_table(mb_of_r):
    r = np.exp(LRB)
    m = np.maximum.accumulate(np.clip(mb_of_r(r), 0, None))
    dm = np.gradient(m, r)
    return np.clip(dm / (4 * np.pi * r ** 2), 0, None)


@njit(cache=True)
def _shoot(s0s, r1s, l0s, lrb, rhob, n, Gc):
    """Integrate many isothermal trajectories. Returns ln rho and mass at r1
    for each (s0, r1, l0) triple."""
    k = s0s.shape[0]
    hout = np.empty(k)
    mout = np.empty(k)
    for j in range(k):
        s0, r1, h = s0s[j], r1s[j], l0s[j]
        K = 4.0 * np.pi * Gc / (s0 * s0)
        lr0 = np.log(r1) - 9.0
        dl = 9.0 / (n - 1)
        u = 0.0
        r = np.exp(lr0)
        m = 4.0 / 3.0 * np.pi * r ** 3 * np.exp(h)
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


@njit(cache=True)
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
    N = 160

    def __init__(self, m200, c, r200, mb_of_r, lsig, w, t_age=10.0, rhob=None, guess=None):
        self.m200, self.c, self.r200 = m200, c, r200
        self.rhob = baryon_density_table(mb_of_r) if rhob is None else rhob
        self.ok = False
        if guess is None or not self._polish(guess, lsig, w, t_age):
            self._solve(lsig, w, t_age)

    def _F(self, x, lsig, w, t_age):
        s0, l0 = x
        if s0 <= 0:
            return None
        r1 = self.r1_of(s0, lsig, w, t_age)
        if r1 is None:
            return None
        hh, mm = _shoot(np.array([s0]), np.array([r1]), np.array([l0]), LRB, self.rhob, self.N, G)
        return np.array([hh[0] - np.log(nfw_rho(r1, self.m200, self.c, self.r200)),
                         np.log(mm[0] / nfw_m(r1, self.m200, self.c, self.r200))])

    def _polish(self, guess, lsig, w, t_age):
        """Newton from a previous solution (warm start); False if it fails."""
        x = np.array(guess, dtype=float)
        for _ in range(12):
            F = self._F(x, lsig, w, t_age)
            if F is None or not np.all(np.isfinite(F)):
                return False
            if np.max(np.abs(F)) < 1e-4:
                break
            J = np.empty((2, 2))
            for j, hstep in enumerate((1e-4 * x[0], 1e-4)):
                xx = x.copy(); xx[j] += hstep
                Fj = self._F(xx, lsig, w, t_age)
                if Fj is None:
                    return False
                J[:, j] = (Fj - F) / hstep
            try:
                dx = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                return False
            # damp large steps (stay on the same branch)
            sc = min(1.0, 0.2 * x[0] / max(abs(dx[0]), 1e-12), 0.5 / max(abs(dx[1]), 1e-12))
            x = x + sc * dx
        else:
            return False
        if np.max(np.abs(self._F(x, lsig, w, t_age))) > 1e-3:
            return False
        s0, l0 = x
        r1 = self.r1_of(s0, lsig, w, t_age)
        rr, mm = _profile(s0, r1, l0, LRB, self.rhob, self.N, G)
        self.sigma0, self.r1, self.rho0 = s0, r1, np.exp(l0)
        self.lr_in, self.lm_in = np.log(rr), np.log(mm)
        self.ok = True
        return True

    def r1_of(self, s0, lsig, w, t_age):
        """r1 where rho_NFW(r1) * sigma/m * <v> * t = 1: solves x(1+x)^2 = K,
        x = r1/r_s, in closed form (one real root of the cubic)."""
        v = 4 / np.sqrt(np.pi) * s0
        k = sigma_m(v, lsig, w) * CM2G_TO_KPC2_MSUN * v * KMS_KPC_GYR * t_age
        rs, rhos = nfw_rs_rhos(self.m200, self.c, self.r200)
        K = rhos * k                         # = rho_s / target density
        x = _cubic_root(K)
        r1 = x * rs
        if r1 < 1e-4 * self.r200:
            return None
        return min(r1, 3 * self.r200)

    def _solve(self, lsig, w, t_age, ns0=16, noff=50, coarse=True):
        """Scan for the cored solution. If a coarse scan finds no matching
        branch, there is none: the scattering region would extend well beyond
        the scale radius, i.e. the halo is past its core-forming phase
        (treated as collapsing, NFW-like, by the caller)."""
        if coarse and not self._has_branch(lsig, w, t_age):
            return
        m200, c, r200 = self.m200, self.c, self.r200
        vmax = 1.2 * np.sqrt(G * m200 / r200)
        s0g = np.linspace(0.15, 1.6, ns0) * vmax
        offs = np.linspace(0, 9, noff)
        S, R1, L, M1, RHO1, ISO = [], [], [], [], [], []
        for s0 in s0g:
            r1 = self.r1_of(s0, lsig, w, t_age)
            if r1 is None:
                continue
            rho1, m1 = nfw_rho(r1, m200, c, r200), nfw_m(r1, m200, c, r200)
            S.append(np.full(len(offs), s0)); R1.append(np.full(len(offs), r1))
            L.append(np.log(rho1) + offs); M1.append(m1); RHO1.append(rho1)
        if not S:
            return                                   # no scattering: plain NFW
        S, R1, L = np.concatenate(S), np.concatenate(R1), np.concatenate(L)
        h, m = _shoot(S, R1, L, LRB, self.rhob, self.N, G)
        k = len(offs)
        pts = []
        for i in range(len(M1)):
            d = h[i * k:(i + 1) * k] - np.log(RHO1[i])
            mm = np.log(m[i * k:(i + 1) * k] / M1[i])
            idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
            for j in idx:
                t = d[i * k + 0 - i * k + j] / (d[j] - d[j + 1]) if False else d[j] / (d[j] - d[j + 1])
                l0 = L[i * k + j] + t * (L[i * k + j + 1] - L[i * k + j])
                pts.append((l0, S[i * k], mm[j] + t * (mm[j + 1] - mm[j])))
        if len(pts) < 2:
            return
        pts.sort()
        pts = np.array(pts)
        cross = np.where(np.sign(pts[:-1, 2]) != np.sign(pts[1:, 2]))[0]
        if len(cross) == 0:
            return
        i = cross[0]
        t = pts[i, 2] / (pts[i, 2] - pts[i + 1, 2])
        x0 = [pts[i, 1] + t * (pts[i + 1, 1] - pts[i, 1]), pts[i, 0] + t * (pts[i + 1, 0] - pts[i, 0])]

        def F(x):
            s0, l0 = x
            r1 = self.r1_of(s0, lsig, w, t_age)
            if r1 is None or s0 <= 0:
                return [10.0, 10.0]
            hh, mm = _shoot(np.array([s0]), np.array([r1]), np.array([l0]), LRB, self.rhob, self.N, G)
            return [hh[0] - np.log(nfw_rho(r1, m200, c, r200)), np.log(mm[0] / nfw_m(r1, m200, c, r200))]
        sol = root(F, x0, method="hybr", options={"xtol": 1e-6})
        if not sol.success or max(abs(np.array(F(sol.x)))) > 1e-3:
            sol_x = x0                               # fall back to the scan estimate
        else:
            sol_x = sol.x
        s0, l0 = sol_x
        r1 = self.r1_of(s0, lsig, w, t_age)
        rr, mm = _profile(s0, r1, l0, LRB, self.rhob, self.N, G)
        self.sigma0, self.r1, self.rho0 = s0, r1, np.exp(l0)
        self.lr_in, self.lm_in = np.log(rr), np.log(mm)
        self.ok = True

    def _has_branch(self, lsig, w, t_age):
        """Cheap test: does the mass mismatch change sign anywhere along the
        density-matching curve? (8 x 24 trajectories, coarse steps)"""
        m200, c, r200 = self.m200, self.c, self.r200
        vmax = 1.2 * np.sqrt(G * m200 / r200)
        offs = np.linspace(0, 9, 24)
        mms = []
        for s0 in np.linspace(0.15, 1.6, 8) * vmax:
            r1 = self.r1_of(s0, lsig, w, t_age)
            if r1 is None:
                continue
            rho1, m1 = nfw_rho(r1, m200, c, r200), nfw_m(r1, m200, c, r200)
            L = np.log(rho1) + offs
            h, m = _shoot(np.full(len(offs), s0), np.full(len(offs), r1), L, LRB, self.rhob, 100, G)
            d = h - np.log(rho1)
            idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
            mms += [np.log(m[j] / m1) for j in idx]
        return len(mms) > 0 and min(mms) < 0.15 and max(mms) > -0.15

    def mass(self, r):
        r = np.asarray(r, dtype=float)
        out = nfw_m(r, self.m200, self.c, self.r200)
        if self.ok:
            inner = np.exp(np.interp(np.log(r), self.lr_in, self.lm_in))
            out = np.where(r < self.r1, inner, out)
        return out

    def vcirc(self, r):
        r = np.asarray(r, dtype=float)
        return np.sqrt(G * self.mass(r) / r)


def nfw_rs_rhos(m200, c, r200):
    rs = r200 / c
    mu = np.log(1 + c) - c / (1 + c)
    return rs, m200 / (4 * np.pi * rs ** 3 * mu)


def nfw_vmax(m200, c, r200):
    rs, _ = nfw_rs_rhos(m200, c, r200)
    r = 2.163 * rs
    return np.sqrt(G * nfw_m(r, m200, c, r200) / r)


_VG = np.linspace(1e-3, 12, 400)          # v / nu grid


def sigma_eff(lsig, w, vmax):
    """Effective cross-section for gravothermal evolution (Yang & Yu 2022):
    sigma(v) averaged with weight v^5 over a Maxwellian of relative speeds,
    1D dispersion nu = 0.64 Vmax."""
    nu = 0.64 * vmax
    v = _VG * nu
    wt = v ** 7 * np.exp(-v ** 2 / (4 * nu ** 2))
    return float(np.sum(sigma_m(v, lsig, w) * wt) / np.sum(wt))


def collapse_tau(m200, c, r200, lsig, w, t_age=10.0, C=0.75):
    """t_age / t_c with t_c = (150/C) / (sigma_eff rho_s r_s) / sqrt(4 pi G rho_s)
    (Balberg, Shapiro & Inagaki 2002; Essig+ 2019)."""
    rs, rhos = nfw_rs_rhos(m200, c, r200)
    se = sigma_eff(lsig, w, nfw_vmax(m200, c, r200)) * CM2G_TO_KPC2_MSUN
    t_c = (150 / C) / (se * rhos * rs) / (np.sqrt(4 * np.pi * G * rhos) * KMS_KPC_GYR)
    return t_age / t_c
