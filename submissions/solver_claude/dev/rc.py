"""Fast rotation-curve harness on the public SPARC galaxies (test 01 logic)."""
import sys
import numpy as np
sys.path[:0] = ['.', 'dev']
from ubench.tests.rotation import load_sparc
from ubench.models.common import v_baryon_sq
from ubench.physics import G, rar_nu, nfw_concentration
from halos import Profile, dc14_shape, dc14_conc_boost, halo_from_baryons

GALS = load_sparc()
A0 = 1.2e-10 * 3.0857e13


def gas_mass(g):
    return max(float(g.Vgas[-1] * abs(g.Vgas[-1]) * g.R[-1] / G), 0.0)


def zs(vfun, gals=GALS):
    out = []
    for g in gals:
        v = vfun(g)
        out.append((np.asarray(v) - g.Vobs) / g.sigma)
    return out


def report(z, label=""):
    allz = np.concatenate(z)
    pct = 100 * np.mean(np.abs(allz) < 2)
    print(f"{label:40s} within2s {pct:5.1f}%  mean z2 {np.mean(allz**2):7.2f}  median|z| {np.median(np.abs(allz)):.2f}")
    return pct


def v_mond(g, ups=0.5, a0=A0):
    vb2 = v_baryon_sq(g, ups, 1.4 * ups)
    gN = np.clip(vb2, 1e-6, None) / g.R
    return np.sqrt(gN * rar_nu(gN / a0) * g.R)


def v_halo_model(g, p, profile="dc14"):
    """p: dict with ups, lN, lM1, beta, gam, dlogc, (gas_in_mb)"""
    u = p["ups"]
    vb2 = v_baryon_sq(g, u, 1.4 * u)
    ms = g.stellar_mass(u)
    mb = ms + gas_mass(g)
    mh = halo_from_baryons(mb, (p["lN"], p["lM1"], p["beta"], p["gam"]))
    c = nfw_concentration(mh) * 10 ** p["dlogc"]
    if profile == "dc14":
        X = np.log10(ms / mh)
        a, b, gm = dc14_shape(X)
        prof = Profile(mh, c * dc14_conc_boost(X), a, b, gm)
    else:
        prof = Profile(mh, c)
    return np.sqrt(np.clip(vb2, 0, None) + prof.vcirc(g.R) ** 2)


from halos import CoreNFW, half_light_radius
_RH = {g.name: half_light_radius(g) for g in GALS}


def v_corenfw(g, p):
    u = p["ups"]
    vb2 = v_baryon_sq(g, u, 1.4 * u)
    ms = g.stellar_mass(u)
    mb = ms + gas_mass(g)
    mh = halo_from_baryons(mb, (p["lN"], p["lM1"], p["beta"], p["gam"]))
    c = nfw_concentration(mh) * 10 ** p["dlogc"]
    prof = CoreNFW(mh, c, 10 ** p["leta"] * _RH[g.name], 1.0)
    return np.sqrt(np.clip(vb2, 0, None) + prof.vcirc(g.R) ** 2)


from sidm import SIDMHalo, collapse_tau
from ubench.physics import nfw_r200


def mb_profile(g, u):
    """Spherical-equivalent baryonic mass M_b(<r) = r v_b^2 / G from the data,
    r^3 inside the first point, constant beyond the last."""
    vb2 = np.clip(v_baryon_sq(g, u, 1.4 * u), 0, None)
    M = np.maximum.accumulate(g.R * vb2 / G)
    R = g.R
    M = np.maximum(M, 1e-3)

    def f(r):
        r = np.asarray(r, dtype=float)
        out = np.exp(np.interp(np.log(r), np.log(R), np.log(M)))
        out = np.where(r < R[0], M[0] * (r / R[0]) ** 3, out)
        return out
    return f


_WARM = {}


def v_sidm(g, p, t_age=10.0):
    u = p["ups"]
    vb2 = v_baryon_sq(g, u, 1.4 * u)
    ms = g.stellar_mass(u)
    mb = ms + gas_mass(g)
    mh = halo_from_baryons(mb, (p["lN"], p["lM1"], p["beta"], p["gam"]))
    c = nfw_concentration(mh) * 10 ** p["dlogc"]
    w = 10 ** p["lw"] if "lw" in p else 1e5
    h = SIDMHalo(mh, c, nfw_r200(mh), mb_profile(g, u), p["lsig"], w, t_age,
                 guess=_WARM.get(g.name))
    if h.ok:
        _WARM[g.name] = (h.sigma0, np.log(h.rho0))
    return np.sqrt(np.clip(vb2, 0, None) + dm_mass_blend(h, g.R, p["lsig"], w, t_age) * G / g.R)


def collapse_weight(tau, ta=0.2, tb=1.0):
    return float(np.clip((tau - ta) / (tb - ta), 0.0, 1.0))


def dm_mass_blend(h, r, lsig, w, t_age=10.0):
    """Core phase (isothermal Jeans with baryons) blended in log mass toward the
    collapsed state (approximated by NFW) as gravothermal evolution proceeds."""
    from sidm import collapse_tau, nfw_m
    tau = collapse_tau(h.m200, h.c, h.r200, lsig, w, t_age)
    s = collapse_weight(tau)
    mc = h.mass(r)
    mn = nfw_m(np.asarray(r, dtype=float), h.m200, h.c, h.r200)
    return np.exp((1 - s) * np.log(mc) + s * np.log(mn))
