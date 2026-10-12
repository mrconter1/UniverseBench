"""
Two-phase (superfluid) dark matter: one dark-matter field whose phase is set by how hot
and how fast-moving it is.

* Cold, slow, isolated systems (galaxies): the dark matter condenses into a superfluid.
  Its phonons couple to baryons and mediate an extra, MOND-like force
  (g = nu(g_N / a0) g_N), so galaxies follow the radial acceleration relation with no
  per-galaxy halo freedom. Light feels the same potential (disformal coupling).
* Hot systems (clusters, where the dark matter's virial speed exceeds the critical
  velocity V_c) and systems moving supersonically through a host's superfluid
  (satellite dwarfs, colliding clusters; Landau criterion) are in the normal phase:
  ordinary collisionless cold dark matter with a cosmological halo.
* Cosmology: on linear scales both phases are pressureless cold dark matter, so the
  expansion history and structure growth are those of LCDM.

See THEORY.md for the equations and constants.
"""
import numpy as np
from scipy.optimize import brentq

from ubench.models.common import CosmologyMixin, GeneralRelativity, v_baryon_sq
from ubench.physics import (G, dispersion_from_mass, halo_mass_from_stars, nfw_concentration,
                            nfw_mass_enclosed, nfw_mass_projected)

A0K = 3.0857e13          # (m/s^2) -> (km/s)^2/kpc
UB_RATIO = 1.4           # bulge / disc mass-to-light ratio (standard, as baselines)


def nu(y):
    """Phonon-force interpolation, McGaugh+2016 form: g/g_N = 1/(1 - exp(-sqrt(y)))."""
    y = np.maximum(np.asarray(y, dtype=float), 1e-30)
    return 1.0 / (-np.expm1(-np.sqrt(y)))


# ----------------------------------------------------------------- lensing helpers
def _hern_proj_frac(X):
    X = np.atleast_1d(np.asarray(X, dtype=float))
    F = np.ones_like(X)
    lo, hi = X < 1 - 1e-9, X > 1 + 1e-9
    F[lo] = np.arccosh(1 / X[lo]) / np.sqrt(1 - X[lo] ** 2)
    F[hi] = np.arccos(1 / X[hi]) / np.sqrt(X[hi] ** 2 - 1)
    out = X ** 2 * (1 - F) / (X ** 2 - 1)
    out[~(lo | hi)] = 1 / 3
    return out


def _hern_a(frac, R):
    frac = min(max(frac, 1e-4), 0.995)
    f = lambda la: _hern_proj_frac(R / np.exp(la))[0] - frac
    return float(np.exp(brentq(f, np.log(R) - 15, np.log(R) + 15)))


def _projected(Menc, R, fac=3000.0, n=3000):
    """Mass inside a cylinder of radius R for a spherical enclosed-mass profile Menc(r),
    which may grow linearly at large r (isothermal-like phonon halo)."""
    r = np.logspace(np.log10(R), np.log10(R * fac), n)
    M = Menc(r)
    dM = np.gradient(M, r)
    w = 1 - np.sqrt(np.clip(1 - (R / r) ** 2, 0, None))
    return float(M[0] + np.trapezoid(dM * w, r) + dM[-1] * R ** 2 / (2 * r[-1]))


class TwoPhaseDM(GeneralRelativity, CosmologyMixin):
    NAME = "Two-phase superfluid dark matter"
    GLOBAL_PARAMS = {
        "ups_disk": (0.5, 0.2, 1.0),        # 3.6um stellar mass-to-light (disc)
        "imf_ell": (1.6, 0.8, 2.0),         # elliptical stellar mass / Chabrier
        "a0_e10": (1.1, 0.6, 2.0),          # phonon acceleration scale, 1e-10 m/s^2
        "logVc": (2.7, 2.2, 3.3),           # critical (phase-transition) velocity, log10 km/s
        "fb_depl": (0.8, 0.3, 1.2),         # baryon retention in normal-phase halos (x cosmic)
        "R_break": (120.0, 40.0, 400.0),    # condensate-breaking radius around a V=200 km/s host, kpc
        "fb_rslope": (0.1, -0.3, 0.5),      # hot gas more extended than dark matter
        # cosmology (both phases are pressureless CDM on linear scales)
        "H0": (67.4, 60.0, 80.0), "ombh2": (0.0224, 0.018, 0.028), "omch2": (0.120, 0.08, 0.16),
        "ns": (0.965, 0.90, 1.02), "ln10As": (3.044, 2.7, 3.4),
        # dark energy: thawing scalar field, CPL w(a) = w0 + wa (1 - a)
        "w0": (-0.9, -1.5, -0.4), "wa": (-0.5, -2.5, 1.0),
    }
    LOCAL_PARAMS = {"ups_scale": (1.0, 0.33, 3.0)}
    P_SHARP = 4.0                           # sharpness of the phase transition

    # ------------------------------------------------------------- phase physics
    def a0(self, gp):
        return gp["a0_e10"] * 1e-10 * A0K

    def v_phonon(self, M_b, gp):
        """Asymptotic (virial) speed of the superfluid phase: V^4 = G M_b a0."""
        return (G * max(M_b, 1.0) * self.a0(gp)) ** 0.25

    def f_super(self, V_internal, gp, V_bulk=0.0):
        """Superfluid fraction. Two ways to destroy superfluidity: (i) the dark matter is
        too hot (internal speed above V_c); (ii) Landau criterion: the system moves
        through the ambient superfluid faster than its own sound speed (~ V_internal)."""
        Vc = 10 ** gp["logVc"]
        heat = (V_internal / Vc) ** self.P_SHARP
        mach = (V_bulk / max(V_internal, 1e-3)) ** self.P_SHARP
        return 1.0 / (1.0 + heat + mach)

    def f_baryon(self, gp, delta=500.0):
        """Baryon fraction enclosed by a normal-phase halo at overdensity delta: cosmic
        value times retention, falling inward because feedback-heated gas is more
        extended than the collisionless dark matter."""
        return (gp["fb_depl"] * gp["ombh2"] / (gp["ombh2"] + gp["omch2"])
                * (delta / 500.0) ** (-gp["fb_rslope"]))

    # ------------------------------------------------------------- galaxies
    def _galaxy_mb(self, g, u):
        Ms = g.stellar_mass(u)
        Mg = max(float(g.Vgas[-1] * abs(g.Vgas[-1]) * g.R[-1] / G), 0.0)
        return Ms + Mg

    def predict(self, g, R, gp, lp):
        u = gp["ups_disk"] * lp["ups_scale"]
        vb2 = v_baryon_sq(g, u, UB_RATIO * u)
        gN = np.clip(vb2, 1e-6, None) / R
        Mb = self._galaxy_mb(g, u)
        fs = self.f_super(self.v_phonon(Mb, gp), gp)
        v2 = gN * R * (1 + fs * (nu(gN / self.a0(gp)) - 1))
        if fs < 0.999:
            m200 = halo_mass_from_stars(g.stellar_mass(u))
            v2 = v2 + (1 - fs) * G * nfw_mass_enclosed(R, m200, nfw_concentration(m200)) / R
        return np.sqrt(np.clip(v2, 0, None))

    # ------------------------------------------------------------- clusters
    def _hot_mass(self, M_b, gN, gp, delta=500.0):
        """Total dynamical mass for baryonic mass M_b (enclosed) in a system that may be
        partly superfluid: phonon boost on the superfluid share, cosmic DM on the
        normal share."""
        fs = self.f_super(self.v_phonon(M_b, gp), gp)
        sf = M_b * nu(gN / self.a0(gp))
        nm = M_b / self.f_baryon(gp, delta)
        return fs * sf + (1 - fs) * nm, fs

    def predict_cluster_mass(self, c, gp):
        gN = G * c.M_bar / c.r_kpc ** 2
        return float(self._hot_mass(c.M_bar, gN, gp, c.delta)[0])

    def predict_bullet_offset(self, s, gp):
        """Each sub-cluster: gas (slowed) at 1, galaxies + collisionless dark matter at 0.
        The 4700 km/s collision is hugely supersonic, so the dark matter is in the
        normal phase and travels with the galaxies."""
        M_b = s.M_gas + s.M_star
        V = self.v_phonon(M_b / self.f_baryon(gp), gp)
        fs = self.f_super(V, gp, V_bulk=4700.0)
        dm = (1 - fs) * M_b * (1 / self.f_baryon(gp) - 1)
        gas_eff = s.M_gas * (1 + fs * (nu(G * M_b / 100.0 ** 2 / self.a0(gp)) - 1))
        return float(gas_eff / (gas_eff + s.M_star + dm))

    # ------------------------------------------------------------- dwarfs
    def breaking_radius(self, V_host, gp):
        """Radius inside which a host's condensate is dense enough to destroy the
        condensate of a satellite orbiting supersonically through it within an orbit.
        Scales like the host's condensate size, R ~ M^(1/5) ~ V^(4/5) (Berezhiani &
        Khoury 2015)."""
        return gp["R_break"] * (V_host / 200.0) ** 0.8

    def predict_dwarf_dispersion(self, d, gp):
        """Satellites inside the breaking radius move supersonically through the dense
        part of the host's condensate (Landau criterion): their own condensate is
        destroyed and they are ordinary cold-dark-matter halos (NFW, abundance
        matching). Satellites farther out keep their own condensate. A condensate is a
        separate phase-coherent domain, so the host's phonon field does not enter it:
        isolated deep-MOND-like dynamics, sigma^4 = (4/81) G M a0."""
        a0 = self.a0(gp)
        V_orbit = np.sqrt(d.g_ext * d.D_host)
        inside = 1.0 / (1.0 + (d.D_host / self.breaking_radius(V_orbit, gp)) ** 8)
        fs = self.f_super(self.v_phonon(d.M_bar, gp), gp, V_bulk=V_orbit * inside)
        s_sf = (4.0 / 81.0 * G * d.M_bar * a0) ** 0.25
        m200 = halo_mass_from_stars(d.M_bar)
        s_n = dispersion_from_mass(d.M_bar / 2 + nfw_mass_enclosed(d.r_half, m200, nfw_concentration(m200)),
                                   d.r_half)
        return float(np.sqrt(fs * s_sf ** 2 + (1 - fs) * s_n ** 2))

    def predict_wide_binary_boost(self, gN_si, g_ext_si, gp):
        """Binaries sit inside the Milky Way's superfluid: the phonon field responds to the
        total acceleration (external field effect)."""
        a0_si = gp["a0_e10"] * 1e-10
        return float(nu(np.hypot(gN_si, g_ext_si) / a0_si))

    # ------------------------------------------------------------- lensing
    def predict_lensing_vcirc(self, M_bar, R_kpc, gp):
        R = np.asarray(R_kpc, dtype=float)
        gN = G * M_bar / R ** 2
        fs = self.f_super(self.v_phonon(M_bar, gp), gp)
        v2 = gN * R * (1 + fs * (nu(gN / self.a0(gp)) - 1))
        if fs < 0.999:
            m200 = halo_mass_from_stars(M_bar)
            v2 = v2 + (1 - fs) * G * nfw_mass_enclosed(R, m200, nfw_concentration(m200)) / R
        return np.sqrt(v2)

    def predict_lens_mass(self, L, gp):
        """Stars as a Hernquist sphere matched to the total and the projected-in-r_Ein
        stellar masses; phonon force on the superfluid share; the effective (phantom)
        density is projected along the line of sight."""
        Mt = gp["imf_ell"] * L.Mstar_chab
        a = _hern_a(L.Mstar_proj_chab / L.Mstar_chab, L.r_Ein)
        a0 = self.a0(gp)
        fs = self.f_super(self.v_phonon(Mt, gp), gp)

        def Menc(r):
            Mn = Mt * r ** 2 / (r + a) ** 2
            return Mn * (1 + fs * (nu(G * Mn / r ** 2 / a0) - 1))
        M = _projected(Menc, L.r_Ein)
        if fs < 0.999:
            m200 = halo_mass_from_stars(L.Mstar_chab)
            M += (1 - fs) * nfw_mass_projected(L.r_Ein, m200, nfw_concentration(m200))
        return M

    # ------------------------------------------------------------- particle
    def dm_particle(self, gp):
        """An ~eV boson condensate, produced non-thermally (cold): no nuclear coupling
        and no annihilation into Standard Model particles at detectable rates."""
        return {"type": "superfluid", "mass_GeV": 1e-9, "sigma_SI_cm2": 0.0,
                "sigma_v_cm3s": 0.0, "wdm_mass_keV": float("inf")}


MODEL = TwoPhaseDM()
