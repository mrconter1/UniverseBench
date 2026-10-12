"""
Action health check for test 37.

A model declares its action in `ACTION`:

  {"type": "GR+minimal"}
      General relativity plus minimally coupled matter fields (particle dark
      matter of any kind, a cosmological constant). Healthy by construction.

  {"type": "horndeski",
   "G2": "...", "G3": "...", "G4": "...", "G5": "...",     # SymPy strings in phi, X
   "backgrounds": [{"z": 0.0, "H": ..., "dH": ..., "phi": ..., "dphi": ...,
                    "ddphi": ..., "rho_plus_p_matter": ...}, ...]}
      A scalar-tensor theory in Horndeski form,
          L = G2(phi,X) - G3(phi,X) box(phi) + G4(phi,X) R + G4_X [...] + G5 ...,
      with X = -(1/2) g^{mu nu} d_mu phi d_nu phi (X = dphi^2/2 on the
      expanding background), in units 8 pi G = 1 (GR is G4 = 1/2), together
      with background snapshots computed by the model's own cosmology (H in
      units of H0, dot = d/dt in units of 1/H0).

Checks, at every snapshot, from the quadratic action of Kobayashi,
Yamaguchi & Yokoyama 2011 (arXiv:1105.5723), scalar sector with minimally
coupled matter removed via dH -> dH + (rho+p)_matter / 2:
  * tensor no-ghost      G_T > 0
  * tensor gradient      F_T > 0
  * GW speed             |c_T^2 - 1| < 1e-15  (GW170817; checked at z < 0.1)
  * scalar no-ghost      G_S > 0
  * scalar gradient      c_S^2 = F_S / G_S > 0
  * (warning only)       c_S^2 <= 1 (no superluminal scalar)
Not covered: non-minimal matter couplings, beyond-Horndeski and vector
theories, strong-coupling or UV questions. A model outside the covered
classes declares {"type": "other", "note": ...} and gets no prediction.
"""
import numpy as np


def _horndeski_conditions(spec):
    import sympy as sp
    phi, dphi, ddphi, H, dH = sp.symbols("phi dphi ddphi H dH", real=True)
    Xs = sp.Symbol("X", real=True)
    ns = {"phi": phi, "X": Xs, "exp": sp.exp, "log": sp.log, "sqrt": sp.sqrt}
    G = [sp.sympify(spec.get(k, "0"), locals=ns) for k in ("G2", "G3", "G4", "G5")]
    G2, G3, G4, G5 = G
    d = lambda f, *v: sp.diff(f, *v)
    K, G3X, G4X, G5X = G2, d(G3, Xs), d(G4, Xs), d(G5, Xs)
    X = Xs
    GT = 2 * (G4 - 2 * X * G4X - X * (H * dphi * G5X - d(G5, phi)))
    FT = 2 * (G4 - X * (ddphi * G5X + d(G5, phi)))
    Theta = (-dphi * X * G3X + 2 * H * G4 - 8 * H * X * G4X - 8 * H * X ** 2 * d(G4, Xs, 2)
             + dphi * d(G4, phi) + 2 * X * dphi * d(G4, phi, Xs)
             - H ** 2 * dphi * (5 * X * G5X + 2 * X ** 2 * d(G5, Xs, 2))
             + 2 * H * X * (3 * d(G5, phi) + 2 * X * d(G5, phi, Xs)))
    Sigma = (X * d(K, Xs) + 2 * X ** 2 * d(K, Xs, 2) + 12 * H * dphi * X * G3X
             + 6 * H * dphi * X ** 2 * d(G3, Xs, 2) - 2 * X * d(G3, phi) - 2 * X ** 2 * d(G3, phi, Xs)
             - 6 * H ** 2 * G4
             + 6 * (H ** 2 * (7 * X * G4X + 16 * X ** 2 * d(G4, Xs, 2) + 4 * X ** 3 * d(G4, Xs, 3))
                    - H * dphi * (d(G4, phi) + 5 * X * d(G4, phi, Xs) + 2 * X ** 2 * d(G4, phi, Xs, 2)))
             + 30 * H ** 3 * dphi * X * G5X + 26 * H ** 3 * dphi * X ** 2 * d(G5, Xs, 2)
             + 4 * H ** 3 * dphi * X ** 3 * d(G5, Xs, 3)
             - 6 * H ** 2 * X * (6 * d(G5, phi) + 9 * X * d(G5, phi, Xs) + 2 * X ** 2 * d(G5, phi, Xs, 2)))
    sub = {Xs: dphi ** 2 / 2}
    GT, FT, Theta, Sigma = [sp.simplify(e.subs(sub)) for e in (GT, FT, Theta, Sigma)]
    GS = Sigma * GT ** 2 / Theta ** 2 + 3 * GT
    Y = GT ** 2 / Theta                                   # F_S = (1/a) d(a Y)/dt - F_T = H Y + dY/dt - F_T
    dY = sp.diff(Y, phi) * dphi + sp.diff(Y, dphi) * ddphi + sp.diff(Y, H) * dH
    FS = H * Y + dY - FT
    syms = (phi, dphi, ddphi, H, dH)
    return {name: sp.lambdify(syms, e, "numpy") for name, e in
            {"G_T": GT, "F_T": FT, "G_S": GS, "F_S": FS}.items()}


def _background_equations(spec):
    """Friedmann constraint of the theory, derived symbolically from its
    minisuperspace action (flat FRW, lapse N, scale factor a, field phi(t)):
        L = N a^3 [ G2 - G3 box(phi) + G4(phi) R ],
        box(phi) = -(1/(N a^3)) d/dt (a^3 phi_dot / N),
        R = 6 (a_ddot/(a N^2) - a_dot N_dot/(a N^3) + a_dot^2/(a^2 N^2)),
    plus minimally coupled matter (dS_m/dN = -a^3 rho). Supported for
    G4 = G4(phi) and G5 = 0, the Horndeski class that survives GW170817.
    Returns a function of a snapshot giving the constraint residual."""
    import sympy as sp
    t = sp.Symbol("t")
    N, a, ph = [sp.Function(n)(t) for n in ("N", "a", "ph")]
    phs, Xs = sp.symbols("phi X", real=True)
    ns = {"phi": phs, "X": Xs, "exp": sp.exp, "log": sp.log, "sqrt": sp.sqrt}
    G2, G3, G4, G5 = [sp.sympify(spec.get(k, "0"), locals=ns) for k in ("G2", "G3", "G4", "G5")]
    if sp.diff(G4, Xs) != 0 or G5 != 0:
        return None
    X = ph.diff(t) ** 2 / (2 * N ** 2)
    sub = {phs: ph, Xs: X}
    box = -1 / (N * a ** 3) * sp.diff(a ** 3 * ph.diff(t) / N, t)
    R = 6 * (a.diff(t, 2) / (a * N ** 2) - a.diff(t) * N.diff(t) / (a * N ** 3)
             + a.diff(t) ** 2 / (a ** 2 * N ** 2))
    L = N * a ** 3 * (G2.subs(sub) - G3.subs(sub) * box + G4.subs(sub) * R)
    rho = sp.Symbol("rho")
    eqN = (sp.diff(L, N) - sp.diff(sp.diff(L, N.diff(t)), t)) - a ** 3 * rho
    Hs, dHs, dphs, ddphs = sp.symbols("H dH dphi ddphi")
    e = eqN.subs({N.diff(t, 2): 0, N.diff(t): 0}).subs(N, 1)
    e = e.subs({a.diff(t, 2): (dHs + Hs ** 2) * a, a.diff(t): Hs * a})
    e = e.subs({ph.diff(t, 2): ddphs, ph.diff(t): dphs, ph: phs})
    e = sp.simplify(e / a ** 3)
    return sp.lambdify((phs, dphs, ddphs, Hs, dHs, rho), e, "numpy")


def check(action):
    """Returns (status, report). status: 'healthy', 'unhealthy' or None
    (no checkable action). report: list of (snapshot z, check, value, ok)."""
    if not action:
        return None, []
    kind = action.get("type")
    if kind == "GR+minimal":
        return "healthy", [(None, "GR + minimally coupled matter", 0.0, True)]
    if kind != "horndeski":
        return None, [(None, f"not checkable: {action.get('note', kind)}", 0.0, False)]
    f = _horndeski_conditions(action)
    bgeq = _background_equations(action)
    report, ok_all = [], True
    for b in action.get("backgrounds", []):
        if bgeq is not None and "rho_matter" in b:
            r1 = float(bgeq(b["phi"], b["dphi"], b.get("ddphi", 0.0), b["H"], b["dH"], b["rho_matter"]))
            r1 /= 3 * b["H"] ** 2
            report.append((b["z"], "background satisfies the theory's Friedmann constraint", abs(r1), abs(r1) < 1e-3))
            ok_all &= abs(r1) < 1e-3
        rpm = b.get("rho_plus_p_matter", b.get("rho_matter", 0.0) + b.get("p_matter", 0.0))
        dH_scalar = b["dH"] + 0.5 * rpm                                 # remove minimally coupled matter
        args = (b["phi"], b["dphi"], b.get("ddphi", 0.0), b["H"], dH_scalar)
        v = {k: float(fn(*args)) for k, fn in f.items()}
        cT2 = v["F_T"] / v["G_T"] if v["G_T"] != 0 else np.nan
        cS2 = v["F_S"] / v["G_S"] if v["G_S"] != 0 else np.nan
        items = [("tensor no-ghost G_T > 0", v["G_T"], v["G_T"] > 0),
                 ("tensor gradient F_T > 0", v["F_T"], v["F_T"] > 0),
                 ("scalar no-ghost G_S > 0", v["G_S"], v["G_S"] > 0),
                 ("scalar gradient c_S^2 > 0", cS2, bool(cS2 > 0))]
        if b["z"] < 0.1:
            items.append(("GW speed |c_T^2-1| < 1e-15 (GW170817)", abs(cT2 - 1), bool(abs(cT2 - 1) < 1e-15)))
        for name, val, ok in items:
            report.append((b["z"], name, val, ok))
            ok_all &= ok
        report.append((b["z"], "warning: c_S^2 <= 1 (subluminal)", cS2, bool(cS2 <= 1 + 1e-12)))
    if not action.get("backgrounds"):
        return None, [(None, "no background snapshots supplied", 0.0, False)]
    return ("healthy" if ok_all else "unhealthy"), report
