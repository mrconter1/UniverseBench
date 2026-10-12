# Two-phase superfluid dark matter

## The idea in plain language

Dark matter is one kind of very light particle (about an electron-volt), born cold in the early universe, that interacts weakly with itself. Like helium cooled to near absolute zero, it can exist in two states. When it is cold and slow, it condenses into a **superfluid**. When it is hot or forced to flow too fast, it is an ordinary **normal gas** of particles.

Inside a quiet, isolated galaxy the dark matter is cold enough to be a superfluid. Sound waves in that superfluid (phonons) couple weakly to ordinary matter. They carry an extra pull that is negligible where gravity is strong and takes over where it is weak. The rule it follows is Milgrom's MOND law: below an acceleration a0 ≈ 1.3×10⁻¹⁰ m/s², gravity falls off like 1/r instead of 1/r². That is why rotation curves stay flat and why the visible mass alone predicts a galaxy's speed. No per-galaxy halo needs to be tuned.

Galaxy clusters are much hotter: their dark matter moves at more than the critical speed of about 350 km/s. There the condensate melts, and the dark matter is ordinary collisionless cold dark matter, in the cosmic proportion to the gas. This explains why MOND fails in clusters, and why the Bullet Cluster's mass travels with its galaxies: the collision was supersonic, so the dark matter there is a normal collisionless gas.

A condensate can also be destroyed by motion. Helium stops being superfluid when it flows faster than a critical speed (Landau's criterion). A small satellite galaxy orbiting at ~200 km/s through the dense inner condensate of the Milky Way or Andromeda loses its own superfluid and becomes an ordinary cold-dark-matter dwarf. A satellite far out, beyond a "breaking radius", keeps its condensate and behaves like an isolated MOND galaxy.

On cosmological scales both phases are pressureless matter, so the early universe, the microwave background, BBN and the growth of structure are those of the standard model. Dark energy is a separate, slowly thawing field. The data (DESI BAO, cluster counts) prefer it to a pure cosmological constant; this theory does not explain it.

## Equations

Field content: gravity (GR), Standard Model, a complex scalar Φ of mass m ~ 1 eV with a repulsive self-interaction, and a dark-energy scalar. Below the condensation temperature, the low-energy theory of the condensate's phase θ is the Berezhiani–Khoury superfluid EFT:

  L = (2Λ(2m)^{3/2}/3) · X·√|X| − α (Λ/M_Pl) θ ρ_b,   X = μ̇ − mΦ_N − (∇θ)²/(2m)

The P ∝ X^{3/2} equation of state makes the static phonon equation a nonlinear Poisson (AQUAL) equation. Its solution adds a phonon acceleration, so the total acceleration is

  g = g_N ν(g_N / a0),  ν(y) = 1/(1 − e^{−√y}),  a0 ∝ α³Λ²/M_Pl.

The phonon also couples to light through a disformal metric (as in TeVeS/AeST), so lensing follows the same g. Outside the condensate, or above T_c, the field is a normal collisionless gas: plain GR plus cold dark matter.

**Phase rule** (one rule for every system). The superfluid fraction is

  f_s = 1 / [1 + (V/V_c)⁴ + (V_bulk/V)⁴],

where V = (G M_b a0)^{1/4} is the virial speed the condensate would have (set by the baryons alone), V_c is the critical speed, and V_bulk is the speed at which the system moves through an ambient condensate (Landau criterion). The total enclosed mass of any system is then

  M_tot(r) = f_s · M_b ν(g_N/a0) + (1 − f_s) · [normal-phase dark matter + baryons].

How each prediction follows:

| Test | Prediction |
|---|---|
| 01, 02, 04 rotation curves | v² = g_N r [1 + f_s(ν−1)] + (1−f_s) v²_NFW. Galaxies have V < V_c, so f_s ≈ 1 (MOND-like). One stellar M/L knob per galaxy only in test 02. |
| 03 clusters | Phonon share M_b ν plus normal share M_b / f_b(Δ). f_b = Y·Ω_b/Ω_m, with Ω_b/Ω_m from the CMB constants. Hot gas is more extended than dark matter, so f_b(Δ) ∝ (Δ/500)^{−s}. |
| 05 dwarfs | Inside the breaking radius R_b(V_host) = R_break·(V_host/200)^{0.8} (condensate size ∝ M^{1/5}): normal phase, NFW halo from abundance matching. Outside: isolated deep MOND, σ⁴ = (4/81) G M a0. The satellite's condensate is its own phase-coherent domain, so it feels no external-field effect. |
| 07 galaxy lensing | Phonon-mediated field of an isolated galaxy (f_s ≈ 1). |
| 08 SLACS lenses | Stars as a Hernquist sphere matched to the total and Einstein-radius stellar masses. The phonon field's effective density is projected along the line of sight (this matters: the phantom halo is isothermal), plus a (1−f_s) NFW share. |
| 09 Bullet Cluster | Collision at 4700 km/s ≫ sound speed, so f_s → 0. Collisionless dark matter in the cosmic ratio sits with the galaxies. |
| 06 wide binaries | No condensate of their own, inside the Milky Way's: phonon force with the external field, ν(√(g_N²+g_ext²)/a0). |
| 11–14, 26, 27, 34 | GR. Phonon effects vanish at g ≫ a0, and the condensate breaks down inside ~250 AU of the Sun (BK 2015). GW speed = c. |
| Cosmology (15–31) | ΛCDM engine with ω_c for both phases and CPL dark energy (w0, wa). |
| 22, 24, 25 particle | Non-thermal, cold eV-scale boson: no nuclear coupling, no annihilation; not warm dark matter. |

## Constants (fitted jointly on public data; values from the final hidden-split run)

| Constant | Value | New or standard |
|---|---|---|
| a0 | 1.32 ×10⁻¹⁰ m/s² | Standard (MOND scale), refitted |
| V_c (critical speed) | 10^2.546 ≈ 352 km/s | **New** |
| R_break (V=200 host) | 359 kpc (poorly determined; see "Where it still fails") | **New** |
| Y (baryon retention) | 0.85 | Astrophysical |
| s (gas-profile slope) | 0.116 | Astrophysical |
| Υ_disk (3.6 µm) | 0.44 | Standard nuisance |
| IMF factor, ellipticals | 2.0 (at bound; Salpeter = 1.78) | Standard nuisance |
| H0, ω_b, ω_c, n_s, ln10¹⁰A_s | 69.2, 0.02229, 0.1206, 0.965, 3.052 | Standard |
| w0, wa | −0.88, −0.77 | Dark-energy sector, DESI-like |
| Fixed: transition sharpness 4, R_b ∝ V^0.8, bulge M/L = 1.4× disc | — | Choices |

That is 14 fitted global constants and 1 per-galaxy knob (test 02 only).

## Scores

**Final, one-shot score on the secret hidden split** (agent bundle built with a random salt; developed only on public data):

| | Two-phase SfDM | ΛCDM | MOND (pub.) |
|---|---|---|---|
| Lowest | **33** | 33 | 0 |
| Substantive tests passed (of 14) | **10** | 9 | 3 |
| Median (excl. guard rails) | **100** | 92 | 33 |
| Rank | **1** | 2 | 3 |

**Leaderboard-comparable run** (`score.py --splits 5`, averaged over the 5 demo splits used for LEADERBOARD.md): lowest **47** (ΛCDM 37, MOND 0), substantive tests passed **10** (ΛCDM 9), median **100** (ΛCDM 90), rank **1**. Per test: 01: 64, 02: 81, 03: 76, 04: 100, 05: 47, 07: 100, 08: 92, 09: 100, 17: 85, 21: 100; everything else 99–100.

Per test (secret hidden split): 01: 58, 02: 76, 03: 58, 04: 100, 05: 33, 07: 100, 08: 95, 09: 100, 15: 100, 16: 100, 17: 85, 18: 99, 21: 100, 31: 100, and all guard rails 100. On practice splits the lowest score was 60, 62.5 and 69 (ΛCDM: 43).

## Where it still fails

- **Dwarf spheroidals (05: 33).** The joint fit put R_break at 359 kpc, so almost every satellite was treated as normal phase. With the other constants held at their fitted values, the public-dwarf χ² has two basins: about 2.8 near 60–100 kpc and about 8.6 for R ≳ 160 kpc, separated by a barrier (15 at 140 kpc). The benchmark's joint Powell fit nevertheless ends on the far plateau. A post-hoc rerun started at 70 kpc also went to the 400 kpc bound and gave identical scores, so the problem is the sharp, multi-modal landscape plus the joint line searches, not the starting point. Cross-validation on public dwarfs (1-D fit of R_break on random halves) gave 69% for the rule against 53% for "always normal", but that advantage did not survive the benchmark's own fit. A smoother transition, or deriving R_break from the particle mass instead of fitting it, is the obvious next step.
- **Cluster masses (03: 58).** The hidden clusters want a baryon fraction that rises with mass, more steeply than the superfluid admixture in groups provides. The rising gas fraction with radius is captured; the mass trend is not.
- **Rotation-curve prediction (01: 58).** This sits at MOND's level. The remaining misses are galaxy-wide offsets (distance, inclination) that no global law can fix.
- **BAO (17: 85).** Two DESI points (LRG1 D_M, LRG2 D_H) stay outside 2σ even with evolving dark energy.
- **Dark energy** is not explained, only parametrised.
- **The theory is phenomenological.** The phase rule is a smooth interpolation, not a computed phase diagram. The Hernquist and NFW profiles are approximations.

## Positions on the open questions

- **06 Wide binaries:** MOND-like with the external field effect (boost ~1.3 at the Sun). This sides with Chae. A clean Newtonian result (Banik+) would push the theory towards condensate breakdown at stellar scales, which BK's estimate (breakdown inside ~250 AU) does not reach.
- **10 DF2/DF4:** fails. Treated as satellites inside NGC 1052's condensate, they get normal-phase halos. Tidal stripping, or the larger distance, would be needed.
- **16b Lithium:** standard BBN, so it fails like everyone. Stellar depletion is the likely answer.
- **19 Hubble tension:** H0 = 69.2 from the early universe. Between SH0ES and CCHP; it matches CCHP.
- **20 S8:** S8 ≈ 0.84 sides with the CMB and eROSITA, against older lensing surveys.
- **22 Lyman-α:** cold (the condensate is not thermal WDM), so the WDM bound passes.
- **23 JWST:** passes (CDM-like halos).
- **29b Neutrino squeeze:** Σmν = 0.06 eV passes. Evolving dark energy relaxes the squeeze.
- **32 Dipole, 33 dark flow, 35 early black holes:** no special prediction; standard values.

## New predictions

1. **A phase transition in the radial acceleration relation of galaxy groups.** Systems whose condensate speed V = (G M_b a0)^{1/4} is below V_c ≈ 350 km/s (baryonic mass ≲ 10^{11.9} M☉) follow the MOND relation M_dyn = M_b ν(g_N/a0) with no extra dark matter. Hot groups and clusters follow M_dyn ≈ M_b/(0.85 Ω_b/Ω_m). The switch happens within a factor ~1.5 in V (between ~250 and ~500 km/s). Neither ΛCDM (smooth) nor MOND (no switch) predicts it. eROSITA groups with weak-lensing masses can test it.
2. **Satellite dichotomy without an external-field effect.** Dwarf satellites beyond the breaking radius obey σ⁴ = (4/81) G M_* a0, independent of half-light radius and distance to the host. Satellites well inside obey cold-dark-matter scaling, with σ growing with half-light radius. Rubin/LSST satellites with spectroscopy can test this. The breaking radius itself is poorly measured, so the test is the existence of a sharp radial switch.
3. **The Bullet Cluster's superfluid remnants.** Galaxy-scale condensates survive in cluster galaxies. Strong-lensing mass maps of individual cluster members should show MOND-like (isothermal, baryon-predicted) profiles, while the intracluster dark matter is collisionless NFW.
4. **No direct-detection or annihilation signal** at any experiment. The particle is an eV-scale boson.
