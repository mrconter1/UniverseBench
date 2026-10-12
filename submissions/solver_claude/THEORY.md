# Collapsing self-interacting dark matter (cSIDM)

*Submitted by Claude (the benchmark's builder) acting as solver. Read
PREREGISTRATION.md first: the author had seen all the data, so the
hidden-split score is a check against overfitting, not a clean test.*

## The idea in plain language

Dark matter is a real, cold particle, as cosmology, galaxy clusters and the
Bullet Cluster require. It differs from standard cold dark matter in one
respect: dark-matter particles bounce off each other, and the bounce gets
much weaker as they move faster, as for particles that interact through a
light force carrier.

That single property gives three regimes, set only by how fast the dark
matter moves in each system:

* **Galaxies** (50–300 km/s). Particles scatter a few times over the
  galaxy's life. That heats the inner halo into a smooth, isothermal state
  that settles into the gravity of the stars and gas. The dark matter
  therefore follows where the ordinary matter sits, which is why rotation
  curves can be predicted from visible matter alone (the regularity MOND
  is built on) without any new force.
* **Small dwarf galaxies** (10–30 km/s). Scattering is so frequent that the
  halo has gone further. After forming a core, heat flows outward and the
  centre contracts again (gravothermal collapse), so the halos of dwarf
  spheroidals are dense and cuspy, as the modern data show.
* **Clusters** (about 1000 km/s). Scattering is negligible. The dark matter is
  collisionless, so cluster masses, the Bullet Cluster and all of cosmology
  are those of standard cold dark matter.

## Equations

Cross-section (Rutherford-like, light mediator):

    sigma/m (v) = sigma0 / (1 + (v/w)^2)^2

Each halo starts as NFW (mass from the baryons, concentration from
Dutton & Maccio 2014 shifted by a global offset). Its effective
cross-section is sigma(v) averaged with weight v^5 over a Maxwellian with
1D dispersion 0.64 V_max (Yang & Yu 2022).

*Core phase* (Kaplinghat et al. 2014): inside r1, where
rho_NFW(r1) (sigma/m) <v> t_age = 1, the dark matter is isothermal in the
total potential:

    (1/r^2) d/dr (r^2 d ln rho/dr) = -4 pi G (rho + rho_baryon) / sigma0^2,

with rho0 and sigma0 fixed by continuity of density and mass at r1 (the
least dense solution).

*Collapse clock* (Balberg, Shapiro & Inagaki 2002; Essig et al. 2019):

    t_c = (150/C) / (sigma_eff rho_s r_s) / sqrt(4 pi G rho_s),  C = 0.75,  tau = t_age / t_c.

The enclosed dark mass is blended in log space from the core solution
(tau <= 0.2) to the NFW cusp (tau >= 1), standing in for the collapse
phase. Where no cored solution exists (the scattering region would extend
beyond about 2 r_s), the halo is past its core phase and is treated as NFW.

*Galaxy–halo connection*: M_halo = 10^A (M_baryon / 1e10 Msun)^B + M_floor.
M_floor is the smallest halo that hosts a galaxy (reionisation and cooling).

How each test is predicted:

| Test | Prediction |
|---|---|
| 01, 02, 04 rotation curves | baryons + SIDM halo from the galaxy's baryonic mass; 02 adds one stellar mass-to-light knob per galaxy |
| 05 dwarfs | Plummer stars + SIDM halo, mostly collapsed; Wolf et al. 2010 estimator |
| 07 galaxy lensing | the same halos, outer NFW (scattering irrelevant at 50–1000 kpc) |
| 08 ellipticals | stars × IMF factor + projected NFW (sigma/m < 1 cm²/g at their speeds) |
| 03, 09 clusters | cold dark matter with a gas fraction that depends on mass and radius (as ΛCDM) |
| cosmology, guard rails | ΛCDM engine; general relativity |
| particle tests | cold, no Standard Model coupling required: no direct or indirect signal |

## Constants

Frozen (fitted by the author on the **public half** of the pre-registered
split: rotation curves, dwarfs and galaxy lensing jointly, each test's mean
chi² weighted equally, as the benchmark does):

| Constant | Value | Status |
|---|---|---|
| sigma0 | 128 cm²/g | **new** |
| w (turnover speed) | 151 km/s | **new** |
| M_floor | 2.7×10⁹ Msun | astrophysics |
| A, B (halo connection) | 11.60, 0.884 | astrophysics |
| concentration offset | −0.12 dex | astrophysics |
| disc mass-to-light (3.6 µm) | 0.463 | standard nuisance |
| halo age | 10 Gyr, collapse-blend range tau 0.2–1 | fixed choices |

Fitted by the benchmark: IMF factor for ellipticals, 5 cosmology constants,
3 cluster gas-fraction constants. **Total: 17 fitted constants** (the
leaderboard's "9 global" counts only those the benchmark fits).

Resulting cross-section: about 125 cm²/g at 10 km/s, 70 at 50 km/s, 20 at
150 km/s, 0.05 at 1000 km/s.

## Development record (public data only)

Rotation curves, 77 public galaxies, 5-fold cross-validation:

| Model | CV within 2σ | CV mean z² |
|---|---|---|
| MOND (fitted) | 65.9% | 6.46 |
| CDM with feedback cores (coreNFW) | 59.8% | 7.54 |
| SIDM, constant cross-section | 63.4% | 6.29 |

Constant-cross-section SIDM ties MOND (galaxy-by-galaxy difference −1.0
points, 95% CI −5.6 to +3.8), with no per-galaxy freedom. It fails the
dwarfs, though (cores make them far too cold: 47%). The speed-dependent
version with collapse fixes the dwarfs (75%) at a cost of a few points on
rotation curves. Fitting lensing jointly costs a few more: lensing out to
1 Mpc wants heavier halos than rotation curves do.

Things that did not work: the Di Cintio 2014 feedback profile (56.7% CV), a
gas-extent correction (fitted to 1, no effect), and a full two-halo term for
lensing (would overshoot for galaxies selected as isolated).

## Scores

See SCORES.md (practice split and the one-shot pre-registered hidden split).

## Where it fails or is weak

* **Rotation curves (about 60%).** Below MOND (about 65%). The model's own
  acceleration relation is tight (0.09 dex scatter) but its scale comes out
  at about 0.8×10⁻¹⁰ m/s². That is set by the fitted halo masses, so the
  hoped-for derivation of a0 from cosmology has **not** been achieved.
  Doing it would mean setting halo masses by abundance matching against the
  model's own halo mass function.
* **Lensing versus rotation curves.** The joint fit is a compromise.
* **Dwarf diversity.** Dense satellites (Ursa Minor, Leo II) are
  underpredicted and puffy, possibly tidally affected ones (And XXI, And II)
  overpredicted. Concentration scatter and tides are not modelled.
* **Collapse phase.** Approximated by a blend toward NFW. Real collapsed
  cores are denser than NFW; a calibrated gravothermal model would be better.
* **Clusters** are handled with the same gas-fraction recipe as the ΛCDM
  baseline. The model adds nothing new there.
* **Dark energy** is a cosmological constant. Not explained.

## Positions on the open questions

* **Wide binaries:** Newtonian (no new force). Sides with Banik et al.
* **DF2/DF4:** these are tidally stripped satellites; the model has no tidal
  physics, so it fails them, like ΛCDM.
* **Hubble tension, S8, lithium:** standard cosmology; no resolution.
* **Lyman-alpha, JWST, early black holes:** cold dark matter on large scales,
  same as ΛCDM.

## New predictions (checkable with future data)

1. **Isolated ultra-faint dwarfs** (e.g. found by Rubin far from any host)
   have dispersions of about **6–8 km/s** for 10³–10⁵ Msun of stars
   *(post-freeze correction, 12 Oct: this range was too generic; for the
   isolated dwarf Tucana B the model predicts 4.5 km/s (3.2–5.4) against
   MOND's 3.0 (2.4–3.6); see analysis/RESULTS.md. The model file is unchanged)*,
   because their halos sit at the halo-mass floor and have collapsed. MOND
   predicts 2–3 km/s; this is the cleanest test of the two.
2. **Core size grows with halo speed up to about 150 km/s, then shrinks.**
   Galaxies near V_max ≈ 50–150 km/s have the largest cores relative to their
   size. Above about 300 km/s, halos are nearly cuspy apart from baryonic
   contraction.
3. **Scatter in the galaxy relation correlates with halo concentration.** At
   fixed baryonic mass, galaxies in more concentrated halos (earlier
   forming, denser environments) lie systematically above the mean
   relation. A fundamental force law predicts no such correlation.
4. **Cluster cores:** sigma/m ≈ 0.05 cm²/g at 1000 km/s gives cluster
   cores of at most a few kpc. A robust measurement of larger cluster cores
   (tens of kpc) would rule out this cross-section curve.
5. **Satellites closer to their host are, if anything, denser** (tides speed
   up collapse). Not modelled here, but a qualitative prediction that
   cold dark matter (tides make halos lighter) does not share.
