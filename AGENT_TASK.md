# Task: explain the dark universe

About 95% of the universe's content is "dark": matter we cannot see and an
energy driving cosmic acceleration. The standard model of cosmology (ΛCDM)
explains most observations with cold dark matter plus a cosmological constant,
but strains in some places (galaxy-scale regularities, several tensions).
Modified-gravity ideas such as MOND explain galaxies well but fail clusters and
the early universe. Nothing yet works everywhere.

Your job: propose **one physical model** (new dark matter, modified gravity,
or something else) and implement it so that it predicts every measurement in
this benchmark, from visible matter and a small set of universal constants.

## How you are scored

* **27 scored tests** (settled data) spanning galaxies, clusters, precision
  gravity, black holes, the early universe and cosmic expansion. The score for
  each is the % of data points you predict within 2 sigma (a perfect theory
  scores about 95; 100% would mean you fitted the noise).
* You are summarized by your **lowest, median and highest** score across the
  scored tests and **ranked by your lowest**: a theory is only as good as its
  worst failure.
* **11 open questions** (disputed or self-contradictory data, e.g. the Hubble
  tension) are reported alongside but never scored or used to fit constants.
  Say which side your theory takes and why.
* Where datasets are large (rotation curves, clusters, dwarf galaxies, lensing
  ellipticals, supernovae), part is **hidden** and only used for final scoring.
  You see the public part. Small tables are fully visible.
* A test your model makes no prediction for scores 0 ("—").
* Each test shows a **noise ceiling**: what the true laws of nature would
  score given its error bars (about 95, with a wider range for small tests).
  A test counts as passed when you reach what a perfect theory gets 95% of
  the time. Ranking: lowest score, then substantive (non-guard-rail) tests passed,
  then median.
* **Guard rails** (◆: solar system, pulsars, gravitational waves, lab, black
  holes, clocks, inflation, neutrinos, curvature, standard sirens, particle
  searches) are gates: failing one sinks your lowest
  score, but they don't inflate the median.

See LEADERBOARD.md for the baselines (ΛCDM, two versions of MOND, visible
matter only, pure Newtonian gravity) and ROADMAP.md for what each test is.

## Rules

1. Write your model in `my_model.py`, ending with `MODEL = YourClass()`.
   Score it with `python score.py my_model.py`.
2. **One set of constants for everything.** `GLOBAL_PARAMS = {name: (start,
   lo, hi)}` are fitted once, jointly, on the public data of every scored test.
   `LOCAL_PARAMS` are per-galaxy knobs, allowed only in the rotation-curve
   extrapolation test (their defaults are used everywhere else).
3. Predict from **inputs only**: never use the observed answer (observed
   speeds, total masses, lensing masses) or look up objects by name.
4. Approach the task as a scientist would: you are trusted to do physics, not
   to exploit simplifications in the tests (the Bullet Cluster test, marked §,
   is deliberately simple; model where your mass actually sits).
5. Use the computer you have as wisely as you can: the practice scorer, the
   literature, simulations, anything that helps you find a better theory.

## What to hand in

1. `my_model.py`: the model, runnable with `python score.py my_model.py`.
2. `THEORY.md`, containing:
   - the physical idea in plain language (a few paragraphs a non-physicist
     could follow),
   - the equations (ideally a Lagrangian or field equations) and how each
     prediction method follows from them,
   - every constant, its fitted value, and whether it is new or standard,
   - your practice scores and where the theory still fails,
   - your position on each open question (which side, and why),
   - at least one **new prediction** that a future observation could confirm
     or rule out, stated precisely enough to be checked.

## Model interface: one method per kind of prediction

| Method | Tests | Returns |
|---|---|---|
| `predict(g, R, gp, lp)` | 01, 02, 04 | rotation speed (km/s) at radii R (kpc) for galaxy `g` (`Vgas`, `Vdisk`, `Vbul` = Newtonian speeds of the visible parts for mass-to-light 1; `ups_fixed` means masses are already real; `stellar_mass(ups)`) |
| `predict_cluster_mass(c, gp)` | 03 | total mass (Msun) inside `c.r_kpc`, given visible mass `c.M_bar`, `c.z`, `c.delta` |
| `predict_dwarf_dispersion(d, gp)` | 05, 10 | stellar velocity dispersion (km/s) given stellar mass `d.M_bar`, gas `d.M_gas` (mostly outside r_half), `d.r_half` (3D, kpc), host field `d.g_ext` ((km/s)²/kpc), `d.D_host`, `d.host` (MW, M31 or field) |
| `predict_wide_binary_boost(gN_si, g_ext_si, gp)` | 06 | g_observed / g_Newton for a wide binary star |
| `predict_lensing_vcirc(M_bar, R_kpc, gp)` | 07 | circular speed (km/s) around an isolated galaxy, 50–1000 kpc |
| `predict_lens_mass(L, gp)` | 08 | projected mass (Msun) inside `L.r_Ein` given stars `L.Mstar_proj_chab`, `L.Mstar_chab` |
| `predict_bullet_offset(s, gp)` | 09 | lensing-peak position between galaxies (0) and gas (1), given `s.M_gas`, `s.M_star` |
| `predict_gravity(quantity, gp)` | 11–14, 26, 27, 34 | value of a precision-gravity quantity (see `data/gravity/measurements.csv`); inherit `models.common.GeneralRelativity` for GR |
| `cosmology(gp)` | 15–23, 28–35 | an object with the interface of `ubench.cosmology.Cosmology` (distances, sound horizon, θ*, growth, power spectrum, σ8, halo mass function, BBN yields, curvature, r, f_NL, N_eff, neutrino masses). Subclass it to change the physics |
| `class_params(gp)` | 15 | dict of CLASS input parameters (classy) for your cosmology; `CosmologyMixin` provides the standard one (CDM, massive neutrino, CPL dark energy, curvature). Override it to express other physics CLASS supports (extra or sterile neutrinos, dark radiation, interacting dark matter–radiation, …). List the constants it uses in `CLASS_DEPENDS` |
| `predict_bbn(gp)` (optional) | 16 | dict of Yp, D_H, Li_H if you have no full cosmology |
| `dm_particle(gp)` | 22, 24, 25 | dict (`mass_GeV`, `sigma_SI_cm2`, `sigma_v_cm3s`, `wdm_mass_keV`) or `None` for "no dark matter particle" |
| optional overrides | 32, 35, 36 | `predict_number_count_dipole(gp)`, `predict_max_bh_mass(z, gp)`, `predict_antihelium_ratio(gp)` |

The baselines in `ubench/models/` implement all of these and are good
templates. Each test module in `ubench/tests/` documents its data and units.
