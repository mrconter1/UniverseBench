# Candidate tests: results log

## A. Is the galaxy relation a relaxation attractor? (12 Oct 2026)

Script: `analysis/rar_age.py` (all 165 SPARC galaxies, RAR with a0 = 1.2e-10,
M/L 0.5 disc / 0.7 bulge). Prediction if A is right: less evolved galaxies
(gas-rich, low surface brightness) scatter MORE around the relation than
evolved ones. A fundamental force law predicts no difference.

| Scatter measure | Split by | Less evolved half | More evolved half | Difference (95% CI) |
|---|---|---|---|---|
| Between galaxies (offsets) | surface brightness | 0.22 dex | 0.14 dex | −0.08 (−0.15..−0.02) |
| Between galaxies (offsets) | baryonic mass | 0.24 | 0.12 | −0.12 (−0.17..−0.06) |
| Between galaxies (offsets) | gas fraction | 0.21 | 0.17 | −0.04 (−0.11..+0.03) |
| Within galaxies (shape) | gas fraction | 0.00 | 0.07 | opposite sign |
| Within galaxies (shape) | surface brightness | 0.00 | 0.07 | opposite sign |

Between galaxies the trend goes the predicted way, but galaxy-wide distance and
inclination errors (not in our files, larger for faint dwarfs) could produce
it. Within galaxies, which is immune to those errors, less evolved galaxies
show **no** intrinsic scatter. The evolved ones show some, most likely from a
fixed stellar M/L in bulge+disc systems.

**Verdict: not supported, not ruled out.** Deciding test: marginalise
per-galaxy distance and inclination uncertainties (SPARC Table 1; Li et al.
2018 method) and repeat the between-galaxy comparison.

## I. Does an "age clock" make a0 rise with redshift? (12 Oct 2026)

Script: `analysis/a0_redshift.py`. The frozen cSIDM model on all SPARC baryon
distributions, halos placed at redshift 1 (galaxy–halo connection held
fixed, an assumption); effective a0 fitted to the model's own accelerations.

| Case | Effective a0 | Ratio to z = 0 |
|---|---|---|
| z = 0 reference | 0.86×10⁻¹⁰ m/s² | 1.00 |
| Clock only (halo age 10 → 4.3 Gyr) | 0.76×10⁻¹⁰ | 0.88 |
| Density only (r200 at ρ_crit(z), Dutton & Maccio c(z)) | 1.04×10⁻¹⁰ | 1.21 |
| Clock + density (redshift 1) | 0.98×10⁻¹⁰ | **1.14** |

The clock pushes a0 *down* (less time to scatter), the opposite of the rough
scaling argued beforehand. **Verdict: the clock idea fails as an explanation
of a rising a0.** It does give a falsifiable discriminator: a0(z=1)/a0(0) is
1.0 for MOND, about 1.1–1.2 for cSIDM, about 1.8 for horizon-type theories
(a0 ∝ H), and about 2.5 as measured by MUSE-DARK III. Nothing tested predicts
2.5, so systematics in the measurement are the first suspect.

**How robust is MUSE-DARK's ~2.5×?** (A&A 2026, doi:10.1051/0004-6361/202659230;
79 galaxies, 0.33 < z < 1.44.) Their linear fit gives a0(0) = 1.0 and a slope
of 1.59 per unit redshift (×10⁻¹⁰ m/s²), so about 2.6× at z = 1. Two
systematics both push the inferred a0 upward. Gas masses at z ~ 1 are not
measured directly (no HI that far), so they come from scaling relations. And
the stellar M/L adopted is lower than locally. Both lower g_bar and so raise
a0. The paper's intrinsic scatter (0.17 dex) is 1.5× SPARC's, and inclination
and gas-mass systematics are not discussed. Treat 2.5× as not yet decisive.
A direct check would refit their galaxies with gas fractions raised by the
z ~ 1 range from molecular-gas surveys.

## C (and E, F): can a MOND-family theory be rescued? (12 Oct 2026)

Model: `candidates/dipolar_dm.py`. A proxy for dipolar dark matter
(gravitational polarization): published-a0 MOND in galaxies and clusters,
cold-dark-matter cosmology. Five demo splits:

| | Dipolar DM proxy | MOND (published) | ΛCDM |
|---|---|---|---|
| Lowest | **0** (Bullet Cluster) | 0 | 37 |
| Substantive tests passed | 8 | 3 | 9 |
| Median | 82 | 26 | 90 |
| Clusters / Bullet Cluster | 11 / 0 | 11 / 0 | 96 / 100 |

A correct cosmology lifts MOND from 3 to 8 passes but cannot fix clusters or
the Bullet Cluster.

**Bullet Cluster filter (corrected wording).** A theory whose extra gravity
is a direct function of where the ordinary matter is *now* fails the Bullet
Cluster: the lensing mass follows the galaxies, not the gas that holds most
of the baryons. This applies to MOND, running G (E), horizon-set inertia (F)
and Verlinde-style emergent gravity, so E and F were not built separately.
It is a strong filter, not a theorem for every MOND-like theory. Escape
routes to test rather than count as casualties:
* **dipolar dark matter (C) proper:** its polarizable medium can in
  principle end up away from the gas after a collision. Only the MOND-plus-
  CDM-cosmology *proxy* was tested here, not where the medium goes;
* **MOND plus a hot dark component** (e.g. ~11 eV sterile neutrinos,
  Angus 2009), which carries real collisionless mass in clusters;
* **gravity with memory (D),** whose field follows where matter *was*.

## B. Dark matter interacting with dark energy (12 Oct 2026)

Model: `candidates/interacting_de.py`. Energy transfer Q = ξ H ρ_DE, w = −1,
exact background, drag-only growth approximation; ξ = 0 reproduces ΛCDM
exactly (checked: E(z), θ*, CMB priors, σ8). One extra constant. Galaxy and
cluster physics as the ΛCDM baseline. Five demo splits:

| | Interacting DE | ΛCDM |
|---|---|---|
| Fitted ξ | −0.072 (dark matter → dark energy) | 0 |
| BAO (DESI DR2) | **92** | 85 |
| Cluster counts | **100** | 75 |
| Microwave background, supernovae, growth rate | 100, 99, 100 | 100, 99, 100 |
| Lowest / passed / median | 37 / 9 / **95** | 37 / 9 / 90 |
| Open question S8 (lensing) | 33 | 67 |
| Open question Lyman-α | 33 | 67 |

**Verdict: mildly favoured, not significant.** One extra constant buys one
more BAO point inside 2σ and removes the cluster-count tension. Note that ξ
being identical on every split is **not** evidence: the cosmology tests are
in-sample (‡), so every split fits exactly the same data. The model also
sides with cluster counts against weak lensing on S8, worsening a tension
that is already a sore point for ΛCDM. Its significance (Δχ² against the
extra constant) was not assessed, and the growth equation is approximate.
