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

First pass verdict: not supported, not ruled out; needed per-galaxy distance
and inclination uncertainties.

### A, decisive version with SPARC Table 1

Script: `analysis/rar_age_v2.py`. Each galaxy's offset gets its own error
budget from Table 1: distance (median 0.060 dex), inclination (0.028), stellar
M/L (0.1 dex on the stellar share) and statistics. The intrinsic
between-galaxy scatter τ is a maximum-likelihood fit. 156 galaxies with
inclination ≥ 30°.

| Split | Less evolved τ | More evolved τ | Difference (95% CI) |
|---|---|---|---|
| All galaxies | | | τ = 0.152 dex overall |
| Baryonic mass | 0.210 | 0.011 | +0.19 (+0.12..+0.25) |
| Hubble type | 0.200 | 0.068 | +0.13 (+0.05..+0.21) |
| Surface brightness | 0.185 | 0.081 | +0.10 (+0.03..+0.18) |
| Gas fraction | 0.171 | 0.129 | +0.05 (−0.05..+0.16) |
| **Quality-1 galaxies only** | | | **τ = 0.058 dex overall** |
| Q1: baryonic mass | 0.078 | 0.000 | +0.07 (0.00..+0.11) |
| Q1: Hubble type | 0.064 | 0.043 | +0.02 (−0.05..+0.09) |
| Q1: surface brightness | 0.065 | 0.049 | +0.02 (−0.06..+0.09) |
| Q1: gas fraction | 0.070 | 0.031 | +0.04 (−0.03..+0.10) |

Low-mass galaxies are only half quality-1 (38 of 78), against 55 of 78 for
massive ones. The signal in the full sample is driven by kinematic data
quality (non-circular motions, pressure support, asymmetries) that the error
bars do not capture. Among the best data the relation is very tight (about
0.06 dex intrinsic) and depends at most marginally on evolutionary state.

**Verdict: not supported.** The tightness itself is a hard target for any
emergent explanation: a model must produce ≤ 0.06 dex of intrinsic scatter
for well-measured galaxies, whatever their history.

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

## Idea 1: can collapsing SIDM be as tight as the real relation? (12 Oct 2026)

Script: `analysis/sidm_scatter.py`. The frozen cSIDM model on the 93 SPARC
quality-1 galaxies (inclination ≥ 30°), with realistic, uncorrelated halo
diversity: 0.15 dex in halo mass at fixed baryonic mass, 0.11 dex in
concentration. Scatter of per-galaxy offsets from one RAR fitted to the
model:

| Model | Between-galaxy scatter |
|---|---|
| Observed (intrinsic, after errors) | **0.058 dex** |
| cSIDM, every galaxy gets the average halo | 0.056 |
| cSIDM, realistic halo scatter (3 draws) | **0.102** (0.093–0.111) |
| cSIDM, concentration scatter only | 0.096 |
| Plain NFW, realistic halo scatter | 0.103 |

**Verdict: as built, cSIDM fails the tightness test.** Self-interactions do
not erase halo diversity here: with realistic scatter the model scatters as
much as plain cold dark matter, nearly twice what is observed. Its
benchmark-level agreement relied on assigning every galaxy the mean halo.
Escape routes, to be tested rather than assumed: halo scatter smaller for disc
galaxies than for halos in general, or scatter correlated with the baryons
(compact discs in concentrated halos), which some galaxy-formation
simulations report. MOND predicts the observed tightness with no extra
assumption, and this is now the strongest single point in its favour from
this project.

## Idea 2: is there an isolated ultra-faint dwarf with measured kinematics? (12 Oct 2026)

Literature search: **no isolated ultra-faint dwarf has a published velocity
dispersion yet.** The best target is **Tucana B** (Sand et al. 2022, ApJL 935,
L17). It is quenched and isolated (nearest galaxy ~500 kpc), at D ≈ 1.4 Mpc,
with L_V = 5(+4/−2)×10⁴ L☉ and R_half = 80 ± 40 pc. No spectroscopy has been
reported. The ultra-faint dwarfs used against MOND so far (e.g. "The Challenge
to MOND from Ultra-faint Dwarf Galaxies", ApJL 2021) are Milky Way
satellites, where the external field effect complicates the test.

**Pre-registered predictions** (script `analysis/tucana_b.py`; stellar
M/L_V 1.5–2.5):

| Picture | σ, central values (km/s) | range over L, M/L, R_half |
|---|---|---|
| MOND, isolated deep regime | 3.0 | 2.4–3.6 |
| cSIDM (frozen model) | 4.5 | 3.2–5.4 |
| CDM, 10⁹ M☉ NFW halo | 4.9 | 3.6–5.8 |
| Stars alone, Newtonian | 0.8 | 0.4–1.7 |

Correction to earlier statements (and to THEORY.md prediction 1, "6–8 km/s"):
for Tucana B specifically, dark matter predicts about 4.5–5 km/s against
MOND's 3. A ~1 km/s measurement with multi-epoch binary correction decides it
at the central size. A sharper half-light radius (space imaging) would remove
the 3.2–3.6 km/s overlap.

## Idea 3 (candidate O): MOND + ~11 eV sterile neutrinos (12 Oct 2026)

Model: `candidates/mond_hdm.py`. Published-a0 MOND in galaxies (the
Tremaine–Gunn phase-space bound keeps 11 eV fermions out of galaxies), plus
neutrino mass in clusters, η × M_baryon with mass and radius trends (three
knobs, like the ΛCDM baseline's baryon fraction), travelling with the
galaxies in the Bullet Cluster. Cosmology: ω_ν = m/94 eV mimics CDM in the
expansion history. A first version filled clusters at the Tremaine–Gunn cap
with a self-consistent σ; it runs away and was dropped.

Five demo splits, **computable tests only** (structure growth 21 and 31 and
extended cosmology 28–30, 34 excluded):

| | MOND + 11 eV ν | ΛCDM (same tests) | MOND alone |
|---|---|---|---|
| Lowest | **41** | 37 | 0 |
| Substantive passed | 7 | 7 | 3 |
| Median | **96** | 93 | 53 |
| Clusters / Bullet | 81 / 100 | 96 / 100 | 11 / 0 |
| Fitted m_ν | 11.1 eV | — | — |
| Neutrino / baryon mass in clusters (η) | ≈3.5 | — | — |

**Verdict: survives everything this benchmark can compute, which is exactly
what makes it dangerous to over-read.** The tests it skips are where it is
expected to fail. Structure growth (σ8, growth rate, cluster counts) needs
MOND N-body simulations; published simulations of this model have reportedly
struggled with the cluster population. The Lyman-α warm-dark-matter bound
(keV scale) excludes 11 eV particles as the small-scale dark matter. And
Tremaine–Gunn may be violated in cluster cores at small radii.
**Benchmark gap:** structure growth for non-CDM models needs either full
Boltzmann-code power spectra or a simulation-calibrated recipe.

## Idea 4 (candidate D): gravity with memory vs the Bullet Cluster (12 Oct 2026)

Order-of-magnitude test, no model fitting. If gravity responds to where
matter *was*, weighted over a memory time τ, the lensing peak of a moving
sub-cluster trails its galaxies by about v τ. Inputs from
`data/bullet_cluster/randall2008_offset.csv` (verified against Randall et al.
2008, arXiv:0704.0261): galaxy centroid 25 ± 29 kpc from the mass peak,
bullet velocity 4700 ± 630 km/s. Script: `analysis/bullet_memory.py`.

The test traps the idea from both sides:
* **long memory** (Gyr, needed to tie a0 to cH0): the lensing source sits
  where the matter was over the past Gyr, near the collision centre, i.e.
  hundreds of kpc behind the galaxies, near the gas. Contradicted by the
  observed alignment of mass and galaxies;
* **short memory:** keeping the peak within 83 kpc (2σ) of the galaxies
  needs τ < 17 Myr at 4700 km/s (< 20 Myr at 4070 km/s; an earlier version
  assumed 3000 km/s and got < 30 Myr). But then gravity effectively follows
  the matter's present position, which places the lensing peak on the gas
  (most of the baryons), the class already excluded.

**Verdict: ruled out for any memory time.** Mashhoon's nonlocal gravity
reduces in its Newtonian limit, as understood here, to a spatial
(instantaneous) kernel, so it falls in the present-matter class and does not
escape either. This removes the last no-dark-matter escape route from the
Bullet filter, leaving only theories with real collisionless mass (dark
matter of any temperature, MOND + hot DM, dipolar medium if it separates
from the gas).

## Idea 5: which well-measured galaxies break the relation? (12 Oct 2026)

Script: `analysis/rar_outliers.py` (93 SPARC quality-1 galaxies, same
per-galaxy error budget as rar_age_v2.py).

* **Mundane outliers:** three of the five largest outliers (UGC06667 +4.3σ,
  UGC08286 +2.5σ, UGC07125 −2.4σ) are edge-on (89–90°), where integrating
  through the disc biases rotation curves. Notable non-edge-on outliers:
  D631-7 (−3.8σ, gas-dominated dwarf, TRGB distance) and NGC0024 (+3.3σ).
* **By distance method:** galaxies with Hubble-flow distances (n=56) show
  intrinsic τ = 0.000 dex (their large distance errors absorb everything);
  galaxies with direct distances (TRGB, Cepheid, Ursa Major cluster,
  supernova; n=37) show **τ = 0.070 dex**.
* 3% of quality-1 galaxies are beyond 3σ (0.3% expected for pure noise).

**Interpretation:** with the most accurate distances the relation is *not*
scatter-free. About 0.07 dex remains, unless stellar M/L varies by more than
the 0.1 dex allowed here, or thick-disc and non-circular-motion corrections
are larger than assumed. The realistic target for any theory is 0.06–0.07 dex,
not zero. Pure MOND needs that much extra galaxy-to-galaxy variation from
systematics; cSIDM with realistic halo scatter (0.10 dex, idea 1) is still
too broad, but less dramatically.

## Idea 6: can a0 be derived instead of fitted? (12 Oct 2026)

Script: `analysis/a0_derived.py`. The fitted galaxy–halo link of cSIDM (A, B,
floor) is replaced by **abundance matching with zero fitted constants**:
GAMA stellar mass function (Baldry et al. 2012 double Schechter) plus
Papastergis et al. 2012 gas fractions (M_b = M* + 1.4 M_HI), both now in
`data/abundance_matching/` with citations and verified against the papers, matched to the Sheth–Tormen halo mass function of the
benchmark's Planck-like cosmology (centrals, no scatter). The abundance-matched
link agrees with the fitted one within ~50% for M_b = 10⁸–10¹¹ M☉.

| Halo link / constants (all 165 SPARC galaxies) | Emergent a0 (10⁻¹⁰ m/s²) | Rotation curves within 2σ |
|---|---|---|
| Data (same M/L) | **1.06** | — |
| Fitted link (frozen cSIDM) | 0.86 | 57.9% |
| Abundance-matched link, **concentration offset still fitted (−0.12 dex)** | 0.95 | 56.9% |
| … cross-section ÷10 | 0.94 | 57.8% |
| … cross-section ×3 | 1.29 | 56.8% |
| … concentration −0.2 dex | 0.36 | 44.7% |
| … concentration +0.2 dex | 2.40 | 48.2% |
| **… standard ΛCDM concentrations: nothing fitted to galaxies** | **1.68** | 54.4% |

For comparison, cH0/2π = 1.04 and the standard MOND value is 1.2.

**Findings.**
1. **With nothing fitted to galaxies (abundance matching + standard ΛCDM
   concentrations) the emergent a0 is 1.68×10⁻¹⁰, about 60% above the
   1.06 measured in the same galaxies.** The 0.95 agreement needs the
   concentration offset fitted to rotation curves, so it is not a
   derivation. What this does show is that the *order of magnitude* of a0
   has an ordinary cosmological origin: halo densities are set by the cosmic
   density when halos formed. It does not explain the precise value or the
   2π in a0 ≈ cH0/2π; any acceleration built from cosmic densities lands
   near cH0 within a small factor.
   This is essentially a **cold-dark-matter result, not an SIDM one**: the
   frozen model's cores barely reach the measured radii, so the cross-section
   hardly matters (÷10 changes a0 by 1%). Once cores grow (×3), it starts to
   (a0 = 1.29).
2. **a0 is nearly independent of the self-interaction cross-section**
   (÷10 changes it by 1%) but **scales steeply with halo concentration**
   (about ×2.5 per +0.2 dex). This is the same lever that broke the
   tightness test (idea 1): 0.11 dex of concentration scatter maps into
   ~0.1 dex of relation scatter. In a dark-matter picture, the galaxy
   relation is as tight as halo concentrations are uniform. The open
   question becomes why disc galaxies' halos would be more uniform in
   concentration than simulations predict.


## Correlated halo–baryon scatter: does it rescue the tightness? (12 Oct 2026)

Script: `analysis/sidm_scatter_corr.py`. As idea 1, but each galaxy's
concentration offset (0.11 dex) is correlated, with coefficient ρ, with its
disc compactness: the residual of log R_disk (SPARC Table 1) at fixed
baryonic mass. ρ > 0 means compact discs sit in more concentrated halos (the
simulation claim). Halo-mass scatter 0.15 dex, independent. 93 quality-1
galaxies, two draws each.

| ρ | frozen cross-section | cross-section ×3 |
|---|---|---|
| −0.9 | 0.107 | 0.126 |
| −0.5 | 0.106 | 0.121 |
| 0 | 0.105 | 0.113 |
| +0.5 | 0.109 | 0.105 |
| +0.9 | 0.111 | 0.107 |
| +1.0 | 0.110 | 0.105 |
| +1.0, no halo-mass scatter | 0.101 | 0.096 |
| no halo scatter at all | 0.056 | 0.059 |
| **target (observed intrinsic)** | **0.058–0.070** | |

**Verdict: correlated scatter does not rescue it in this model.** Even a
perfect correlation with disc compactness leaves ~0.10 dex. With the larger
cross-section, the simulation-motivated sign helps slightly (0.113 → 0.105),
far from enough. Reaching the observed tightness needs concentration scatter
below roughly 0.05 dex for disc galaxies, about half the ΛCDM prediction.
Caveat: disc scale length is the only compactness proxy available;
correlation with formation history (not in SPARC) could behave differently,
and real simulations include feedback effects not in this model. **The
tightness of the galaxy relation remains a genuine clue** against the
simplest dark-matter accounts, including this one.

## Test 15 upgraded to full Planck spectra via CLASS (12 Oct 2026)

The scored CMB test now compares each model's CLASS-computed lensed TT/TE/EE
spectra with Planck 2018 plik-lite (616 decorrelated points, own likelihood
implementation reproducing the reference χ² = 589.28 exactly). The compressed
summary, which presupposes ΛCDM perturbations, is now a fit aid only.

| Model (split 0) | CMB spectra within 2σ | χ² (616 points) | old compressed test |
|---|---|---|---|
| ΛCDM baseline (fitted) | **96.3%** | ≈ 590 | 100 |
| MOND + 11 eV sterile neutrinos (fitted; m_ν → 11.04 eV) | **74.4%** | **≈ 2580** | 100 |
| Visible matter only (default constants) | 11.5% | ≈ 2×10⁶ | 0 |

**MOND + hot DM is excluded by the Planck spectra** (Δχ² ≈ 2000), although
the compressed summary had given it a perfect score. The worst misses are in
the high-ℓ TT peaks (up to 9.6σ), and τ is pushed 2.8σ above the low-ℓ EE
value. Caveat: the sterile neutrino is fully thermalised (ΔN_eff ≈ 1); a
colder, non-thermal variant would add less radiation and should be tested
before the family is closed. A Δχ² of ~2000 is far beyond what that tweak
typically buys. This is the ΛCDM bias in compressed CMB summaries that
motivated the upgrade.

Interacting dark energy (candidate B) is not expressible in vanilla CLASS; it
now makes no prediction for test 15 until an interacting-DE CLASS fork is
wired in (rather than silently being scored as ΛCDM).
