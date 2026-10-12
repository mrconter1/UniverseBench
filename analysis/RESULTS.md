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

## Idea 10 (user): one light field behind evolving dark energy and DM self-interactions (12 Oct 2026)

Script: `analysis/idea10_sigma_drift.py`. If the self-interaction strength
drifts as the field rolls, does the emergent a0 drift with it? Frozen cSIDM,
abundance-matched halo link:

| σ/m × | 0.1 | 0.3 | 0.5 | 1 | 2 | 3 | 10 |
|---|---|---|---|---|---|---|---|
| emergent a0 (10⁻¹⁰ m/s²) | 0.94 | 0.82 | 0.80 | 0.95 | 1.18 | 1.29 | 1.33 |

Near the fitted value, d ln a0 / d ln(σ/m) ≈ 0.2–0.3, and the response is not
monotonic. A 30% drift in the interaction strength shifts a0 by ~7%, which is
small against the halo-density effect (×1.14 by z = 1, idea I) and far below
what current a0(z) measurements resolve. **Verdict: the predicted signal
exists but is weak in this model; a0(z) is not a sharp test of it.** Sharper
consequences of such a field are a long-range fifth force between dark matter
particles (constrained by tidal streams, e.g. Sagittarius) and a DM–DE
energy exchange (candidate B).

## Idea 8: a heavier electron at recombination (Hubble tension) (12 Oct 2026)

Model: `candidates/varying_me.py`. ΛCDM plus m_e(early)/m_e(today), switching
at z = 50. Test 15 uses CLASS's varying-constants module exactly; the fast
engine rescales z* ∝ m_e and r_d ∝ m_e^−0.617 (calibrated with CLASS,
exponent −0.61 to −0.63 across cosmologies). Cosmology tests only, split 0.
**No local H0 data were used in the fit** (the Hubble tension is an open
question).

| | ΛCDM | varying m_e |
|---|---|---|
| fitted m_e ratio | — | **1.024** |
| H0 | 67.7 | **70.6** |
| CMB spectra within 2σ | 96.3% | 96.6% |
| BAO / supernovae | 92 / 99 | 92 / 99 |
| Open question 19 (SH0ES 73.0, CCHP 70.4) | 1 of 2 | 1 of 2 (now at CCHP) |
| S8 (open) | 67 | 0 (ω_c rises to 0.127) |
| **χ² CMB + BAO + SN** | **1310.9** | **1325.2** (CMB 613 vs 592, BAO 10.5 vs 19.6) |

**Verdict: no real preference. The benchmark's fitting rule creates one.**
Every test counts equally in the fit (mean χ² per point), so DESI's 13 BAO
numbers weigh as much as Planck's 616 spectral points, and the fit trades
21 units of CMB χ² for 9 of BAO. On the summed likelihood, a heavier early
electron is disfavoured (Δχ² = +14 with an extra constant). It remains a
working mechanism: it does raise H0 to the CCHP value while keeping the CMB
at 96.6%, and it pays for it in S8.

**Benchmark lesson:** per-test equal weighting is right for *ranking* theories
across very different observables, but it distorts *parameter inference*
when tests differ greatly in size. The cosmology constants should be fitted
on the summed χ² of the cosmology tests (a proper joint likelihood), with
equal weighting kept for scoring. **Changed later the same day:** the
cosmology data sets now form one summed χ² block in the fit
(`ubench/core.py`, `JOINT_IDS`).

## Idea 7 (user): fermionic keV dark matter with a quantum core instead of Sgr A*'s black hole (12 Oct 2026)

Script: `analysis/fermion_core.py`; data: `data/sgr_a/sgra_constraints.csv`
(GRAVITY 2018, 2020). For each fermion mass, the structure equations of
general relativity (TOV) for a fully degenerate spin-½ Fermi gas holding Sgr
A*'s 4.26 million solar masses. Zero temperature is the most compact case on
the stable branch, so this is the best case for the idea.

| Fermion mass | Core radius | Mass inside the flare orbit (8 r_g) | Flare period (obs. 45 ± 15 min) |
|---|---|---|---|
| 56 keV (proponents' favourite) | 115 AU (just inside S2's 120 AU pericentre) | 1.5 × 10⁻⁷ | 89 days (+8500σ) |
| 300 keV (Crespi+ 2025 "high compactness") | 1.1 AU | 14% | 131 min (+5.7σ) |
| 345 keV | 0.69 AU | 47% | 73 min (+1.8σ) |
| 360 keV | 0.58 AU | 65% | 62 min (+1.1σ) |
| ≥ 383 keV | — | — | core cannot exist: collapses to a black hole |

**Verdict: excluded except in a narrow window (~345–383 keV) where the core
sits at the edge of collapse.** S2 cannot tell (every core is inside its
orbit, which is why S2-based papers find the models indistinguishable), but
GRAVITY's flares orbiting at 6–10 gravitational radii need almost all the
mass inside ~8 r_g. No configuration has a photon sphere (the most compact
degenerate ball is still ~9 r_g across), so none casts a black-hole shadow;
reproducing the EHT image would need the accretion flow to mimic one.
Fine-tuning plus the missing shadow make this < 1%.

## Idea 9 (user): cosmologically coupled black holes as dark energy (12 Oct 2026)

Model: `candidates/coupled_bh.py`; data: `data/star_formation/madau_dickinson2014.csv`.
With coupling k = 3, each generation of stellar black holes adds a constant
physical density from its birth, so the dark-energy history follows from the
cosmic star-formation history with no free shape (Croker+ 2024). The black
holes' birth mass comes out of the baryons: with Madau–Dickinson, **46% of
baryons are converted by today** (Croker+ fit ~30% with a higher early
star-formation rate). Net effect: the dark energy grows with time, i.e.
phantom, w = −1.03 today, −1.21 at z = 1, −1.41 at z = 2. For CLASS, the
CPL (w0, wa) = (−1.003, −0.361) matching H(z) to 0.3%.

Cosmology tests, split 0, same constants as ΛCDM (old equal-weight fit):

| | ΛCDM | Coupled BHs |
|---|---|---|
| H0 | 67.7 | 69.7 |
| Cluster counts | 75 | 100 |
| χ² CMB / BAO / SN | 592.3 / 19.6 / 698.9 | 596.4 / 18.7 / 711.6 |
| **χ² sum** | **1310.9** | **1326.6** |

**Verdict: disfavoured** (Δχ² = +16 at equal constants, mostly supernovae,
which dislike phantom dark energy), on top of gravitational-wave population
bounds k < 1.1–2.1 at 2σ (arXiv:2307.02474) and the need for a large
missing-baryon fraction.

## Inverse design: what thermal state must dark matter have in each galaxy? (12 Oct 2026)

Scripts: `analysis/inverse_design.py`, `inverse_design_laws.py`,
`inverse_design_laws2.py`, `inverse_eos.py` (87 quality-1 SPARC galaxies).

1. **No dark matter with one universal equation of state can make the
   galaxy relation.** Scaling argument: a self-gravitating P = Kρ^γ halo
   gives v⁴ ∝ M only for γ = 0. With data: requiring the RAR exactly, the
   dark-matter "temperature" P/ρ each galaxy needs scatters by 0.45 dex at
   fixed density. A local coupling to the baryonic field alone: 0.36 dex.
2. **What the data want is a per-galaxy thermostat:** nearly isothermal
   (0.10 dex variation within a galaxy), temperature set by the galaxy's
   baryonic mass, σ² ∝ M_b^0.57 (MOND: 0.5), 0.073 dex between galaxies
   (0.051 adding the local field). Caveat: the outer boundary condition
   (isothermal tail) builds part of this in.
3. **The thermostat alone fails.** Isothermal halos in each galaxy's own
   baryonic potential with σ⁴ = G M_b a0/4 and simple rules for the central
   density (constant, power law, core at the MOND radius, core ∝ disc size,
   no core) give χ² 25,000–37,000 against MOND's 9,200 with no free
   constants, and 0.10–0.20 dex scatter. The halo must also respond to the
   baryons locally in its inner parts.

## Is the tightness of the galaxy relation really a problem for cold dark matter? (12 Oct 2026)

Course check: the leading hypothesis (cold dark matter + galaxy formation)
had never faced the scatter test. Literature (all values idealised, no mock
observational errors added to the simulations):

| Simulation | Intrinsic RAR scatter |
|---|---|
| MUGS2 (Keller & Wadsley 2017) | 0.06 dex |
| ZOMG centrals (Garaldi+ 2018) | 0.048 (satellites 0.19) |
| NIHAO (Dutton+ 2019) | 0.079 overall; ≈ 0.04 for M* 10⁹–10¹¹, ≈ 0.11 below |
| EAGLE/APOSTLE (Ludlow+ 2017) | 0.09, but a0 ≈ 3 × 10⁻¹⁰ (2.5× too high) |
| MassiveBlack-II (Tenneti+ 2018) | ≈ 0.1 |
| Semi-analytic (Paranjape & Sheth 2021) | 0.06 spirals; 0.03–0.04 at fixed halo mass and concentration |

Observed intrinsic estimates: 0.06 (Dutton+ 2019, 90% range 0–0.09), 0.034 ±
0.014 (Desmond 2023, → 0 with a 10% baryonic systematic). Ours, robustness
(`analysis/rar_scatter_robustness.py`, same estimator, one assumption at a time):

| Variation | τ (dex) |
|---|---|
| Baseline, quality 1 (68% bootstrap) | 0.058 (0.039–0.069) |
| Distance errors ×0.5 / ×1.5 | 0.081 / 0.047 |
| M/L uncertainty 0.05 / 0.15 dex | 0.068 / 0.050 |
| +5% / +10% baryonic-velocity systematic | 0.052 / 0.025 |
| a0 1.0 or 1.4; M/L 0.4 or 0.6 | 0.057–0.059 |
| Quality 1 with direct distances (n = 37) | 0.070 (0.047–0.087) |
| Quality 2 only | 0.169 |

Our own halo models, same galaxies, noise-free (`analysis/sidm_scatter.py`):
NFW halos from abundance matching give 0.068 between galaxies with no halo
scatter and ~0.10 with 0.15 dex halo-mass + 0.11 dex concentration scatter;
cSIDM is the same (0.10). Pooled point-by-point widths: NFW 0.11, cSIDM 0.13.

**Verdict: the scatter is not (yet) a discriminator.** The observed intrinsic
scatter is uncertain by about a factor of two (0.03–0.08 dex under
reasonable error models), and full hydrodynamic CDM simulations land at
0.04–0.09 because galaxy formation narrows the halos a given galaxy can live
in. What our simple halo models lacked was that narrowing, not new physics.
The sharper CDM problems in this literature are the **normalisation and shape**
(EAGLE's a0 2.5× too high, MB-II preferring a power law), the absence of
residual correlations with galaxy size in the data (Stiskalek & Desmond
2023), and dwarfs. The decisive missing analysis is a like-for-like one:
mock-observed simulated galaxies (SPARC-like errors and selection) through
the same estimator.

**Lensing (Brouwer+ 2021, KiDS-1000):** the relation continues MOND-like to
g_bar ~ 10⁻¹⁵ m/s² for isolated lenses; MICE (dark-matter-only plus halo
occupation) also reproduces it, BAHAMAS (hydro) is offset; early- and
late-type galaxies differ at ≥ 6σ, which a single force law struggles with
and CDM expects. Not a clean MOND-vs-CDM discriminator either.

## New model: a dark force tied to dark energy (12 Oct 2026), stopped on paper

Proposal: the field driving dark energy exerts a k-mouflage-type force with
a nonlinear kinetic term whose scale is set by its cosmological rolling, so
an acceleration scale appears, a0 ≈ β √(3(1+w)Ω_DE) cH0 (β ≈ 0.25 for
w ≈ −0.75). Before building it, the inverse-design result was applied:
* **Force on dark matter only:** the baryons still reach the dark matter only
  through ordinary gravity, so the halo cannot know each galaxy's baryonic
  mass any better than a cold or self-interacting halo does. Same class as
  cSIDM: it can reshape halos but not supply the per-galaxy link.
* **Force on baryons too, conformal coupling:** a scalar force of this kind
  does not bend light, so lensing would see only the dark matter while
  rotation curves see the extra force. KiDS lensing (Brouwer+ 2021) and
  rotation curves agree on the same relation, which rules this out.
* **Coupled so that it also lenses:** this is relativistic MOND (AeST-like),
  whose extra pull follows the baryons, so it meets the Bullet Cluster filter.

**Verdict: no viable coupling, so it was not built.** What survives is the
a0–w relation as a curiosity. And since the scatter turned out not to be a
discriminator (previous section), the motivation for this kind of mechanism
is weaker than it looked this morning.
