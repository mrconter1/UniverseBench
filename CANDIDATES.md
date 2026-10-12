# Candidate theories to test

Written 12 Oct 2026 after the first solver run. Percentages are **gut
feelings** (Claude's), not results. They answer "what is actually going on?"
and add up to 100 across the main list (updated the evening of 12 Oct, see below).

## Main list

| # | Explanation | Gut | Why |
|---|---|---|---|
| 1 | Cold dark matter; the galaxy regularity comes from galaxy formation (feedback, gas) | 40% | Simplest story that fits cosmology, clusters and dwarfs; nobody has shown it gives such a tight galaxy relation |
| 2 | A dark sector with its own forces; self-interaction is the first symptom (cSIDM) | 15% | Matches MOND on galaxies and passes everything else with no per-galaxy freedom (submissions/solver_claude) |
| 3 | A different particle: ultra-light (wave) dark matter, warm, or a mixture | 10% | Untested here; wave dark matter naturally gives dense dwarf cores |
| 4 | Hybrid: real dark matter plus a medium or force acting in some conditions (superfluid, "two horizons", MOND + hot DM) | 8% | Explains why galaxies "know" about the cosmic horizon (a0 ≈ cH0); dwarfs showed no sign of the superfluid switch |
| 5 | Modified gravity alone (MOND, relativistic MOND) | 5.5% | Excellent on galaxies; clusters, the Bullet Cluster, the CMB and dwarf halos count against it |
| 6 | Emergent or horizon-based gravity (Verlinde-like) | 3% | Explains a0 ≈ cH0; behaves like MOND where it fails |
| 7 | Gravity with memory (nonlocal gravity, Mashhoon) | 0.5% | Galaxy relation for free; likely killed by the CMB |
| 8 | Something on no one's list | 18% | Historically, the answer is often not among the candidates argued about |

Cross-cutting (before test A; now ~15%, see table): given particle dark matter, ~40% that the galaxy relation is best
understood as a **relaxation attractor**, with a0 set by cosmic age. Dark
energy: ~60% a cosmological constant, 35% an evolving field, 5% linked to
dark matter.

## Assumption-breaking list

Each row denies one assumption behind the standard picture. Its percentages
are **not additive** with the main list above (which sums to 100): most rows
are sub-cases of a main-list entry (B lies within #2, C and O within #4,
D within #7, E and F within #5).

| # | Assumption denied | Idea | Status | Gut | Cheap to test here? |
|---|---|---|---|---|---|
| A | Galaxies are equilibrium snapshots | Galaxy relation as relaxation attractor; a0 set by age | **Tested with SPARC Table 1: not supported** (the signal was data quality; best data give 0.06 dex intrinsic scatter, nearly independent of evolution), see analysis/RESULTS.md | ~15% given particle DM | Done |
| B | Dark matter and dark energy are unrelated | Coupled quintessence / interacting dark energy | **Tested: mildly favoured** (ξ ≈ −0.07; BAO 85→92, cluster counts 75→100, median 90→95; not significant) | 5% | Done (approximate growth) |
| C | Gravity has no medium | Gravitational polarization / dipolar dark matter (Blanchet) | **Proxy tested: lowest 0** (MOND + CDM cosmology). The full theory is untested: where does the medium go in a collision? | 1.5% | Proxy done; collision check open |
| D | Gravity acts on the present only | Nonlocal gravity with memory | **Ruled out** by the Bullet Cluster for any memory time (long: peak lags galaxies; short: follows gas), see analysis/RESULTS.md | <0.5% | Done |
| E | G is the same at all scales | Running G from quantum-gravity renormalization | Exists; **fails the Bullet Cluster by construction** (extra gravity follows present baryons) | 0.5% | Closed by C |
| F | Inertia is fixed | Inertia set by Unruh vs horizon temperature | Exists; **same Bullet Cluster failure** | 0.5% | Closed by C |
| G | Three space dimensions at all scales | Extra dimensions, braneworld or mirror sectors | Established, mostly constrained | 1.5% | Partly |
| H | The universe is smooth on average | Backreaction or local void (dark energy, Hubble tension) | Established, contested | 3% (dark-energy side) | Partly |
| I | Galaxies' clocks are cosmic time | SIDM relaxation sets a0 through halo age; a0 rises with redshift | **Tested: fails** (clock lowers a0; cSIDM predicts 1.1–1.2× at z=1), see analysis/RESULTS.md | — | Done |
| J | The cosmic horizon doesn't act locally | Transition radius = √(r_Schwarzschild × R_Hubble) | Partly in Verlinde-type work | 1% | Yes |
| K | Mass is positive | Negative-mass fluid (Farnes 2018) | Tested and widely criticized | <0.5% | Partly |
| L | Spacetime is smooth | Topological defects as dark matter | Established, constrained | <0.5% | Hard |
| M | Dark matter is particles | Primordial black holes (as a fraction) | Established, constrained | minor | Partly |
| O | Dark matter is cold | MOND plus a hot dark component (~11 eV sterile neutrinos, Angus 2009): real mass in clusters, MOND in galaxies | Survived the compressed tests (lowest 41), but **the full Planck spectra exclude it** (χ² ≈ 2580 vs ≈ 590 for ΛCDM, 74% vs 96% within 2σ); a non-thermal (colder) sterile variant remains untested | 0.5% | CMB done; non-thermal variant open |
| N | There is one universe | Gravity leaking from parallel branches | No worked theory found | ~0 | No |

About 22 candidates in total, of which about 10 can be tested in an afternoon.
The radical ones need a full CMB calculation, which the benchmark lacks: it
uses compressed summaries that assume near-standard physics. Closing that gap
would let them be judged fairly.

## Second set (user, 12 Oct 2026)

| # | Idea | Status |
|---|---|---|
| 6 | Dark matter heated by baryons (Famaey, Khoury & Penco 2017/2019) | **examined**: a0 is an input and the scatter is not predicted; benchmark proxy = MOND in equilibrated discs, CDM elsewhere (best case) |
| 7 | Fermionic (keV) dark matter with dense quantum cores replacing central black holes (Argüelles, Ruffini et al.) | **tested: excluded** except 345–383 keV at the edge of collapse (GRAVITY flares need the mass inside 8 r_g; 56 keV gives 1.5e-7 of it); no photon sphere for any mass |
| 8 | Heavier electron at recombination (Hubble tension) | **tested**: H0 → 70.6 with the CMB intact, but summed χ² worse than ΛCDM (+14); the apparent preference came from equal test weighting |
| 9 | Cosmologically coupled black holes as dark energy | **tested: disfavoured** (phantom history, Δχ² ≈ +16 at equal constants, mostly supernovae; GW bounds k < 2; 46% of baryons consumed) |
| 10 | One rolling field drives dark energy and DM self-interaction (user's own) | **tested**: a0 responds weakly (d ln a0/d ln σ ≈ 0.2–0.3); not a sharp test in this model |

## Evening update (12 Oct 2026)

Inverse design showed no dark matter with one universal equation of state can
produce the galaxy relation; a per-galaxy link to the baryonic mass is
needed. But the tightness itself proved softer than assumed: the observed
intrinsic scatter is 0.03–0.08 dex depending on the error model, and full
hydrodynamic cold-dark-matter simulations give 0.04–0.09 (analysis/RESULTS.md).
Gut update: #1 (CDM + galaxy formation) 30 → 40%, #2 (cSIDM) 20 → 15%,
#4 (hybrids) 11.5 → 8%, #5 (modified gravity) 6 → 5.5%, #8 (unknown) 19 → 18%. The
dark-energy-tied dark force was stopped on paper (no viable coupling).

## Test queue (cheapest and most informative first)

1. **A:** does the scatter around the galaxy relation shrink for older,
   gas-poor, quieter SPARC galaxies? (15 min)
2. **I:** does the cSIDM clock reproduce MUSE-DARK's ~2.5× rise of a0 by
   redshift 1? (15 min)
3. **C, E, F** as quick benchmark models: does any MOND-family variant rescue
   clusters? (about 1 h)
4. **B** in the cosmology engine: the most plausible dark matter–dark energy
   link. (about 1 h)

Also queued from the solver run: an independent check by a fresh agent that has
not seen the data, deriving a0 by abundance matching instead of fitting it
(**done**: a0 emerges within a factor 1.6 from cosmology, see analysis/RESULTS.md),
scoring over 20 or more splits, and group-scale lensing (KiDS-GAMA, Viola et
al. 2015) as a new open question.

## Sources

- IAC, [Can a negative-mass cosmology explain dark matter and dark energy?](https://iac.es/en/science-and-technology/publications/can-negative-mass-cosmology-explain-dark-matter-and-dark-energy)
- Blanchet et al., [Dipolar dark matter and cosmology](https://ar5iv.arxiv.org/html/1312.6991)
- Famaey & McGaugh, [MOND phenomenology and the concordance cosmological scenario](https://arxiv.org/pdf/0807.1200)
- Smolin, [MOND as a regime of quantum gravity](https://arxiv.org/pdf/1704.00780)
- Rodrigues et al., [Galaxy rotation curves from GR with renormalization-group corrections](https://ar5iv.arxiv.org/html/0911.4967)
- Rahvar & Mashhoon, [Observational tests of nonlocal gravity](https://arxiv.org/pdf/1401.4819)
- Li et al., [The intrinsic scatter of the radial acceleration relation](https://ar5iv.labs.arxiv.org/html/1908.06105)
- [The radial acceleration relation at the EDGE of galaxy formation](https://arxiv.org/pdf/2510.06905)
- [Braneworld effects as an alternative to dark matter](https://ar5iv.labs.arxiv.org/html/astro-ph/0512494)
