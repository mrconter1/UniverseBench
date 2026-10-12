# UniverseBench

*95% of the universe is dark and unexplained. Can AI explain it?*

UniverseBench scores theories of gravity, dark matter and dark energy against
real measurements, from the solar system to the cosmic microwave background,
with **one set of constants for everything** and hidden data where datasets
are large. It is built for AI agents: propose a theory, implement its
predictions, and see where it holds and where it breaks.

**38 tests: 27 scored, 11 open questions.**

| Stage | Scored tests | Open questions (reported, not scored) |
|---|---|---|
| 1 Galaxies and clusters | 01 rotation curves (predict), 02 rotation curves (extrapolate), 03 cluster masses, 04 Milky Way, 05 dwarf galaxies, 07 lensing around galaxies, 08 lensing ellipticals, 09 Bullet Cluster | 06 wide binary stars, 10 galaxies lacking dark matter |
| 2 Precision gravity | 11 solar system, 12 binary pulsars, 13 gravitational waves, 14 laboratory gravity, 26 black holes, 27 clocks | |
| 3 Cosmology | 15 microwave background, 16 nucleosynthesis, 17 BAO (DESI DR2), 18 supernovae (Pantheon+), 21 cluster counts, 28 inflation, 29 neutrinos, 30 curvature and topology, 31 growth rate, 34 standard sirens | 16b lithium, 19 Hubble tension, 20 S8, 22 Lyman-alpha, 23 JWST early galaxies, 29b neutrino-mass squeeze, 32 cosmic dipole, 33 dark flow, 35 early supermassive black holes |
| 4 Particle searches | 24 direct detection, 25 indirect detection, 36 cosmic antimatter | |

Every number comes from a published paper and was checked against its source
(abstracts or tables, two independent extractions or a built-in checksum where
transcribed); machine-readable data (DESI, SDSS, Pantheon+) were copied
unchanged from the collaborations' public repositories. Each file in `data/`
cites its source. Details per test: ROADMAP.md and the docstrings in
`ubench/tests/`.

## Scoring

* Per test: the % of data points predicted within 2 sigma, shown next to the
  test's **noise ceiling** (what a perfect theory would score given the error
  bars: about 95, with a 5th-95th percentile range for small tests). A test is
  passed when the score reaches that range.
* Per model: ranked by the **lowest** score over the scored tests (a theory is
  only as good as its worst failure), then by **substantive tests passed**, then by the
  **median** over substantive tests.
* **Guard rails** (◆: precision gravity, black holes, clocks, inflation,
  neutrinos, curvature, standard sirens, particle searches) are gates every
  viable theory must pass: they count
  toward the lowest score but not the median, so automatic passes don't
  inflate it.
* Open questions: disputed or mutually contradictory data. Even a correct
  theory could not pass them all, so they are reported, never scored, and never
  used to fit constants.
* Global constants are fitted once, jointly, on the public data of all scored
  tests. Per-galaxy knobs are only allowed in rotation-curve extrapolation.
* Large datasets (rotation curves, clusters, dwarfs, lensing ellipticals,
  supernovae) are split by a secret salt; scoring uses the hidden half. Small
  tables cannot be hidden and are scored in-sample (marked ‡).

## Baselines

`ubench/models/`: **ΛCDM** (cold dark matter + cosmological constant), **MOND**
(published a0, and a0 fitted to all tests), **standard gravity with visible
matter only**, and **pure Newtonian gravity** (no relativity). Current results:
LEADERBOARD.md.

## Layout

```
score.py                scorecard: python score.py [my_model.py] [--splits 5 --write]
make_agent_bundle.py    builds agent_bundle/ with public data only
AGENT_TASK.md           instructions and model interface for an AI agent
ROADMAP.md              every test: data, role, status
LEADERBOARD.md          latest baseline scores
ubench/core.py          hidden splits, Test base class, joint fitting, run()
ubench/physics.py       constants, disc gravity, halos, MOND interpolation
ubench/cosmology.py     cosmology engine (validated against Planck 2018 to <=0.03%)
ubench/tests/           one module per test family
ubench/models/          baseline theories
data/<test>/            data with source notes
```

## Workflow

    pip install numpy scipy
    export BENCH_SALT='your-private-secret'    # decides the hidden half
    python make_agent_bundle.py                # give agent_bundle/ to an agent
    # ...agent writes my_model.py, scoring itself on a practice split...
    python score.py agent_bundle/my_model.py   # final score on hidden data

A full five-split scorecard takes about 15 minutes on two CPU cores.

## Known simplifications

- Galaxies: no SPARC quality cuts; ΛCDM halos are plain NFW without baryonic
  feedback or tidal stripping; abundance matching uses Moster+ 2013.
- Dwarfs: single-epoch dispersions of the faintest systems may be inflated by
  binary stars; gas-rich field dwarfs also rotate a little, and the test uses
  their dispersions. Host fields assume flat host rotation curves.
- Clusters: stellar mass is 15% of the gas mass; the ΛCDM cluster model uses an
  empirical baryon fraction with mass and radius trends.
- The Bullet Cluster test is a crude mass-weighted peak position, not a lensing
  map. The Milky Way's stellar masses come from a fit that included a halo.
- Cosmology uses compressed summaries (CMB distance priors, fσ8, S8) derived
  assuming near-standard physics, an Eisenstein–Hu power spectrum (σ8 within
  0.5% of Planck) and fitting formulas for recombination and BBN. A radically
  different theory may need the full likelihoods.
- Several conversions of bounds to uncertainties are ours; each is noted in the
  relevant CSV.

## Related work

- Gravity-Bench (Koblischke et al., ICML 2025, arXiv:2501.18411): agents
  discover gravitational physics in simulated environments, some with
  deliberately non-standard laws. UniverseBench is the complement: real data
  and a real open problem.
- GAMBIT and Cobaya: global fits of physics models to many datasets, built for
  human researchers. UniverseBench adds hidden data, an agent-facing task, and
  one set of constants across all scales.
