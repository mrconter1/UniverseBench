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
ubench/cosmology.py     fast cosmology engine (validated against Planck 2018 to <=0.03%)
ubench/class_engine.py  cached wrapper around CLASS (full CMB spectra for test 15)
ubench/cmb_lite.py      Planck 2018 plik-lite likelihood (own implementation)
ubench/tests/           one module per test family
ubench/models/          baseline theories
data/<test>/            data with source notes
submissions/            agent submissions, re-scores, solver run
CANDIDATES.md           candidate theories still to test, with gut-feeling odds
```

## Workflow

    pip install numpy scipy classy numba     # classy = CLASS Boltzmann code (test 15); numba optional, speeds up some models
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
- The CMB is scored on Planck's full binned spectra through CLASS (since
  12 Oct 2026). Other cosmology tests still use compressed summaries (fσ8, S8)
  derived assuming near-standard physics, an Eisenstein–Hu power spectrum (σ8 within
  0.5% of Planck) and fitting formulas for recombination and BBN. A radically
  different theory may need the full likelihoods.
- Several conversions of bounds to uncertainties are ours; each is noted in the
  relevant CSV.

## Project timeline

UniverseBench was built and first tried in about 18 hours (11–12 Oct 2026),
by a person working with Claude, plus a separate Claude session acting as
the first agent. Commit hashes refer to this repository.

**11 Oct, evening: the benchmark.** The idea grew out of a discussion of
Terence Tao's "Math 2.0" talk: could real data play the role a proof checker
plays in mathematics, so that AI agents could attack open physics problems
and be scored? Starting from galaxy rotation curves and cluster masses
(bc8503c), tests were added one family at a time: the Milky Way, dwarf
galaxies, wide binaries, galaxy lensing, elliptical lenses, the Bullet
Cluster, precision gravity, a cosmology engine validated against Planck 2018,
supernovae, BAO, black holes, clocks, extended cosmology and particle
searches. That made 38 tests (c963e24). Contested data became unscored "open
questions" (3d26276). Noise ceilings and guard-rail gates fixed a ranking
quirk where visible matter alone outranked MOND (51f525a).

**12 Oct, night: the pilot.** An agent in another session, given only the
public half of the data, built a two-phase superfluid dark matter model and
ranked first on the old benchmark. It also found two benchmark bugs: a 5%
error in the a0 unit conversion, and a fitter that could not handle
threshold-like constants (683aedd).

**12 Oct, morning: verification.** Re-scoring the frozen pilot model showed
that scipy's bounded Powell search could end worse than it started; fixed by
keeping the best point ever evaluated plus coordinate sweeps (fa591dd). With
the fix, the pilot model scored lowest 64 over five splits but 33 on a fresh
secret split. The ranking turned on a single hidden cluster, so small hidden
sets make "lowest score" noisy (25e0d74, submissions/superfluid_dm/RESCORE.md).
The pilot agent's follow-up on modern dwarf data withdrew its dwarf mechanism.
It also showed that 5 of the 25 dwarf dispersions in test 05 were outdated,
so test 05 was rebuilt on the Local Volume Database: 78 dwarfs, upper limits
kept, tidal systems dropped (47010a3).

**12 Oct, midday: Claude as solver.** Claude pre-registered a secret hidden
split (hash committed first, 3463747), developed only on its public half, and
built a collapsing self-interacting dark matter model (cSIDM). One
speed-dependent cross-section gives three regimes: cored galaxy halos that
track the baryons, collapsed (cuspy) dwarf halos, and collisionless clusters.
Scored once on the pre-registered split it ranked first: lowest 56, 10 of 14
substantive tests passed, median 94. ΛCDM scored 40, 9 and 92 on the same
split; MOND 0. Caveats: the author had seen all the data, there are 17 fitted
constants in total, and rotation curves only tie MOND (12d5f6b, 126ffcf,
submissions/solver_claude/).

**12 Oct, afternoon: what next.** A systematic list of candidate theories,
including ones outside established ideas, with gut-feeling odds and which can
be tested cheaply: CANDIDATES.md.

| Submission | Idea | Hidden-split result | Folder |
|---|---|---|---|
| Pilot agent | Two-phase superfluid dark matter | lowest 33 on a fresh split (rank 2) | submissions/superfluid_dm |
| Claude (solver) | Collapsing self-interacting dark matter | lowest 56, rank 1 (pre-registered split) | submissions/solver_claude |

## Data credits

All data are public and are redistributed here with their sources cited in
each `data/<test>/` file. Main sources: SPARC (Lelli, McGaugh & Schombert
2016), the Local Volume Database (Pace 2024, CC0), Gaia/APOGEE Milky Way
curve (Eilers et al. 2019), Chandra clusters (Vikhlinin et al. 2006), SLACS
(Auger et al. 2009), KiDS lensing (Mistele et al. 2024), Planck 2018, DESI
DR2, Pantheon+, SDSS/eBOSS, eROSITA and SPT, and the precision-gravity, black
hole, clock and particle-search papers listed per test. Please cite the
original papers when using these data.

## Related work

- Gravity-Bench (Koblischke et al., ICML 2025, arXiv:2501.18411): agents
  discover gravitational physics in simulated environments, some with
  deliberately non-standard laws. UniverseBench is the complement: real data
  and a real open problem.
- GAMBIT and Cobaya: global fits of physics models to many datasets, built for
  human researchers. UniverseBench adds hidden data, an agent-facing task, and
  one set of constants across all scales.
