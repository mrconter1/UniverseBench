# Roadmap: from galaxies to the whole universe

The goal is one benchmark where a candidate theory of gravity / dark matter
is scored against every major observation, with a hidden half of the data
in each test. One rule holds across everything: a theory gets **one set of
global parameters for all tests**. It cannot use one value of a constant for
galaxies and another for clusters.

Each test is tagged:

- **Role**: what it discriminates. "Guard rail" tests are ones the
  standard picture passes easily and which rule out whole classes of ideas.
- **Status**: settled (clean, trusted data) or contested (live anomaly or
  disputed analysis; weighted lower and used to test anomaly triage).
- **Effort**: S (table of numbers), M (needs modelling), L (needs a physics
  code run per theory).

## Stage 1: galaxies and clusters

| # | Test | Data | Role | Status | Effort | State |
|---|---|---|---|---|---|---|
| 1 | Rotation curves, predict | SPARC, 165 galaxies | Core MOND-vs-dark-matter test | Settled | S | Done |
| 2 | Rotation curves, extrapolate outer parts | SPARC | Punishes per-galaxy overfitting | Settled | S | Done |
| 3 | Galaxy cluster masses from X-ray gas | Chandra relaxed clusters | Guard rail: MOND is known to fall short by about 2x | Settled | S | Done |
| 4 | Milky Way rotation curve | Gaia (Eilers+ 2019 curve) | Our own galaxy, very precise | Settled | S | Done |
| 5 | Dwarf spheroidal galaxies | Local Group velocity dispersions | Dark-matter-dominated; tests MOND's external field effect | Mostly settled | M | **Next** |
| 6 | Wide binary stars | Gaia DR3 | Gravity at ultra-low acceleration | **Contested** | M | Planned |
| 7 | Weak lensing around galaxies | KiDS lensing acceleration relation | Rotation-curve test extended 10x further out | Settled-ish | M | Planned |
| 8 | Strong lensing plus dynamics of elliptical galaxies | SLACS | Mass in galaxy centers, other galaxy types | Settled | M | Planned |
| 9 | Bullet Cluster and other merging clusters | Lensing vs X-ray maps | Guard rail: mass separated from gas | Settled | M | Planned |
| 10 | Galaxies apparently lacking dark matter | NGC 1052-DF2/DF4 | Hard case for both camps | **Contested** | S | Planned |

## Stage 2: gravity guard rails

| # | Test | Data | Role | Status | Effort | State |
|---|---|---|---|---|---|---|
| 11 | Solar system | Planetary ephemerides, Cassini, lunar laser ranging | Any change to gravity must vanish here | Settled | S | Planned |
| 12 | Binary pulsars | Hulse-Taylor, double pulsar | Strong-field gravity, gravitational radiation | Settled | S | Planned |
| 13 | Gravitational waves | GW170817 speed, LIGO/Virgo waveform tests | Killed many modified-gravity theories | Settled | S | Planned |
| 14 | Laboratory gravity | Torsion balances, equivalence principle | Short-range and composition tests | Settled | S | Planned |

## Stage 3: cosmology (needs a physics code per theory)

| # | Test | Data | Role | Status | Effort | State |
|---|---|---|---|---|---|---|
| 15 | Cosmic microwave background | Planck, ACT, SPT | Strongest evidence for non-baryonic matter | Settled | L | Planned |
| 16 | Big Bang nucleosynthesis | Primordial He, D, Li | Independent count of ordinary matter; lithium anomaly | Settled + contested Li | S | Planned |
| 17 | Baryon acoustic oscillations | DESI | Expansion history; hints of evolving dark energy | Partly contested | M | Planned |
| 18 | Supernova distances | Pantheon+, DES-SN | Expansion history, dark energy | Settled | M | Planned |
| 19 | Local expansion rate | SH0ES, TRGB | Hubble tension | **Contested** | S | Planned |
| 20 | Weak lensing clumpiness (S8) | DES, KiDS | Growth of structure; mild tension | **Contested** | L | Planned |
| 21 | Galaxy cluster counts | SPT, eROSITA | Growth of structure | Settled-ish | L | Planned |
| 22 | Lyman-alpha forest | Quasar spectra | Small-scale structure; warm/fuzzy dark matter | Settled-ish | L | Planned |
| 23 | Early massive galaxies | JWST | Possible tension, messy astrophysics | **Contested** | L | Planned |

## Stage 4: optional particle constraints

| # | Test | Data | Role | Status | Effort | State |
|---|---|---|---|---|---|---|
| 24 | Direct detection limits | LZ, XENONnT | Rules out particle dark matter models | Settled | S | Optional |
| 25 | Collider and indirect searches | LHC, Fermi gamma rays | Same | Settled | M | Optional |

## Benchmark infrastructure (alongside the tests)

- Joint fitting of global parameters across all tests (enforces the one-constant rule). Done.
- Scorecard of every model on every test, '—' for no prediction. Done.
- Weighted total score: settled tests full weight, contested tests reduced.
- Noise ceiling per test (what a perfect theory would score).
- Planned baselines: warm/fuzzy/self-interacting dark matter, superfluid dark
  matter, MOND + sterile neutrinos, evolving dark energy, early dark energy,
  emergent gravity, relativistic MOND (AeST).
- Per-test report of which known anomalies a model resolves.
- Agent bundle per stage, with hidden data never leaving the scorer.
