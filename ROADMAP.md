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
| 5 | Dwarf spheroidal galaxies | Local Group velocity dispersions | Dark-matter-dominated; tests MOND's external field effect | Mostly settled | M | Done |
| 6 | Wide binary stars | Gaia DR3 | Gravity at ultra-low acceleration | **Contested** | M | Done |
| 7 | Weak lensing around galaxies | KiDS lensing acceleration relation | Rotation-curve test extended 10x further out | Settled-ish | M | Done |
| 8 | Strong lensing plus dynamics of elliptical galaxies | SLACS | Mass in galaxy centers, other galaxy types | Settled | M | Done |
| 9 | Bullet Cluster and other merging clusters | Lensing vs X-ray maps | Guard rail: mass separated from gas | Settled | M | Done |
| 10 | Galaxies apparently lacking dark matter | NGC 1052-DF2/DF4 | Hard case for both camps | **Contested** | S | Done |

## Stage 2: gravity guard rails

| # | Test | Data | Role | Status | Effort | State |
|---|---|---|---|---|---|---|
| 11 | Solar system | Planetary ephemerides, Cassini, lunar laser ranging | Any change to gravity must vanish here | Settled | S | Done |
| 12 | Binary pulsars | Hulse-Taylor, double pulsar | Strong-field gravity, gravitational radiation | Settled | S | Done |
| 13 | Gravitational waves | GW170817 speed, LIGO/Virgo waveform tests | Killed many modified-gravity theories | Settled | S | Done |
| 14 | Laboratory gravity | Torsion balances, equivalence principle | Short-range and composition tests | Settled | S | Done |
| 26 | Black holes | EHT shadows of M87* and Sgr A*; star S2 around Sgr A* (GRAVITY); merger ringdown | Strongest-field tests of gravity | Settled | S | Planned |
| 27 | Clocks and time dilation | Gravity Probe A, optical clocks at different heights | Gravitational redshift of time | Settled | S | Planned |

## Stage 3: cosmology (tables of numbers; predictions from ubench/cosmology.py)

| # | Test | Data | Role | Status | Effort | State |
|---|---|---|---|---|---|---|
| 15 | Cosmic microwave background | Planck, ACT, SPT | Strongest evidence for non-baryonic matter | Settled | L | Done |
| 16 | Big Bang nucleosynthesis | Primordial He, D, Li | Independent count of ordinary matter; lithium anomaly | Settled + contested Li | S | Done |
| 17 | Baryon acoustic oscillations | DESI | Expansion history; hints of evolving dark energy | Partly contested | M | Done |
| 18 | Supernova distances | Pantheon+, DES-SN | Expansion history, dark energy | Settled | M | Done |
| 19 | Local expansion rate | SH0ES, TRGB | Hubble tension | **Contested** | S | Done |
| 20 | Weak lensing clumpiness (S8) | DES, KiDS | Growth of structure; mild tension | **Contested** | L | Done |
| 21 | Galaxy cluster counts | SPT, eROSITA | Growth of structure | Settled-ish | L | Done |
| 22 | Lyman-alpha forest | Quasar spectra | Small-scale structure; warm/fuzzy dark matter | Settled-ish | L | Done |
| 23 | Early massive galaxies | JWST | Possible tension, messy astrophysics | **Contested** | L | **In progress** |

## Stage 3b: extended cosmology (queued from the list of unsolved problems in physics)

| # | Test | Data | Role | Status | Effort | State |
|---|---|---|---|---|---|---|
| 28 | Inflation signatures | Primordial gravitational waves r < 0.036 (BICEP/Keck), non-Gaussianity f_NL and running of n_s (Planck) | Tests the earliest moments; constrains alternatives to inflation | Settled | S | Queued |
| 29 | Neutrino masses and extra species | Sum of masses from DESI + CMB vs oscillation minimum (~0.06 eV); KATRIN direct limit; N_eff from CMB and BBN | Cosmology vs particle physics; sterile neutrinos | Partly contested | S | Queued |
| 30 | Curvature and topology | Existing CMB/BAO/supernova distances with curvature allowed in the engine; Planck limits on a finite (multiply connected) universe | Shape of the universe | Settled | M | Queued |
| 31 | Growth rate of structure | f sigma8 from redshift-space distortions (BOSS/eBOSS, DESI); models may vary dark energy w0, wa | Gravity on cosmic scales; dynamical dark energy | Settled | M | Queued |
| 32 | Cosmological principle | Quasar and radio-galaxy number-count dipole vs CMB kinematic dipole; large-scale homogeneity scale | Is the universe isotropic and homogeneous? | **Contested** | S | Queued |
| 33 | Dark flow | Bulk flows from the kinetic Sunyaev-Zel'dovich effect (Planck) and peculiar-velocity surveys | Large-scale motions beyond LCDM expectations | **Contested** | S | Queued |
| 34 | Gravitational-wave standard sirens | GW170817 distance vs electromagnetic distance (leakage into extra dimensions; GW friction); siren H0 | Extra dimensions; gravity propagation over cosmic distances | Settled | S | Queued |
| 35 | Early supermassive black holes | Billion-solar-mass quasars at z > 7; JWST broad-line AGN at z > 5 | How fast cosmic structure grew | **Contested** | M | Queued |
| 36 | Cosmic antimatter | Absence of antimatter domains (gamma-ray background), AMS-02 antinuclei limits | Guard rail for exotic early-universe models | Settled | S | Queued (optional) |

Not queued (outside the cosmology focus): Lorentz-invariance tests with gamma-ray bursts, the flyby and Pioneer anomalies, the neutron lifetime puzzle (noted as a BBN systematic), and particle-physics anomalies.
Not testable with data (left out on purpose): theory of everything, values of constants, quantum gravity mechanism, black hole information paradox, holography, problem of time, the cosmological constant and coincidence *problems* (why, not what), and the size of the universe beyond the observable part.

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
