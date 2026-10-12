# Roadmap: every test in UniverseBench

All 40 tests are implemented (38 plus the action check and whole-galaxy offsets 01b, Oct 2026). **Scored** tests use settled data a correct
theory could pass; **open questions** use disputed or mutually contradictory
data and are reported but never scored or used to fit constants. One rule
holds everywhere: a theory gets one set of global constants for all tests.

‡ marks tests with too few numbers to hide (fitted and scored in-sample);
the others have a secret public/hidden split. Guard rails (precision gravity,
black holes, clocks, inflation, neutrinos, curvature, standard sirens, particle
searches) are gates: standard physics passes them automatically, and they exist
to catch exotic theories. They count toward a model's lowest score but not its
median or its count of substantive tests passed. Every test also reports a
noise ceiling (what a perfect theory would score).

## Stage 1: galaxies and clusters

| # | Test | Data | Role | Scored? |
|---|---|---|---|---|
| 01 | Rotation curves: predict | SPARC, 165 galaxies (Lelli+ 2016) | Core MOND-vs-dark-matter test | Scored |
| 01b | Rotation curves: whole-galaxy offsets | SPARC quality 1, one offset per galaxy with its distance, inclination and M/L errors (Table 1) | Test 01 treats points as independent and barely discriminates; this scores the error that actually dominates | Scored |
| 02 | Rotation curves: extrapolate | SPARC; fit inner 60%, predict outer 40% | Punishes per-galaxy overfitting | Scored |
| 03 | Galaxy cluster masses | 13 Chandra clusters (Vikhlinin+ 2006) | MOND falls short by ~2x | Scored |
| 04 | Milky Way | Gaia/APOGEE curve (Eilers+ 2019), McMillan 2017 baryons | Pure out-of-sample galaxy | Scored ‡ |
| 05 | Dwarf spheroidals | 78 dwarfs from the Local Volume Database (Pace 2024): 66 satellites, 12 field; 9 upper limits | Dark-matter-dominated; MOND external field effect | Scored |
| 06 | Wide binary stars | Chae 2023 vs Banik+ 2024 (Gaia DR3) | Gravity at ultra-low acceleration | Open question |
| 07 | Weak lensing around galaxies | KiDS lensing (Mistele+ 2024), flat to ~1 Mpc | Gravity far outside galaxies | Scored ‡ |
| 08 | Strong-lensing ellipticals | 70 SLACS lenses (Auger+ 2009) | Mass in other galaxy types | Scored |
| 09 | Bullet Cluster | Clowe+ 2006 aperture masses | Mass separated from gas | Scored ‡ |
| 10 | Galaxies lacking dark matter | NGC 1052-DF2/DF4 | Distances disputed | Open question |

## Stage 2: precision gravity

| # | Test | Data | Role | Scored? |
|---|---|---|---|---|
| 11 | Solar system | Cassini gamma, MESSENGER beta, lunar laser ranging | Any change to gravity must vanish here | Scored ‡ |
| 12 | Binary pulsars | Hulse-Taylor, double pulsar | Gravitational radiation | Scored ‡ |
| 13 | Gravitational waves | GW170817 speed, GWTC-3 graviton mass | Ruled out many modified gravities | Scored ‡ |
| 14 | Laboratory gravity | MICROSCOPE, Eot-Wash | Equivalence principle, inverse-square law | Scored ‡ |
| 26 | Black holes | EHT shadows of Sgr A* and M87*, star S2 (GRAVITY) | Strongest-field tests | Scored ‡ |
| 27 | Clocks | Galileo satellites, Tokyo Skytree optical clocks | Gravitational time dilation | Scored ‡ |
| 37 | Healthy action | The model's declared action: no ghosts or gradient instabilities, GW speed = c (GW170817), declared backgrounds obey the Friedmann constraint derived from the action | Ties a theory to a consistent field theory; stops scale-by-scale stitching | Scored ‡ (guard rail) |

## Stage 3: cosmology

| # | Test | Data | Role | Scored? |
|---|---|---|---|---|
| 15 | Cosmic microwave background | Planck 2018 plik-lite TT/TE/EE spectra (616 decorrelated points) through the CLASS Boltzmann code; the old compressed summary (theta*, R, ...) is now a fit aid only | Strongest evidence for non-baryonic matter; no longer presupposes ΛCDM perturbations | Scored ‡ |
| 16 | Big Bang nucleosynthesis | PDG 2025: helium, deuterium | Independent count of ordinary matter | Scored ‡ |
| 16b | Lithium problem | PDG 2025 lithium-7 | Standard BBN overpredicts ~3x | Open question |
| 17 | Baryon acoustic oscillations | DESI DR2 (13 numbers + covariance) | Expansion history | Scored ‡ |
| 18 | Supernova distances | 1580 Pantheon+ supernovae | Expansion history, dark energy | Scored |
| 19 | Hubble constant (local) | SH0ES 2022, CCHP 2024 | Hubble tension | Open question |
| 20 | Clumpiness S8 (lensing) | DES Y3, KiDS-Legacy, HSC Y3 | S8 tension | Open question |
| 21 | Galaxy cluster counts | eROSITA eRASS1, SPT 2024 | Growth of structure | Scored ‡ |
| 22 | Lyman-alpha forest | PRIYA/eBOSS small-scale power; WDM mass > 5.7 keV | Small scales; tension reported | Open question |
| 23 | JWST early galaxies | Labbe+ 2023 via Boylan-Kolchin 2023 | Massive galaxies too early? | Open question |

## Stage 3b: extended cosmology

| # | Test | Data | Role | Scored? |
|---|---|---|---|---|
| 28 | Inflation signatures | BICEP/Keck r < 0.036; Planck f_NL | The earliest moments | Scored ‡ |
| 29 | Neutrinos | N_eff (Planck), KATRIN, oscillation minimum | Neutrino masses and species | Scored ‡ |
| 29b | Neutrino-mass squeeze | DESI DR2 + CMB, Planck mass bounds | Cosmology vs oscillations | Open question |
| 30 | Curvature and topology | Planck + BAO Omega_K; Planck topology limits | Shape of the universe | Scored ‡ |
| 31 | Growth rate | SDSS/eBOSS DR16 f sigma8 | Gravity on cosmic scales | Scored ‡ |
| 32 | Cosmic dipole | CatWISE quasar dipole (Secrest+ 2021) | Is the universe isotropic? | Open question |
| 33 | Dark flow | Planck bulk-flow limit | Large-scale motions | Open question |
| 34 | Standard sirens | GW170817 H0 and spacetime dimensions | Extra dimensions, GW propagation | Scored ‡ |
| 35 | Early supermassive black holes | J0313-1806 at z = 7.64 (Wang+ 2021) | How fast structure grew | Open question |

## Stage 4: particle guard rails

| # | Test | Data | Role | Scored? |
|---|---|---|---|---|
| 24 | Direct detection | LZ 2024 | Limits on dark matter scattering | Scored ‡ |
| 25 | Indirect detection | Fermi-LAT dwarfs | Limits on dark matter annihilation | Scored ‡ |
| 36 | Cosmic antimatter | BESS-Polar antihelium | No primordial antimatter domains | Scored ‡ |

## Possible extensions

- Fresh-data protocol: freeze submitted models and score them on new releases
  (Gaia DR4, Euclid, DESI DR3, LiteBIRD) they could not have seen.
- Full likelihoods (CMB spectra, lensing maps) for theories far from LCDM.
- More baselines: warm/fuzzy/self-interacting dark matter, superfluid dark
  matter, MOND + sterile neutrinos, evolving dark energy, early dark energy,
  emergent gravity, relativistic MOND (AeST).

Not included (outside the cosmology focus): Lorentz-invariance tests with
gamma-ray bursts, the flyby and Pioneer anomalies, the neutron lifetime puzzle
and particle-physics anomalies. Not testable with data: the theory of
everything, values of constants, the quantum gravity mechanism, the black hole
information paradox, holography, the problem of time, the cosmological
constant and coincidence *problems*, and the universe beyond what we observe.
