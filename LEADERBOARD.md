# UniverseBench leaderboard

## Benchmark (settled tests)

| # | Test | Perfect theory | Cold dark matter (ΛCDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|---|
| 01 | Rotation curves: predict | 95 (95–96) | 57 | 65 | 57 | 29 | 29 |
| 02 | Rotation curves: extrapolate | 95 (94–97) | 37 | 79 | 66 | 19 | 19 |
| 03 | Galaxy clusters | 96 (82–100) | 96 | 11 | 42 | 0 | 0 |
| 04 | Milky Way | 96 (90–100) | 87 | 100 | 16 | 0 | 0 |
| 05 | Dwarf galaxies | 95 (83–100) | 47 | 37 | 38 | 14 | 14 |
| 07 | Lensing around galaxies ‡ | 95 (88–100) | 50 | 100 | 100 | 0 | 0 |
| 08 | Elliptical galaxies (lensing) | 96 (88–100) | 93 | 73 | 80 | 71 | 71 |
| 09 | Bullet Cluster § ‡ | 95 (50–100) | 100 | 0 | 0 | 0 | 0 |
| 11 | Solar system ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | 25 |
| 12 | Binary pulsars ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 0 |
| 13 | Gravitational waves ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | — |
| 14 | Laboratory gravity ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 |
| 26 | Black holes ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | — |
| 27 | Clocks and time dilation ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 0 |
| 15 | Microwave background ‡ | 96 (67–100) | 100 | — | — | 0 | — |
| 16 | Big Bang nucleosynthesis ‡ | 95 (50–100) | 100 | 100 | 100 | 0 | — |
| 17 | Baryon acoustic oscillations ‡ | 96 (85–100) | 85 | — | — | 8 | — |
| 18 | Supernova distances | 95 (94–97) | 99 | — | — | 98 | — |
| 21 | Galaxy cluster counts ‡ | 96 (75–100) | 75 | — | — | 0 | — |
| 28 | Inflation signatures ◆ ‡ | 98 | 100 | — | — | 100 | — |
| 29 | Neutrino masses and species ◆ ‡ | 98 | 100 | — | — | 100 | — |
| 30 | Curvature and topology ◆ ‡ | 98 | 100 | — | — | 100 | — |
| 31 | Growth rate of structure ‡ | 96 (75–100) | 100 | — | — | 0 | — |
| 34 | Standard sirens ◆ ‡ | 95 (50–100) | 100 | — | — | 100 | — |
| 24 | Direct detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 |
| 25 | Indirect detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 |
| 36 | Cosmic antimatter ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 |
| | **Lowest** | | **37** | **0** | **0** | **0** | **0** |
| | **Substantive tests passed** (of 14) | | **9** | **3** | **2** | **1** | **0** |
| | Median (excl. guard rails) | | 90 | 24 | 27 | 0 | 0 |
| | Highest | | 100 | 100 | 100 | 100 | 100 |
| | Knobs | | 10 global, 3 per galaxy | 3 global, 1 per galaxy | 4 global, 1 per galaxy | 6 global, 1 per galaxy | 2 global, 1 per galaxy |
| | **Rank** | | 1 | 2 | 3 | 4 | 5 |

Scores: % of data points predicted within 2 sigma. 'Perfect theory' is the noise ceiling: the expected score of the true laws of nature given each test's error bars (with the 5th-95th percentile range for small tests). A test is passed when the score reaches the 5th percentile of that range. The benchmark uses settled tests only; models are ranked by their lowest score (a theory is only as good as its worst failure), then substantive tests passed, then median. ◆ = guard rail: a gate every viable theory must pass, counted in the lowest score but not the median. § = simplified test (09: mass-weighted peak position from aperture masses, not a fit to the lensing map; a faithful model should compute where its lensing mass sits). '—' = no prediction (counts as 0). ‡ = too few numbers to hide: fitted and scored in-sample. Averaged over 5 random public/hidden splits.

## Open questions

| # | Open question (not scored) | Cold dark matter (ΛCDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|
| 06 | Wide binary stars | 50 | 50 | 10 | 50 | 50 |
| 10 | Galaxies lacking dark matter | 0 | 50 | 0 | 100 | 100 |
| 16b | Lithium problem | 0 | 0 | 0 | 0 | — |
| 19 | Hubble constant (local) | 50 | — | — | 0 | — |
| 20 | Clumpiness S8 (lensing) | 67 | — | — | 0 | — |
| 22 | Lyman-alpha forest | 67 | — | — | 33 | — |
| 23 | JWST early galaxies | 100 | — | — | 0 | — |
| 29b | Neutrino mass squeeze | 100 | — | — | 100 | — |
| 32 | Cosmic dipole (isotropy) | 0 | — | — | 0 | — |
| 33 | Dark flow (bulk motions) | 100 | — | — | 100 | — |
| 35 | Early supermassive black holes | 0 | — | — | 100 | — |

Open questions: tests whose data are disputed or contradict each other, so even a correct theory could not pass them all. Reported, never scored, never used to fit constants.
