# UniverseBench leaderboard

## Benchmark (settled tests)

| # | Test | Cold dark matter (ΛCDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|
| 01 | Rotation curves: predict | 57 | 64 | 58 | 29 | 29 |
| 02 | Rotation curves: extrapolate | 37 | 79 | 66 | 19 | 19 |
| 03 | Galaxy clusters | 96 | 12 | 42 | 0 | 0 |
| 04 | Milky Way | 87 | 100 | 27 | 0 | 0 |
| 05 | Dwarf galaxies | 47 | 37 | 38 | 14 | 14 |
| 07 | Lensing around galaxies ‡ | 50 | 100 | 100 | 0 | 0 |
| 08 | Elliptical galaxies (lensing) | 93 | 76 | 80 | 71 | 71 |
| 09 | Bullet Cluster ‡ | 100 | 0 | 0 | 0 | 0 |
| 11 | Solar system ‡ | 100 | 100 | 100 | 100 | 25 |
| 12 | Binary pulsars ‡ | 100 | 100 | 100 | 100 | 0 |
| 13 | Gravitational waves ‡ | 100 | 100 | 100 | 100 | — |
| 14 | Laboratory gravity ‡ | 100 | 100 | 100 | 100 | 100 |
| 15 | Microwave background ‡ | 100 | — | — | 40 | — |
| 16 | Big Bang nucleosynthesis ‡ | 100 | 100 | 100 | 0 | — |
| 17 | Baryon acoustic oscillations ‡ | 85 | — | — | 8 | — |
| 18 | Supernova distances | 99 | — | — | 98 | — |
| 21 | Galaxy cluster counts ‡ | 75 | — | — | 0 | — |
| | **Lowest** | **37** | **0** | **0** | **0** | **0** |
| | **Median** | 96 | 76 | 58 | 19 | 0 |
| | **Highest** | 100 | 100 | 100 | 100 | 100 |
| | Knobs | 10 global, 3 per galaxy | 3 global, 1 per galaxy | 4 global, 1 per galaxy | 6 global, 1 per galaxy | 2 global, 1 per galaxy |
| | **Rank** | 1 | 2 | 3 | 4 | 5 |

Scores: % of data points predicted within 2 sigma (a perfect theory scores about 95). The benchmark score uses settled tests only: models are ranked by their lowest settled-test score (a theory is only as good as its worst failure), then median. '—' = no prediction (counts as 0). ‡ = too few numbers to hide: fitted and scored in-sample. Averaged over 5 random public/hidden splits.

## Open questions

| # | Open question (not scored) | Cold dark matter (ΛCDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|
| 06 | Wide binary stars | 50 | 50 | 50 | 50 | 50 |
| 10 | Galaxies lacking dark matter | 0 | 50 | 0 | 100 | 100 |
| 16b | Lithium problem | 0 | 0 | 0 | 0 | — |
| 19 | Hubble constant (local) | 50 | — | — | 0 | — |
| 20 | Clumpiness S8 (lensing) | 67 | — | — | 0 | — |
| 22 | Lyman-alpha forest | 67 | — | — | 33 | — |
| 23 | JWST early galaxies | 100 | — | — | 0 | — |

Open questions: tests whose data are disputed or contradict each other, so even a correct theory could not pass them all. Reported, never scored, never used to fit constants.
