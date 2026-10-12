# Re-score of the frozen superfluid submission

`my_model.py` MD5 e2ebbbcb78880983f7eef319a4563372 (unchanged from submission).
Benchmark state: commit fa591dd plus the cached-scorer fix (a0 unit fix
683aedd, and the fitter fix: scipy's bounded Powell search could end worse
than its start, which had pushed R_break to 359 kpc).

## Summary

| Run | Superfluid: lowest / passed / median | ΛCDM: lowest / passed / median | Rank of superfluid |
|---|---|---|---|
| Agent's hidden split (old benchmark) | 33 / 10 / 100 | 33 / 9 / 92 | 1 |
| Agent's 5 demo splits (old benchmark) | 47 / 10 / 100 | 37 / 9 / 90 | 1 |
| **5 demo splits, fixed benchmark** | **64 / 10 / 100** | 37 / 9 / 90 | **1** |
| **Fresh secret salt, fixed benchmark** | **33 / 10 / 100** | 37 / 9 / 88 | **2** |

Fitted breaking radius per split: 44, 59, 60, 61, 63, 71 kpc (was 359 kpc).
Dwarf test 05 rose from 47 to 65 (5 splits) and from 50 to 62 on demo-salt.

## What this shows

1. **A score gain without physics behind it.** With the fitter fixed, the
   condensate-breaking rule lands in its basin and test 05 rises by about
   18 points. The agent's own follow-up on modern data (Local Volume
   Database) found no orbital switch in the dwarfs and withdrew this
   mechanism. A constant tuned for a mechanism the data don't support
   still earned points, on outdated dispersions (5 of the 25 differ from
   modern values by more than 2 sigma).
2. **The ranking is not stable.** On the fresh salt, 1 of 3 hidden clusters
   is missed, so the cluster score is 33 and ΛCDM ranks first. Across the six
   splits, the superfluid model's cluster score ranges from 33 to 92. The
   lowest-score ranking is decided by hidden sets of 3 to 7 objects.
3. **Rotation curves** stay at MOND's level (61 to 64), as the agent reported.

Contamination: after this submission, the agent analysed all benchmark
dwarfs, hidden ones included. Test 05 is not a held-out check for any later
version of this model.

## Five demo splits (fixed benchmark)

| # | Test | Perfect theory | Two-phase superfluid dark matter | Cold dark matter (ΛCDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|---|---|
| 01 | Rotation curves: predict | 95 (95–96) | 64 | 57 | 65 | 57 | 29 | 29 |
| 02 | Rotation curves: extrapolate | 95 (94–97) | 80 | 37 | 79 | 66 | 19 | 19 |
| 03 | Galaxy clusters | 96 (82–100) | 76 | 96 | 11 | 42 | 0 | 0 |
| 04 | Milky Way | 96 (90–100) | 100 | 87 | 100 | 16 | 0 | 0 |
| 05 | Dwarf galaxies | 95 (83–100) | 65 | 47 | 37 | 38 | 14 | 14 |
| 07 | Lensing around galaxies ‡ | 95 (88–100) | 100 | 50 | 100 | 100 | 0 | 0 |
| 08 | Elliptical galaxies (lensing) | 96 (88–100) | 92 | 93 | 73 | 80 | 71 | 71 |
| 09 | Bullet Cluster § ‡ | 95 (50–100) | 100 | 100 | 0 | 0 | 0 | 0 |
| 11 | Solar system ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | 100 | 25 |
| 12 | Binary pulsars ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 0 |
| 13 | Gravitational waves ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | — |
| 14 | Laboratory gravity ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 100 |
| 26 | Black holes ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | 100 | — |
| 27 | Clocks and time dilation ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 0 |
| 15 | Microwave background ‡ | 96 (67–100) | 100 | 100 | — | — | 0 | — |
| 16 | Big Bang nucleosynthesis ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 0 | — |
| 17 | Baryon acoustic oscillations ‡ | 96 (85–100) | 85 | 85 | — | — | 8 | — |
| 18 | Supernova distances | 95 (94–97) | 99 | 99 | — | — | 98 | — |
| 21 | Galaxy cluster counts ‡ | 96 (75–100) | 100 | 75 | — | — | 0 | — |
| 28 | Inflation signatures ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 29 | Neutrino masses and species ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 30 | Curvature and topology ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 31 | Growth rate of structure ‡ | 96 (75–100) | 100 | 100 | — | — | 0 | — |
| 34 | Standard sirens ◆ ‡ | 95 (50–100) | 100 | 100 | — | — | 100 | — |
| 24 | Direct detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 25 | Indirect detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 36 | Cosmic antimatter ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| | **Lowest** | | **64** | **37** | **0** | **0** | **0** | **0** |
| | **Substantive tests passed** (of 14) | | **10** | **9** | **3** | **2** | **1** | **0** |
| | Median (excl. guard rails) | | 100 | 90 | 24 | 27 | 0 | 0 |
| | Highest | | 100 | 100 | 100 | 100 | 100 | 100 |
| | Knobs | | 14 global, 1 per galaxy | 10 global, 3 per galaxy | 3 global, 1 per galaxy | 4 global, 1 per galaxy | 6 global, 1 per galaxy | 2 global, 1 per galaxy |
| | **Rank** | | 1 | 2 | 3 | 4 | 5 | 6 |

## Fresh secret salt (fixed benchmark)

| # | Test | Perfect theory | Cold dark matter (ΛCDM) | Two-phase superfluid dark matter | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|---|---|
| 01 | Rotation curves: predict | 95 (95–96) | 51 | 61 | 62 | 52 | 32 | 32 |
| 02 | Rotation curves: extrapolate | 95 (94–97) | 37 | 82 | 81 | 71 | 15 | 15 |
| 03 | Galaxy clusters | 95 (83–100) | 100 | 33 | 33 | 33 | 0 | 0 |
| 04 | Milky Way | 96 (90–100) | 87 | 100 | 100 | 8 | 0 | 0 |
| 05 | Dwarf galaxies | 96 (85–100) | 38 | 77 | 54 | 69 | 23 | 23 |
| 07 | Lensing around galaxies ‡ | 95 (88–100) | 50 | 100 | 100 | 100 | 0 | 0 |
| 08 | Elliptical galaxies (lensing) | 96 (89–100) | 89 | 89 | 79 | 82 | 75 | 75 |
| 09 | Bullet Cluster § ‡ | 95 (50–100) | 100 | 100 | 0 | 0 | 0 | 0 |
| 11 | Solar system ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | 100 | 25 |
| 12 | Binary pulsars ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 0 |
| 13 | Gravitational waves ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | — |
| 14 | Laboratory gravity ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 100 |
| 26 | Black holes ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | 100 | — |
| 27 | Clocks and time dilation ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 0 |
| 15 | Microwave background ‡ | 96 (67–100) | 100 | 100 | — | — | 0 | — |
| 16 | Big Bang nucleosynthesis ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 0 | — |
| 17 | Baryon acoustic oscillations ‡ | 96 (85–100) | 85 | 85 | — | — | 8 | — |
| 18 | Supernova distances | 95 (94–97) | 99 | 99 | — | — | 98 | — |
| 21 | Galaxy cluster counts ‡ | 96 (75–100) | 75 | 100 | — | — | 0 | — |
| 28 | Inflation signatures ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 29 | Neutrino masses and species ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 30 | Curvature and topology ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 31 | Growth rate of structure ‡ | 96 (75–100) | 100 | 100 | — | — | 0 | — |
| 34 | Standard sirens ◆ ‡ | 95 (50–100) | 100 | 100 | — | — | 100 | — |
| 24 | Direct detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 25 | Indirect detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 36 | Cosmic antimatter ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| | **Lowest** | | **37** | **33** | **0** | **0** | **0** | **0** |
| | **Substantive tests passed** (of 14) | | **9** | **10** | **3** | **2** | **1** | **0** |
| | Median (excl. guard rails) | | 88 | 100 | 44 | 21 | 0 | 0 |
| | Highest | | 100 | 100 | 100 | 100 | 100 | 100 |
| | Knobs | | 10 global, 3 per galaxy | 14 global, 1 per galaxy | 3 global, 1 per galaxy | 4 global, 1 per galaxy | 6 global, 1 per galaxy | 2 global, 1 per galaxy |
| | **Rank** | | 1 | 2 | 3 | 4 | 5 | 6 |
