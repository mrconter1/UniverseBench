# Scores: collapsing SIDM

Frozen model: `my_model.py`, MD5 df615fef6d3e06c028508e8d21281e1a (commit 12d5f6b, made before
the hidden run).

## Summary

| Run | cSIDM (lowest / passed / median) | ΛCDM | MOND (published) | cSIDM rank |
|---|---|---|---|---|
| Practice split (public data only) | 61 / 10 / 91 | 32 / 7 / 84 | 0 / 3 / 21 | 1 |
| **Pre-registered hidden split (one shot)** | **56 / 10 / 94** | 40 / 9 / 92 | 0 / 3 / 29 | **1** |

Galaxy-scale tests, hidden split: rotation curves 56 (MOND 57, ΛCDM 52);
extrapolation 76 (MOND 79, ΛCDM 40); dwarfs 58 (ΛCDM 47, MOND 42); galaxy
lensing 7 of 8 (MOND 8 of 8, ΛCDM 4 of 8); Milky Way 100; clusters 100;
Bullet Cluster 100; ellipticals 88 (ΛCDM 98).

Dwarfs dropped from 75% on the public half to 58% on the hidden half: the
frozen galaxy constants were fitted on the public dwarfs, and dwarf
diversity is not modelled.

The weakest test (lowest score) is rotation-curve prediction, where the
model ties MOND rather than beating it.

## Salt reveal

The pre-registered salt (SHA-256 126d2695...f419 in PREREGISTRATION.md) is:

    solver-d25bec698ee055604b636bead82dd761

Check: `printf '%s' 'solver-d25bec698ee055604b636bead82dd761' | sha256sum`.

## Caveats (repeat of PREREGISTRATION.md)

The author had seen all the data before building the model, so this is not
a clean held-out test. Counting the frozen galaxy constants, the model has
17 fitted constants (the scorecard's "9 global" counts only those the
benchmark fits).

## Hidden-split scorecard

| # | Test | Perfect theory | Collapsing self-interacting dark matter (cSIDM) | Cold dark matter (ΛCDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|---|---|
| 01 | Rotation curves: predict | 95 (95–96) | 56 | 52 | 57 | 50 | 33 | 33 |
| 02 | Rotation curves: extrapolate | 95 (94–97) | 76 | 40 | 79 | 70 | 19 | 19 |
| 03 | Galaxy clusters | 96 (85–100) | 100 | 100 | 15 | 38 | 0 | 0 |
| 04 | Milky Way | 96 (90–100) | 100 | 87 | 100 | 13 | 0 | 0 |
| 05 | Dwarf galaxies | 96 (90–100) | 58 | 47 | 42 | 42 | 21 | 21 |
| 07 | Lensing around galaxies ‡ | 95 (88–100) | 88 | 50 | 100 | 100 | 0 | 0 |
| 08 | Elliptical galaxies (lensing) | 95 (91–100) | 88 | 98 | 79 | 84 | 74 | 74 |
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
| 18 | Supernova distances | 96 (94–97) | 99 | 99 | — | — | 98 | — |
| 21 | Galaxy cluster counts ‡ | 96 (75–100) | 75 | 75 | — | — | 0 | — |
| 28 | Inflation signatures ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 29 | Neutrino masses and species ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 30 | Curvature and topology ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 31 | Growth rate of structure ‡ | 96 (75–100) | 100 | 100 | — | — | 0 | — |
| 34 | Standard sirens ◆ ‡ | 95 (50–100) | 100 | 100 | — | — | 100 | — |
| 24 | Direct detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 25 | Indirect detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 36 | Cosmic antimatter ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| | **Lowest** | | **56** | **40** | **0** | **0** | **0** | **0** |
| | **Substantive tests passed** (of 14) | | **10** | **9** | **3** | **2** | **1** | **0** |
| | Median (excl. guard rails) | | 94 | 92 | 29 | 26 | 0 | 0 |
| | Highest | | 100 | 100 | 100 | 100 | 100 | 100 |
| | Knobs | | 9 global, 1 per galaxy | 10 global, 3 per galaxy | 3 global, 1 per galaxy | 4 global, 1 per galaxy | 6 global, 1 per galaxy | 2 global, 1 per galaxy |
| | **Rank** | | 1 | 2 | 3 | 4 | 5 | 6 |

## Practice scorecard (public data only)

| # | Test | Perfect theory | Collapsing self-interacting dark matter (cSIDM) | Cold dark matter (ΛCDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|---|---|
| 01 | Rotation curves: predict | 95 (94–97) | 61 | 49 | 66 | 66 | 28 | 28 |
| 02 | Rotation curves: extrapolate | 95 (93–97) | 68 | 32 | 82 | 79 | 25 | 25 |
| 03 | Galaxy clusters | 95 (88–100) | 62 | 62 | 12 | 50 | 0 | 0 |
| 04 | Milky Way | 96 (90–100) | 100 | 87 | 100 | 76 | 0 | 0 |
| 05 | Dwarf galaxies | 96 (88–100) | 71 | 71 | 29 | 35 | 6 | 6 |
| 07 | Lensing around galaxies ‡ | 95 (88–100) | 88 | 50 | 100 | 100 | 0 | 0 |
| 08 | Elliptical galaxies (lensing) | 95 (88–100) | 94 | 82 | 59 | 71 | 59 | 59 |
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
| 18 | Supernova distances | 95 (94–97) | 100 | 100 | — | — | 98 | — |
| 21 | Galaxy cluster counts ‡ | 96 (75–100) | 75 | 75 | — | — | 0 | — |
| 28 | Inflation signatures ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 29 | Neutrino masses and species ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 30 | Curvature and topology ◆ ‡ | 98 | 100 | 100 | — | — | 100 | — |
| 31 | Growth rate of structure ‡ | 96 (75–100) | 100 | 100 | — | — | 0 | — |
| 34 | Standard sirens ◆ ‡ | 95 (50–100) | 100 | 100 | — | — | 100 | — |
| 24 | Direct detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 25 | Indirect detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 36 | Cosmic antimatter ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| | **Lowest** | | **61** | **32** | **0** | **0** | **0** | **0** |
| | **Substantive tests passed** (of 14) | | **10** | **7** | **3** | **2** | **1** | **0** |
| | Median (excl. guard rails) | | 91 | 84 | 21 | 43 | 0 | 0 |
| | Highest | | 100 | 100 | 100 | 100 | 100 | 100 |
| | Knobs | | 9 global, 1 per galaxy | 10 global, 3 per galaxy | 3 global, 1 per galaxy | 4 global, 1 per galaxy | 6 global, 1 per galaxy | 2 global, 1 per galaxy |
| | **Rank** | | 1 | 2 | 3 | 4 | 5 | 6 |
