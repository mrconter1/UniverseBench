# Scores: solver model 2, cold dark matter + galaxy formation

Frozen model: `my_model.py`, MD5 30c05b8fff96f0c79cfc7a5c88d9ea9a (frozen before the hidden run,
no changes after the practice run).

## Result (one-shot hidden split, 12 Oct 2026)

On every galaxy test except the Milky Way, plain cold dark matter with
simulation-calibrated halo shapes (DC14) **ties the self-interacting model**,
and it roughly doubles standard ΛCDM where ΛCDM was weakest:

| Test | New model | cSIDM | ΛCDM | MOND |
|---|---|---|---|---|
| Rotation curves, predict | 58 | 58 | 55 | 58 |
| Rotation curves, extrapolate | **75** | 72 | 36 | 77 |
| Dwarf galaxies | 64 | 64 | 58 | 36 |
| Lensing around galaxies | 88 | 88 | 50 | 100 |
| Ellipticals | 87 | 87 | 90 | 74 |
| Milky Way (pure out-of-sample) | **53** | 100 | 87 | 100 |

Clusters, the Bullet Cluster and cosmology are the same for all three dark
matter models.

**Ranking caveat, stated plainly:** the scorecard ranks the new model first
(lowest 50, 9 passes) and cSIDM third (lowest 0), but cSIDM's 0 is test 37
(healthy action), which did not exist when cSIDM was frozen, so it declares
no action. It is GR plus minimally coupled dark matter and would pass.
Without that bookkeeping, cSIDM has lowest 50, **10 passes, median 92**,
against the new model's 50 / 9 / 87: **cSIDM edges it, entirely through the
Milky Way.** Both models' lowest score (50) is cluster counts, shared with
ΛCDM under the new joint cosmology fit.

**Milky Way:** the fitted halo connection gives our galaxy a ~1.7 × 10¹² Msun
halo, contracted (c × 1.7) by the DC14 recipe, so the rotation curve is
over-predicted at 15–25 kpc. Not tuned: the Milky Way is never used to fit.

## Salt reveal

Pre-registered SHA-256 (PREREGISTRATION.md):
`fdd10d33244173f762aefa04a77b018db43bb0a8eab75e8606c3a427ee122a96`. Salt:

    solver-galform-0f817a36566eb5988728d1799469865d

Check: `printf '%s' 'solver-galform-0f817a36566eb5988728d1799469865d' | sha256sum`.

## Practice split (demo-salt), for reference

Rotation 55, extrapolation 76, Milky Way 66, dwarfs 64, lensing 88, ellipticals 95.

## Hidden-split scorecard


| # | Test | Perfect theory | Cold dark matter + galaxy formation (DC14 halos) | Cold dark matter (ΛCDM) | Collapsing self-interacting dark matter (cSIDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|---|---|---|
| 01 | Rotation curves: predict | 95 (95–96) | 58 | 55 | 58 | 58 | 50 | 30 | 30 |
| 02 | Rotation curves: extrapolate | 96 (94–97) | 75 | 36 | 72 | 77 | 66 | 26 | 26 |
| 03 | Galaxy clusters | 95 (86–100) | 100 | 100 | 100 | 21 | 50 | 0 | 0 |
| 04 | Milky Way | 96 (90–100) | 53 | 87 | 100 | 100 | 5 | 0 | 0 |
| 05 | Dwarf galaxies | 96 (89–100) | 64 | 58 | 64 | 36 | 42 | 17 | 17 |
| 07 | Lensing around galaxies ‡ | 95 (88–100) | 88 | 50 | 88 | 100 | 100 | 0 | 0 |
| 08 | Elliptical galaxies (lensing) | 95 (87–100) | 87 | 90 | 87 | 74 | 77 | 74 | 74 |
| 09 | Bullet Cluster § ‡ | 95 (50–100) | 100 | 100 | 100 | 0 | 0 | 0 | 0 |
| 11 | Solar system ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | 100 | 100 | 25 |
| 12 | Binary pulsars ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 100 | 0 |
| 13 | Gravitational waves ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 100 | — |
| 14 | Laboratory gravity ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 26 | Black holes ◆ ‡ | 96 (75–100) | 100 | 100 | 100 | 100 | 100 | 100 | — |
| 27 | Clocks and time dilation ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 100 | 0 |
| 37 | Healthy action ◆ | 100 | 100 | 100 | — | — | — | 100 | — |
| 15 | Microwave background (Planck spectra) ‡ | 96 (94–97) | 96 | 96 | 96 | — | — | 0 | — |
| 16 | Big Bang nucleosynthesis ‡ | 95 (50–100) | 100 | 100 | 100 | 100 | 100 | 0 | — |
| 17 | Baryon acoustic oscillations ‡ | 96 (85–100) | 85 | 85 | 85 | — | — | 8 | — |
| 18 | Supernova distances | 95 (94–97) | 99 | 99 | 99 | — | — | 98 | — |
| 21 | Galaxy cluster counts ‡ | 96 (75–100) | 50 | 50 | 50 | — | — | 0 | — |
| 28 | Inflation signatures ◆ ‡ | 98 | 100 | 100 | 100 | — | — | 100 | — |
| 29 | Neutrino masses and species ◆ ‡ | 98 | 100 | 100 | 100 | — | — | 100 | — |
| 30 | Curvature and topology ◆ ‡ | 98 | 100 | 100 | 100 | — | — | 100 | — |
| 31 | Growth rate of structure ‡ | 96 (75–100) | 100 | 100 | 100 | — | — | 0 | — |
| 34 | Standard sirens ◆ ‡ | 95 (50–100) | 100 | 100 | 100 | — | — | 100 | — |
| 24 | Direct detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 25 | Indirect detection ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| 36 | Cosmic antimatter ◆ ‡ | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| | **Lowest** | | **50** | **36** | **0** | **0** | **0** | **0** | **0** |
| | **Substantive tests passed** (of 14) | | **9** | **8** | **10** | **3** | **2** | **1** | **0** |
| | Median (excl. guard rails) | | 87 | 89 | 92 | 29 | 24 | 0 | 0 |
| | Highest | | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| | Knobs | | 15 global, 1 per galaxy | 11 global, 3 per galaxy | 10 global, 1 per galaxy | 3 global, 1 per galaxy | 4 global, 1 per galaxy | 7 global, 1 per galaxy | 2 global, 1 per galaxy |
| | **Rank** | | 1 | 2 | 3 | 4 | 5 | 6 | 7 |

Scores: % of data points predicted within 2 sigma. 'Perfect theory' is the noise ceiling: the expected score of the true laws of nature given each test's error bars (with the 5th-95th percentile range for small tests). A test is passed when the score reaches the 5th percentile of that range. The benchmark uses settled tests only; models are ranked by their lowest score (a theory is only as good as its worst failure), then substantive tests passed, then median. ◆ = guard rail: a gate every viable theory must pass, counted in the lowest score but not the median. § = simplified test (09: mass-weighted peak position from aperture masses, not a fit to the lensing map; a faithful model should compute where its lensing mass sits). '—' = no prediction (counts as 0). ‡ = too few numbers to hide: fitted and scored in-sample. 

| # | Open question (not scored) | Cold dark matter + galaxy formation (DC14 halos) | Cold dark matter (ΛCDM) | Collapsing self-interacting dark matter (cSIDM) | MOND (a0 fixed at 1.2, as published) | MOND (a0 fitted to all tests) | Standard gravity, visible matter only | Pure Newtonian gravity, visible matter only |
|---|---|---|---|---|---|---|---|---|
| 06 | Wide binary stars | 50 | 50 | 50 | 50 | 0 | 50 | 50 |
| 10 | Galaxies lacking dark matter | 0 | 0 | 0 | 50 | 0 | 100 | 100 |
| 16b | Lithium problem | 0 | 0 | 0 | 0 | 0 | 0 | — |
| 19 | Hubble constant (local) | 50 | 50 | 50 | — | — | 0 | — |
| 20 | Clumpiness S8 (lensing) | 67 | 67 | 67 | — | — | 0 | — |
| 22 | Lyman-alpha forest | 67 | 67 | 67 | — | — | 33 | — |
| 23 | JWST early galaxies | 100 | 100 | 100 | — | — | 0 | — |
| 29b | Neutrino mass squeeze | 100 | 100 | 100 | — | — | 100 | — |
| 32 | Cosmic dipole (isotropy) | 0 | 0 | 0 | — | — | 0 | — |
| 33 | Dark flow (bulk motions) | 100 | 100 | 100 | — | — | 100 | — |
| 35 | Early supermassive black holes | 0 | 0 | 0 | — | — | 100 | — |

Open questions: tests whose data are disputed or contradict each other, so even a correct theory could not pass them all. Reported, never scored, never used to fit constants.
