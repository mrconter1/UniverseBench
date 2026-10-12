# Candidate tests: results log

## A. Is the galaxy relation a relaxation attractor? (12 Oct 2026)

Script: `analysis/rar_age.py` (all 165 SPARC galaxies, RAR with a0 = 1.2e-10,
M/L 0.5 disc / 0.7 bulge). Prediction if A is right: less evolved galaxies
(gas-rich, low surface brightness) scatter MORE around the relation than
evolved ones. A fundamental force law predicts no difference.

| Scatter measure | Split by | Less evolved half | More evolved half | Difference (95% CI) |
|---|---|---|---|---|
| Between galaxies (offsets) | surface brightness | 0.22 dex | 0.14 dex | −0.08 (−0.15..−0.02) |
| Between galaxies (offsets) | baryonic mass | 0.24 | 0.12 | −0.12 (−0.17..−0.06) |
| Between galaxies (offsets) | gas fraction | 0.21 | 0.17 | −0.04 (−0.11..+0.03) |
| Within galaxies (shape) | gas fraction | 0.00 | 0.07 | opposite sign |
| Within galaxies (shape) | surface brightness | 0.00 | 0.07 | opposite sign |

Between galaxies the trend goes the predicted way, but galaxy-wide distance and
inclination errors (not in our files, larger for faint dwarfs) could produce
it. Within galaxies, which is immune to those errors, less evolved galaxies
show **no** intrinsic scatter. The evolved ones show some, most likely from a
fixed stellar M/L in bulge+disc systems.

**Verdict: not supported, not ruled out.** Deciding test: marginalise
per-galaxy distance and inclination uncertainties (SPARC Table 1; Li et al.
2018 method) and repeat the between-galaxy comparison.

## I. Does an "age clock" make a0 rise with redshift? (12 Oct 2026)

Script: `analysis/a0_redshift.py`. The frozen cSIDM model on all SPARC baryon
distributions, halos placed at redshift 1 (galaxy–halo connection held
fixed, an assumption); effective a0 fitted to the model's own accelerations.

| Case | Effective a0 | Ratio to z = 0 |
|---|---|---|
| z = 0 reference | 0.86×10⁻¹⁰ m/s² | 1.00 |
| Clock only (halo age 10 → 4.3 Gyr) | 0.76×10⁻¹⁰ | 0.88 |
| Density only (r200 at ρ_crit(z), Dutton & Maccio c(z)) | 1.04×10⁻¹⁰ | 1.21 |
| Clock + density (redshift 1) | 0.98×10⁻¹⁰ | **1.14** |

The clock pushes a0 *down* (less time to scatter), the opposite of the rough
scaling argued beforehand. **Verdict: the clock idea fails as an explanation
of a rising a0.** It does give a falsifiable discriminator: a0(z=1)/a0(0) is
1.0 for MOND, about 1.1–1.2 for cSIDM, about 1.8 for horizon-type theories
(a0 ∝ H), and about 2.5 as measured by MUSE-DARK III. Nothing tested predicts
2.5, so stellar-mass systematics at high redshift are the first suspect.
