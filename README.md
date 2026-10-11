# Galaxy Benchmark

Can a theory of gravity or dark matter predict, from visible matter alone,
how 165 real galaxies rotate and how massive 13 galaxy clusters are, on data
it has never seen? See ROADMAP.md for the 25 planned tests.

Data: SPARC rotation curves, Lelli, McGaugh & Schombert (2016), AJ 152, 157
(175 galaxies; 165 with at least 6 points are used). Mirrored from
github.com/villano-lab/galactic-spin-W1.
Clusters: Vikhlinin et al. (2006), ApJ 640, 691, Tables 1 and 3 (Chandra;
13 clusters, 23 mass measurements at r2500 and r500).

## Files

- `benchmark.py`: data loading, secret hidden split, the tests, scoring
- `clusters.py`: cluster data and test C
- `models.py`: baselines (Newton, MOND, cold dark matter with NFW halos)
- `score.py`: leaderboard on hidden galaxies; pass candidate model files
- `make_agent_bundle.py`: builds `agent_bundle/` with public galaxies only
- `AGENT_TASK.md`: instructions for an AI agent

## Workflow

    pip install numpy scipy
    export BENCH_SALT='your-private-secret'    # decides the hidden half
    python make_agent_bundle.py                # give agent_bundle/ to an agent
    # ...agent writes my_model.py...
    python score.py agent_bundle/my_model.py   # score on hidden galaxies

## Scoring

Per-point uncertainty is the quoted error plus a 5% floor for distance and
inclination uncertainty. Score = % of hidden points within 2 sigma
(perfect is about 95). Test A predicts whole curves with global
parameters only; test B tunes per-galaxy knobs on the inner 60% and
predicts the outer 40%. Test C predicts cluster total masses from gas plus
stars (15% quoted error floor). Total = average. Global constants are fitted
jointly on the public data of all tests.

## Known simplifications

- No SPARC quality or inclination cuts (the galaxy table was unavailable).
- Stellar mass for the dark matter model is estimated by integrating surface
  brightness over measured radii, which underestimates it somewhat.
- NFW halos have no baryonic feedback (cores), which would help the dark
  matter model in dwarf galaxies.
- Cluster stellar mass is a fixed 15% of gas mass, not measured per cluster.
- The dark matter cluster model uses an empirical baryon fraction with mass
  and radius trends rather than simulated halos.
- Only 13 clusters, so the cluster score moves noticeably between splits.
- No solar system, lensing or cosmology tests yet.
