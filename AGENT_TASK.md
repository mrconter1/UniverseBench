# Task: explain how galaxies rotate

You are given rotation curves for a set of real galaxies (SPARC catalog) in
`data/all/`. Each file lists, at several radii R (kpc):

- `Vobs`, `errV`: the measured orbital speed of gas and stars (km/s) and its error
- `Vgas`, `Vdisk`, `Vbul`: the speed that the *visible* gas, stellar disk and
  bulge would produce alone under ordinary Newtonian gravity. Disk and bulge
  assume a mass-to-light ratio of 1; scale their squares by your chosen
  mass-to-light ratio.
- `SBdisk`, `SBbul`: surface brightness (Lsun/pc^2)

You are also given galaxy clusters in `data/clusters/` (loaded by
`clusters.py`): for each cluster and radius, the visible (gas + stellar)
mass inside that radius and the measured total mass.

Visible matter alone cannot explain the observed speeds or cluster masses.
Your job is to propose a physical model, dark matter, modified gravity or
something new, that predicts both from the visible matter.

## Rules

1. Implement your model in a file `my_model.py` with the same interface as the
   classes in `models.py`, ending with `MODEL = YourClass()`.
2. Your model will be scored on **galaxies you have not seen**. You only have
   the public half. Overfitting the public galaxies will hurt you.
3. `GLOBAL_PARAMS` are shared by every galaxy and are tuned automatically on
   the public galaxies. `LOCAL_PARAMS` are per-galaxy knobs; their defaults
   are used in the pure prediction test.
4. Three tests, averaged into a total:
   - **Predict**: global parameters only; predict whole rotation curves.
   - **Extrapolate**: your per-galaxy knobs are tuned on the inner 60% of each
     hidden galaxy; predict the outer 40%. Extra knobs only help if they
     capture real physics.
   - **Clusters**: implement `predict_cluster_mass(c, gp)` returning total
     mass (Msun) inside `c.r_kpc`, given visible mass `c.M_bar`. A model
     without this method scores 0 on this test.
5. **One set of constants for everything.** Global parameters are fitted
   jointly on galaxies and clusters. A constant that appears in both must take
   one value.
6. The model must be physically meaningful: no lookup tables of galaxy or
   cluster names, no using observed speeds or total masses inside your
   predictions, and it must reduce to ordinary Newtonian gravity where
   gravity is strong (the solar system).
7. Use `benchmark.run(MODEL)` on the public data to develop. The score you see
   there uses a practice split of the public data and is NOT your final score.

## Baselines to beat (hidden data)

| Model | Predict | Extrapolate | Clusters | Total |
|---|---|---|---|---|
| Newton, visible matter only | ~31 | ~13 | 0 | ~15 |
| MOND (one universal a0) | ~58 | ~75 | ~43 | ~57 |
| Cold dark matter (NFW halo) | ~46 | ~41 | ~95 | ~61 |

MOND wins on galaxies; dark matter wins on clusters. Nothing yet wins both.

A perfect model would score about 95.

Report: your model's physical idea in plain language, its equations, how
many parameters it uses, and what new observation could prove it wrong.
