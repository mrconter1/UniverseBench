# Pre-registration: solver model 2, cold dark matter + galaxy formation

Written 12 Oct 2026, before any score of this model on its hidden split.

* **Hypothesis:** ordinary cold dark matter, Λ and general relativity; the
  galaxy-scale regularities come from the halo response to galaxy formation
  as measured in hydrodynamic simulations (DC14 profile, Di Cintio et al.
  2014b), plus a fitted halo connection. Motivation: analysis/RESULTS.md,
  "Is the tightness of the galaxy relation really a problem for cold dark
  matter?".
* **Development:** on the public practice split (`demo-salt`), which has been
  seen before. Any change after the hidden run is reported as post-hoc.
* **Hidden run:** one shot on a fresh salt whose SHA-256 is
  `fdd10d33244173f762aefa04a77b018db43bb0a8eab75e8606c3a427ee122a96`.
  The salt is revealed in SCORES.md after the run.
* **Comparison:** the frozen cSIDM solver model, ΛCDM and MOND on the same
  hidden split.
* **Caveat:** as for the first solver run, the author has seen all the data;
  this is not a clean held-out test, only a guard against tuning on the
  scored points.
