# Pre-registration: Claude as solver

Written before any model development (12 Oct 2026).

## Final scoring split

The final hidden split is set by a secret salt stored outside the repo.
SHA-256 of the salt:

    126d2695b7445768f92ceb6774b4434a167c404f9655e28d38ceb167a671f419

The salt is revealed only after the final model is frozen, so anyone can
check that the reported hidden scores used this split.

## Protocol

* Development happens only in an agent bundle built with this salt, which
  holds the public half of every split test. Fitting and practice scoring
  use only that bundle.
* The final model is frozen (file checksum recorded) and then scored once on
  the hidden half, and again over the five demo splits for the leaderboard.
* Everything tried, including failures, is logged in this folder.

## Known contamination (honest statement)

The developer also built the benchmark. It has seen every data file,
including the hidden halves, and earlier in this project it analysed all
benchmark dwarfs (old and modern) by name. The hidden-split score is
therefore **not** a clean held-out test, only a check against accidental
overfitting to the public half. The clean tests of anything built here are
its predictions for data that do not exist yet (see THEORY.md once written)
and scoring on future data releases.

## Starting hypothesis (from the pilot and its follow-ups)

Dark matter is a real, cold particle (cosmology, clusters, Bullet Cluster
and dwarf satellites need it). The galaxy regularity (radial acceleration
relation) emerges where dark matter has thermalised in the baryonic
potential through velocity-dependent self-interactions (SIDM; Kaplinghat
et al. 2014, Ren et al. 2019), not from a new force. The first candidate is
therefore velocity-dependent SIDM with semi-analytic isothermal cores
matched to NFW halos, plus standard galaxy-formation inputs. Alternatives
will be tried if this fails.
