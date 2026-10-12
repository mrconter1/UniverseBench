"""
Test 37 (guard rail): is the theory's action healthy, and do its declared
backgrounds obey it?

A model declares `ACTION` (see ubench/action_check.py). Particle dark matter
in general relativity declares {"type": "GR+minimal"} and passes by
construction. Scalar-tensor theories declare their Horndeski functions and
background snapshots, which are checked for ghosts, gradient instabilities,
the GW170817 speed of gravity, and consistency with the Friedmann constraint
derived from the action itself. A model with no checkable action (e.g. MOND
without a declared relativistic completion) makes no prediction (0).

Assumes faithful agents: the action must be the one the model's predictions
come from. The background-consistency check ties the two together on
cosmological scales; galaxy-scale consistency is not checked automatically.
"""
import numpy as np

from ..core import Test
from ..action_check import check


class ActionHealth(Test):
    ID = "37_action_health"
    LABEL = "Healthy action"
    STAGE = 2
    HIDDEN_SPLIT = False
    USED_IN_FIT = False
    GUARD_RAIL = True
    N_BOUND = 10 ** 6            # pass/fail checks: a perfect theory passes them all

    def residuals(self, model, gp, part, salt):
        act = getattr(model, "ACTION", None)
        status, report = check(act)
        if status is None:
            return None
        return np.array([0.0 if ok else 1e6 for (_, name, _, ok) in report if not name.startswith("warning")])
