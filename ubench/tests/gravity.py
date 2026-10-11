"""
Tests 11-14: precision gravity guard rails.

General relativity passes all of these; they rule out whole classes of
modified gravity. Each is a short table of measured quantities
(data/gravity/measurements.csv) that a model must predict:

  11 Solar system: PPN gamma (Cassini), PPN beta (MESSENGER), Nordtvedt
     eta and dG/dt / G (lunar laser ranging).
  12 Binary pulsars: orbital decay of the Hulse-Taylor pulsar and the
     double pulsar's gravitational-wave emission, relative to GR.
  13 Gravitational waves: speed of gravity (GW170817) and graviton mass
     (GWTC-3).
  14 Laboratory: equivalence principle (MICROSCOPE) and the inverse-square
     law at 40 micrometres (Eot-Wash).
  26 Black holes: the size of Sgr A*'s shadow (EHT), M87*'s shadow mass vs
     its stellar-dynamics mass, and the orbit of the star S2 around Sgr A*
     (relativistic precession and gravitational redshift, GRAVITY).
  27 Clocks: gravitational time dilation measured with Galileo satellites
     on eccentric orbits and optical clocks on Tokyo Skytree.

Bounds are converted to Gaussian uncertainties (see the CSV notes). Too few
numbers to hide: scored in-sample.

Model method: predict_gravity(quantity, gp) -> value (same units as the
CSV), or None if the model makes no prediction. Models built on GR can
inherit models.common.GeneralRelativity.
"""
import csv

import numpy as np

from ..core import Test, data_path


def load_gravity():
    path = data_path("gravity", "measurements.csv")
    return list(csv.DictReader(l for l in open(path) if not l.startswith("#")))


class _GravityTable(Test):
    HIDDEN_SPLIT = False
    STAGE = 2
    DATA_DIRS = ("gravity",)
    NUMBER = ""

    def load(self):
        self.rows = [r for r in load_gravity() if r["test"] == self.NUMBER]

    def residuals(self, model, gp, part, salt):
        if not hasattr(model, "predict_gravity"):
            return None
        z = []
        for r in self.rows:
            v = model.predict_gravity(r["quantity"], gp)
            if v is None:
                return None
            z.append((v - float(r["value"])) / float(r["sigma"]))
        return np.array(z)


class SolarSystem(_GravityTable):
    ID, LABEL, NUMBER = "11_solar_system", "Solar system", "11"


class BinaryPulsars(_GravityTable):
    ID, LABEL, NUMBER = "12_binary_pulsars", "Binary pulsars", "12"


class GravitationalWaves(_GravityTable):
    ID, LABEL, NUMBER = "13_gravitational_waves", "Gravitational waves", "13"


class LabGravity(_GravityTable):
    ID, LABEL, NUMBER = "14_lab_gravity", "Laboratory gravity", "14"


class BlackHoles(_GravityTable):
    ID, LABEL, NUMBER = "26_black_holes", "Black holes", "26"


class Clocks(_GravityTable):
    ID, LABEL, NUMBER = "27_clocks", "Clocks and time dilation", "27"
