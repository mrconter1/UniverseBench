"""Baseline models. A model is any object with NAME, GLOBAL_PARAMS
({name: (start, lo, hi)}), optional LOCAL_PARAMS, and one prediction method
per test (see each test module's docstring for the method it needs)."""
from .lcdm import LCDM
from .mond import MOND, MONDFixed
from .visible import VisibleOnly


def baselines():
    return [VisibleOnly(), MOND(), MONDFixed(), LCDM()]
