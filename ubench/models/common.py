"""Pieces shared by the baseline models."""
import numpy as np

UPS = {"ups_disk": (0.5, 0.1, 1.5)}           # global mass-to-light (3.6 um)
LOCAL_UPS = {"ups_scale": (1.0, 0.33, 3.0)}   # per-galaxy correction


def v_baryon_sq(g, ups_disk, ups_bul):
    """Newtonian speed^2 from visible matter. Signed squares: SPARC's gas
    term can be negative in central holes (net outward pull). Galaxies
    with ups_fixed (the Milky Way) already carry real stellar masses."""
    if getattr(g, "ups_fixed", False):
        ups_disk = ups_bul = 1.0
    s = lambda v: v * np.abs(v)
    return s(g.Vgas) + ups_disk * s(g.Vdisk) + ups_bul * s(g.Vbul)
