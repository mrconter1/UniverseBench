"""
Build data/dwarfs/lvdb_dwarfs.csv from the Local Volume Database
(Pace 2024, OJAp; github.com/apace7/local_volume_database, CC0).

    git clone https://github.com/apace7/local_volume_database
    python data/dwarfs/build_lvdb.py path/to/local_volume_database

Selection (all decided before scoring any model):
  * confirmed galaxies in dwarf_mw.csv, dwarf_m31.csv, dwarf_local_field.csv
    with a stellar velocity dispersion, either measured or a one-sided
    upper limit (95% in the database). Upper limits are KEPT: dropping them
    removes exactly the faintest systems, where the measured dispersion is
    smallest, and biases every sample upward.
  * finite stellar mass and half-light radius.
  * measured dispersions with relative error <= 0.5 (larger ones carry no
    information at 2 sigma).
  * excluded: discs and compact ellipticals that are not dispersion-
    supported dwarfs (LMC, SMC, M32, NGC 205, and anything with
    M* > 10^8.5 Msun), Sagittarius (being disrupted), and systems with
    reported tidal disruption or stripping (TIDAL below).
No cut on the dispersion itself.
"""
import csv
import math
import subprocess
import sys
from pathlib import Path

EXCLUDE = {"LMC", "SMC", "Sagittarius", "M 32", "NGC 205"}
MAX_LOG_MSTAR = 8.5
# Reported tidal disruption or stripping (tidal tails, S-shaped isophotes,
# velocity gradients): Hercules (Kuepper+ 2017), Willman 1 (Willman+ 2011),
# Antlia II and Crater II (Torrealba+ 2016/2019; Ji+ 2021), Andromeda XXVII
# (Collins+ 2013), NGC 147 (Crnojevic+ 2014), Tucana III (Drlica-Wagner+
# 2015), Bootes III (Carlin+ 2009), Ursa Major II (Munoz+ 2010),
# Segue 2 (Kirby+ 2013, possibly stripped).
TIDAL = {"Hercules", "Willman 1", "Antlia II", "Crater II", "Andromeda XXVII",
         "NGC 147", "Tucana III", "Bootes III", "Ursa Major II", "Segue 2"}
FILES = [("dwarf_mw.csv", "MW"), ("dwarf_m31.csv", "M31"), ("dwarf_local_field.csv", "field")]


def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def main(lvdb):
    lvdb = Path(lvdb)
    commit = subprocess.run(["git", "-C", str(lvdb), "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip() or "unknown"
    out, dropped = [], []
    for fn, host in FILES:
        for r in csv.DictReader(open(lvdb / "data" / fn)):
            if r["confirmed_galaxy"] != "1":
                continue
            name = r["name"]
            s, s_ul = num(r["vlos_sigma"]), num(r["vlos_sigma_ul"])
            lm, rh = num(r["mass_stellar"]), num(r["rhalf_sph_physical"])
            if s is None and s_ul is None:
                continue
            if lm is None or rh is None:
                dropped.append((name, "no stellar mass or radius")); continue
            if name in EXCLUDE or lm > MAX_LOG_MSTAR:
                dropped.append((name, "not a dispersion-supported dwarf")); continue
            if name in TIDAL:
                dropped.append((name, "tidal disruption reported")); continue
            em, ep = num(r["vlos_sigma_em"]), num(r["vlos_sigma_ep"])
            if s is not None:
                if em is None or ep is None:
                    dropped.append((name, "no dispersion error")); continue
                if 0.5 * (em + ep) / s > 0.5:
                    dropped.append((name, "dispersion error > 50%")); continue
            d_host = num(r["distance_host"])
            if d_host is None and host != "field":
                d_host = num(r["distance_gc" if host == "MW" else "distance_m31"])
            if d_host is None and host != "field":
                dropped.append((name, "no host distance")); continue
            near, d_near = None, None
            if host == "field":
                # nearest giant for the external field: MW or M31 (or the
                # listed host, e.g. NGC 3109 for Antlia B)
                cands = [("MW", num(r["distance_gc"])), ("M31", num(r["distance_m31"]))]
                if r.get("host") and d_host is not None:
                    cands.append((r["host"], d_host))
                cands = [c for c in cands if c[1] is not None]
                near, d_near = min(cands, key=lambda c: c[1]) if cands else (None, None)
            mhi = num(r["mass_HI"])
            out.append({
                "name": name, "host": host,
                "nearest_giant": near or host,
                "D_host": round(d_host if host != "field" else d_near, 1),
                "sigma": "" if s is None else s, "sigma_em": "" if s is None else em,
                "sigma_ep": "" if s is None else ep, "sigma_ul": "" if s_ul is None else s_ul,
                "log_Mstar": lm, "log_MHI": "" if mhi is None else mhi,
                "r_half_pc": round(4.0 / 3.0 * rh, 1),     # 3D half-light radius (Wolf+ 2010)
                "ref_vlos": r["ref_vlos"],
            })
    dest = Path(__file__).with_name("lvdb_dwarfs.csv")
    with open(dest, "w", newline="") as f:
        f.write(f"# Local Volume Database (Pace 2024, OJAp 7, 111; CC0), commit {commit}\n")
        f.write("# Built by build_lvdb.py; see its docstring for the selection. sigma = stellar line-of-sight\n")
        f.write("# velocity dispersion (km/s) with errors em/ep, or sigma_ul = 95% upper limit (one-sided).\n")
        f.write("# r_half_pc = 3D half-light radius = 4/3 x projected (Wolf+ 2010). log_MHI = HI mass if detected.\n")
        f.write("# D_host = distance to host (satellites) or to the nearest giant (field dwarfs), kpc.\n")
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} dwarfs written ({sum(1 for o in out if o['sigma_ul'] != '')} upper limits); "
          f"{len(dropped)} dropped:")
    for n, why in dropped:
        print(f"   {n:20s} {why}")


if __name__ == "__main__":
    main(sys.argv[1])
