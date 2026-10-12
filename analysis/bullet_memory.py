"""Idea 4: memory-time limit from the Bullet Cluster offset (data file, not memory).
If gravity responds to where matter was over a memory time tau, the
sub-cluster's lensing peak trails its galaxies by about v * tau."""
import csv
r = {k: float(v) for k, v in next(csv.DictReader(l for l in open("data/bullet_cluster/randall2008_offset.csv") if not l.startswith("#"))).items()}
lim = r["offset_kpc"] + 2 * r["offset_err_kpc"]                    # 2 sigma upper limit on the lag
KPC_KM, MYR_S = 3.0857e16, 3.156e13
for v in (r["velocity_kms"] - r["velocity_err_kms"], r["velocity_kms"], 3000.0):
    print(f"v = {v:6.0f} km/s: lag < {lim:.0f} kpc (2 sigma)  =>  memory time tau < {lim * KPC_KM / v / MYR_S:.0f} Myr")
