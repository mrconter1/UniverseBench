"""
Build the folder you give to an AI agent: public galaxies only, plus the
model interface and instructions. The hidden galaxies never leave your
machine.

    BENCH_SALT='my-private-secret' python make_agent_bundle.py
"""
import os
import shutil

import benchmark as B

out = os.path.join(B.HERE, "agent_bundle")
shutil.rmtree(out, ignore_errors=True)
os.makedirs(os.path.join(out, "data", "all"))

gals = B.load_all()
pub, hid = B.split(gals)
for g in pub:
    shutil.copy(os.path.join(B.HERE, "data", "all", f"{g.name}_rotmod.dat"),
                os.path.join(out, "data", "all"))
for f in ("models.py", "clusters.py", "AGENT_TASK.md"):
    shutil.copy(os.path.join(B.HERE, f), out)

# clusters: copy only rows for public clusters
os.makedirs(os.path.join(out, "data", "clusters"))
src_csv = os.path.join(B.HERE, "data", "clusters", "vikhlinin2006.csv")
cpub, chid = B.split_clusters(B.clusters.load_clusters())
pub_names = {c.name for c in cpub}
with open(src_csv) as fi, open(os.path.join(out, "data", "clusters", "vikhlinin2006.csv"), "w") as fo:
    for line in fi:
        if line.startswith(("#", "name,")) or line.split(",")[0] in pub_names:
            fo.write(line)
print(f"clusters: {len(pub_names)} public; {len({c.name for c in chid})} kept hidden")
# The agent's copy splits its public galaxies into a practice train/test
# split with a fixed, non-secret salt.
src = open(os.path.join(B.HERE, "benchmark.py")).read()
src = src.replace('os.environ.get("BENCH_SALT", "demo-salt")', '"agent-practice-split"')
open(os.path.join(out, "benchmark.py"), "w").write(src)

print(f"agent_bundle/: {len(pub)} public galaxies; {len(hid)} kept hidden")
