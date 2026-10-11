"""
Build the folder you give to an AI agent: the package, the instructions,
and only PUBLIC data. Hidden data never leaves your machine.

    BENCH_SALT='my-private-secret' python make_agent_bundle.py

Inside the bundle, scoring uses a fixed practice split of the public data,
so an agent can develop and self-check without ever seeing hidden data.
"""
import os
import shutil

from ubench import core
from ubench.tests import all_tests

ROOT = core.ROOT
out = os.path.join(ROOT, "agent_bundle")
shutil.rmtree(out, ignore_errors=True)
shutil.copytree(os.path.join(ROOT, "ubench"), os.path.join(out, "ubench"),
                ignore=shutil.ignore_patterns("__pycache__"))
for f in ("score.py", "AGENT_TASK.md", "LEADERBOARD.md", "ROADMAP.md"):
    if os.path.exists(os.path.join(ROOT, f)):
        shutil.copy(os.path.join(ROOT, f), out)

# the agent's copy uses a fixed, non-secret practice split
core_py = os.path.join(out, "ubench", "core.py")
src = open(core_py).read().replace('os.environ.get("BENCH_SALT", "demo-salt")',
                                   '"agent-practice-split"')
open(core_py, "w").write(src)

dest = os.path.join(out, "data")
os.makedirs(dest)
for t in all_tests():
    t.ensure_loaded()
    t.export_public(dest, core.SPLIT_SALT)
    print(f"  {t.ID}: exported {'public part' if t.HIDDEN_SPLIT else 'all'} of its data")
print(f"agent_bundle/ written (hidden split salt from BENCH_SALT)")
