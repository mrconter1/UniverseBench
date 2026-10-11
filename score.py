"""
Scorecard: every model on every test.

    python score.py                         # baselines, one split (BENCH_SALT)
    python score.py my_model.py             # + a candidate (defines MODEL = ...)
    python score.py --splits 5 --write      # average over 5 public demo splits,
                                            # write LEADERBOARD.md
    python score.py --tests 01,03,05        # only some tests (by number)

Scores: % of hidden data points predicted within 2 sigma (a perfect theory
scores about 95). Models are summarized by lowest / median / highest score
across tests and ranked by lowest, then median.
"""
import argparse
import importlib.util
import json
import time

import numpy as np

from ubench import core
from ubench.models import baselines
from ubench.tests import all_tests

ap = argparse.ArgumentParser()
ap.add_argument("candidates", nargs="*")
ap.add_argument("--splits", type=int, default=1)
ap.add_argument("--write", action="store_true")
ap.add_argument("--tests", default="")
args = ap.parse_args()

models = baselines()
for path in args.candidates:
    spec = importlib.util.spec_from_file_location("candidate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    models.append(mod.MODEL)

tests = all_tests()
if args.tests:
    keep = {s.zfill(2) for s in args.tests.split(",")}
    tests = [t for t in tests if t.ID[:2] in keep]
salts = [core.SPLIT_SALT] if args.splits == 1 else [f"public-split-{i}" for i in range(args.splits)]

rows = []
for m in models:
    t0 = time.time()
    runs = [core.run(m, tests, salt=s) for s in salts]
    row = {"model": m.NAME, "knobs": f"{runs[0]['n_global_params']} global, "
                                      f"{runs[0]['n_local_params']} per galaxy", "runs": runs}
    for t in tests:
        r0 = runs[0]["tests"][t.ID]
        row[t.ID] = None if r0.get("no_prediction") else \
            float(np.mean([r["tests"][t.ID]["within_2sigma_pct"] for r in runs]))
    # the score uses SETTLED tests only; contested ones go to the open-questions report
    vals = [row[t.ID] or 0.0 for t in tests if t.STATUS == "settled"] or [0.0]
    row["LOWEST"], row["MEDIAN"], row["HIGHEST"] = min(vals), float(np.median(vals)), max(vals)
    rows.append(row)
    print(f"  scored {m.NAME} in {time.time() - t0:.0f}s", flush=True)
rows.sort(key=lambda r: (-r["LOWEST"], -r["MEDIAN"]))

fmt = lambda v: "—" if v is None else f"{v:.0f}"
mark = lambda t: t.LABEL + (" ‡" if not t.HIDDEN_SPLIT and t.USED_IN_FIT else "")
settled = [t for t in tests if t.STATUS == "settled"]
contested = [t for t in tests if t.STATUS != "settled"]

# Wide table: tests as rows, models as columns (reads better with many tests)
head = ["#", "Test"] + [r["model"] for r in rows]
lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
for t in settled:
    lines.append("| " + " | ".join([t.ID.split("_")[0], mark(t)] + [fmt(r[t.ID]) for r in rows]) + " |")
for key in ("LOWEST", "MEDIAN", "HIGHEST"):
    cells = [f"**{r[key]:.0f}**" if key == "LOWEST" else f"{r[key]:.0f}" for r in rows]
    lines.append("| | **" + key.capitalize() + "** | " + " | ".join(cells) + " |")
lines.append("| | Knobs | " + " | ".join(r["knobs"] for r in rows) + " |")
lines.append("| | **Rank** | " + " | ".join(str(i) for i in range(1, len(rows) + 1)) + " |")
table = "\n".join(lines)

oq = ["| # | Open question (not scored) | " + " | ".join(r["model"] for r in rows) + " |",
      "|" + "---|" * (len(rows) + 2)]
for t in contested:
    oq.append("| " + " | ".join([t.ID.split("_")[0], t.LABEL] + [fmt(r[t.ID]) for r in rows]) + " |")
open_table = "\n".join(oq)

note = ("Scores: % of data points predicted within 2 sigma (a perfect theory scores "
        "about 95). The benchmark score uses settled tests only: models are ranked by "
        "their lowest settled-test score (a theory is only as good as its worst failure), "
        "then median. '—' = no prediction (counts as 0). ‡ = too few numbers to hide: "
        "fitted and scored in-sample. "
        + (f"Averaged over {len(salts)} random public/hidden splits." if len(salts) > 1 else ""))
oq_note = ("Open questions: tests whose data are disputed or contradict each other, so "
           "even a correct theory could not pass them all. Reported, never scored, never "
           "used to fit constants.")
print("\n" + table + "\n\n" + note + "\n\n" + open_table + "\n\n" + oq_note)

if args.write:
    with open("LEADERBOARD.md", "w") as f:
        f.write("# UniverseBench leaderboard\n\n## Benchmark (settled tests)\n\n" + table
                + "\n\n" + note + "\n\n## Open questions\n\n" + open_table + "\n\n" + oq_note + "\n")
with open("results.json", "w") as f:
    json.dump([r["runs"] for r in rows], f, indent=1, default=str)
