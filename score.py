"""
Scorecard: every model on every test, scored on hidden data.

    BENCH_SALT='my-private-secret' python score.py              # baselines
    BENCH_SALT='my-private-secret' python score.py my_model.py  # + a candidate
    python score.py --splits 5 --write                          # baselines averaged
                                                                # over 5 public splits,
                                                                # written to LEADERBOARD.md

A candidate file must define MODEL = <instance with the models.py interface>.
"""
import argparse
import importlib.util
import json

import numpy as np

import benchmark as B
import clusters
import milkyway
import models as M

ap = argparse.ArgumentParser()
ap.add_argument("candidates", nargs="*")
ap.add_argument("--splits", type=int, default=1,
                help="average over N public demo splits (ignores BENCH_SALT)")
ap.add_argument("--write", action="store_true", help="write LEADERBOARD.md")
args = ap.parse_args()

candidates = list(M.BASELINES)
for path in args.candidates:
    spec = importlib.util.spec_from_file_location("candidate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    candidates.append(mod.MODEL)

salts = [B.SPLIT_SALT] if args.splits == 1 else [f"public-split-{i}" for i in range(args.splits)]
gals, cl, mw = B.load_all(), clusters.load_clusters(), milkyway.MilkyWay()

rows = []
for m in candidates:
    runs = [B.run(m, gals, cl, mw, salt=s, verbose=False) for s in salts]
    row = {"model": m.NAME, "knobs": f"{runs[0]['n_global_params']} global, "
                                      f"{runs[0]['n_local_params']} per galaxy"}
    for t in B.TESTS:
        r = runs[0][t + "_hidden"]
        row[t] = None if r.get("no_prediction") else \
            float(np.mean([x[t + "_hidden"]["within_2sigma_pct"] for x in runs]))
    vals = [row[t] or 0.0 for t in B.TESTS]
    row["LOWEST"], row["MEDIAN"], row["HIGHEST"] = min(vals), float(np.median(vals)), max(vals)
    row["runs"] = runs
    rows.append(row)
# rank by weakest test, then median
rows.sort(key=lambda r: (-r["LOWEST"], -r["MEDIAN"]))

fmt = lambda v: "—" if v is None else f"{v:.0f}"
head = ["#", "Model"] + [B.TEST_LABELS[t] for t in B.TESTS] + \
       ["Lowest", "Median", "Highest", "Knobs"]
lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
for i, r in enumerate(rows, 1):
    cells = [str(i), r["model"]] + [fmt(r[t]) for t in B.TESTS] + \
            [f"**{r['LOWEST']:.0f}**", f"{r['MEDIAN']:.0f}", f"{r['HIGHEST']:.0f}", r["knobs"]]
    lines.append("| " + " | ".join(cells) + " |")
table = "\n".join(lines)

note = (f"Scores: % of hidden data points predicted within 2 sigma (a perfect theory "
        f"scores about 95). Lowest / median / highest are taken across all tests; models "
        f"are ranked by their lowest score (a theory is only as good as its worst failure), "
        f"then by median. '—' = model makes no prediction for that test (counts as 0). "
        f"{'Averaged over ' + str(len(salts)) + ' random public/hidden splits.' if len(salts) > 1 else ''}")
print(table + "\n\n" + note)

if args.write:
    with open("LEADERBOARD.md", "w") as f:
        f.write("# UniverseBench leaderboard\n\n" + table + "\n\n" + note + "\n")
with open("results.json", "w") as f:
    json.dump([r["runs"] for r in rows], f, indent=1, default=str)
