#!/usr/bin/env python3
"""E33 adoption: score the homology-start-assist runs (PyTD a5b24c3+) against the releases.

usage: python 16_score_adoption.py <out_dir>

Benchmark (84 homology inputs): tools pytd_v7 (ab initio), pytd_v7_homology (3a793ff, homology for
selection only) and pytd_v7h_homology (new code, same hits) in $P/work/bench_v3/runs/<org>/<variant>/.
Simulated assemblies human d1, d8, d32, athal d1: pytd_v7 (ab initio release run), and the new code
with the E33 step-3 candidate searches, assist off (`off`, = homology for selection only) and on (`on`)
in $P/work/start_classifier/hstart_runs/<dataset>/{off,on}; predictions on excluded contigs removed.
Scoring as score_predictions.py (stop within 2 nt, exact = start also within 2 nt, one credit per
truth CDS; F1 = 2 tp / (n_truth + n_pred)); per truth type F1 = 2 tp_exact / (n_truth of the type +
predictions on transcripts of that type).
Writes bench_scores.tsv, sim_scores.tsv and bench_summary.tsv (median and range of the differences).
"""
import collections
import os
import sys

import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
sys.path.insert(0, f"{P}/TransDecoderPaper/benchmark_v2/scripts")
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]


def score(truth, preds):
    n_pred = sum(len(v) for v in preds.values())
    row = dict(n_truth=len(truth), n_pred=n_pred, tp_stop=0, tp_exact=0)
    per = {t: [0, 0, 0] for t in TYPES}
    for tid, t in truth.items():
        best = None
        for p in preds.get(tid, []):
            if p[2] != t[2]:
                continue
            p5, p3 = ends(*p)
            t5, t3 = ends(*t[:3])
            if abs(p3 - t3) <= TOL:
                lvl = "exact" if abs(p5 - t5) <= TOL else "stop"
                if best != "exact":
                    best = lvl
        row["tp_stop"] += best is not None
        row["tp_exact"] += best == "exact"
        if t[3] in per:
            a = per[t[3]]
            a[0] += 1
            a[1] += len(preds.get(tid, []))
            a[2] += best == "exact"
    row["f1_stop"] = 2 * row["tp_stop"] / (row["n_truth"] + n_pred)
    row["f1_exact"] = 2 * row["tp_exact"] / (row["n_truth"] + n_pred)
    for t, (n, npr, tp) in per.items():
        row[f"n_{t}"], row[f"tp_exact_{t}"] = n, tp
        row[f"f1_exact_{t}"] = 2 * tp / (n + npr) if n + npr else float("nan")
    return row


rows = []
R = f"{P}/work/bench_v3/runs"
for org in sorted(os.listdir(R)):
    if not os.path.isdir(f"{R}/{org}") or org == "logs":
        continue
    for v in sorted(os.listdir(f"{R}/{org}")):
        if not os.path.exists(f"{R}/{org}/{v}/pytd_v7_homology/preds.tsv"):
            continue
        truth = read_truth(f"{P}/data/bench_v3/inputs/{org}/{v}.truth.tsv")
        for tool in ("pytd_v7", "pytd_v7_homology", "pytd_v7h_homology"):
            pf = f"{R}/{org}/{v}/{tool}/preds.tsv"
            if os.path.exists(pf):
                rows.append(dict(org=org, variant=v, tool=tool, **score(truth, read_preds(pf))))
b = pd.DataFrame(rows)
b.to_csv(f"{out}/bench_scores.tsv", sep="\t", index=False, float_format="%.5f")

# differences per input, summarised
w = b.pivot_table(index=["org", "variant"], columns="tool",
                  values=["f1_stop", "f1_exact"] + [f"f1_exact_{t}" for t in TYPES])
summ = []
for ref in ("pytd_v7_homology", "pytd_v7"):
    for m in ["f1_stop", "f1_exact"] + [f"f1_exact_{t}" for t in TYPES]:
        if (m, "pytd_v7h_homology") not in w.columns or (m, ref) not in w.columns:
            continue
        d = 100 * (w[(m, "pytd_v7h_homology")] - w[(m, ref)]).dropna()
        for grp, dd in [("all", d)] + [(g, d[d.index.get_level_values(1).str.split(".").str[0] == g])
                                       for g in ("complete", "complete_fullutr", "fragments")]:
            if len(dd):
                summ.append(dict(vs=ref, metric=m, variant_group=grp, n_inputs=len(dd), median=dd.median(),
                                 min=dd.min(), max=dd.max(), n_down=int((dd < -1e-9).sum())))
pd.DataFrame(summ).to_csv(f"{out}/bench_summary.tsv", sep="\t", index=False, float_format="%.3f")

# simulated assemblies
man = pd.read_csv(f"{P}/TransDecoderPaper/data/start_classifier/step2_manifest.tsv", sep="\t").set_index("dataset")
srows = []
for ds in ("sim.human.d1", "sim.human.d8", "sim.human.d32", "sim.athal.d1"):
    m = man.loc[ds]
    excl = {l.split("\t")[0] for l in open(m.excluded)}
    truth = read_truth(m.truth)
    for arm, pf in (("pytd_v7", f"{m.run_dir}/preds.tsv"),
                    ("assist_off", f"{P}/work/start_classifier/hstart_runs/{ds}/off/preds.tsv"),
                    ("assist_on", f"{P}/work/start_classifier/hstart_runs/{ds}/on/preds.tsv")):
        if not os.path.exists(pf):
            continue
        preds = {k: v for k, v in read_preds(pf).items() if k not in excl}
        srows.append(dict(dataset=ds, arm=arm, **score(truth, preds)))
pd.DataFrame(srows).to_csv(f"{out}/sim_scores.tsv", sep="\t", index=False, float_format="%.5f")
print("done", len(b), "bench rows,", len(srows), "sim rows")
