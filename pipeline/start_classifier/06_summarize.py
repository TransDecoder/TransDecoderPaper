#!/usr/bin/env python3
"""E33 step 2: compact summary of 05_evaluate.py output.

usage: python 06_summarize.py <step2_dir>
Writes <step2_dir>/summary_by_group.tsv (per input group x method: mean decision accuracy, AUC, and
the change in exact F1 points vs the release overall and per type; min/max over inputs; inputs where
overall exact F1 fell) and <step2_dir>/summary_sim.tsv (each simulated assembly x method).
"""
import sys

import numpy as np
import pandas as pd

d = sys.argv[1]
dm = pd.read_csv(f"{d}/decision_metrics.tsv", sep="\t")
f1 = pd.read_csv(f"{d}/exact_f1.tsv", sep="\t")
x = f1.merge(dm, on=["dataset", "method"])
rel = x[x.method == "release"].set_index("dataset")
cols = ["f1_exact", "f1_complete", "f1_5prime_partial", "f1_3prime_partial", "f1_internal"]
for c in cols:
    x["d_" + c] = 100 * (x[c] - x.dataset.map(rel[c]))
x["d_acc"] = x.acc - x.dataset.map(rel.acc)


def group(ds):
    s, org, inp = ds.split(".")
    if s == "sim":
        return f"sim_{org}"
    return inp.split("_t")[0] if inp.startswith("fragments_mixed") else inp


x["group"] = x.dataset.map(group)
rows = []
for (g, m), h in x.groupby(["group", "method"]):
    r = dict(group=g, method=m, n_inputs=len(h), acc=h.acc.mean(), auc=h.auc.mean(), d_acc=h.d_acc.mean(),
             d_f1_min=h.d_f1_exact.min(), d_f1_max=h.d_f1_exact.max(), n_f1_down=(h.d_f1_exact < -1e-9).sum())
    for c in cols:
        r["d_" + c] = h["d_" + c].mean()
    rows.append(r)
s = pd.DataFrame(rows)
s.to_csv(f"{d}/summary_by_group.tsv", sep="\t", index=False, float_format="%.4f")
sim = x[x.group.str.startswith("sim")][["dataset", "method", "acc", "auc", "n_changed", "tp_exact", "f1_exact"] +
                                       ["d_" + c for c in cols]]
sim.to_csv(f"{d}/summary_sim.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 400); pd.set_option("display.max_columns", 30)
print(s.round(3).to_string(index=False))
print(sim.round(4).to_string(index=False))
