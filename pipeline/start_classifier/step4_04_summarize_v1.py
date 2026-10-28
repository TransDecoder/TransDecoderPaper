#!/usr/bin/env python3
"""E33 step 4: choose rule variants on simulated human d1 + d8 only; report them held out.

usage: python step4_04_summarize_v1.py <variants_f1.tsv.gz> <out_dir>

Baseline = the a5b24c3 rules on the a5b24c3 inputs: variant 'k1|1e-05|best|top1|ga'.
Tuning score = mean change in exact F1 (points) over sim.human.d1 and sim.human.d8.
Selections (all on the tuning sets):
  best_overall    the variant with the largest tuning score (ties: fewer decisions changed on the tuning sets)
  best per family, varying one dimension from the baseline and holding the others at the a5b24c3 values:
    evalue x search (top1, best HSP, ga), depth/scheme (fast, 1e-5, best HSP, ga), hsp (fast, 1e-5, top1, ga),
    pfam (k1 hits); plus the best variant whose tuning score is the largest among variants with the
    a5b24c3 Pfam (ga) setting.
Held out: sim.human.d32, sim.athal.d1, the benchmark fragments and complete inputs (per input, median,
range, number down).
Writes tuning_ranking.tsv (all variants, tuning score), selected_variants.tsv (held-out deltas, overall and
per type), family_marginals.tsv (every family member on every dataset).
"""
import os
import sys

import numpy as np
import pandas as pd

vf, odir = sys.argv[1:3]
os.makedirs(odir, exist_ok=True)
V = pd.read_csv(vf, sep="\t")
BASE = "k1|1e-05|best|top1|ga"
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
mets = ["f1_exact"] + [f"f1_{t}" for t in TYPES]
b = V[V.variant == BASE].set_index("dataset")[mets]
for mt in mets:
    V[f"d_{mt}"] = 100 * (V[mt] - V.dataset.map(b[mt]))
TUNE = ["sim.human.d1", "sim.human.d8"]
cand = V[~V.variant.isin(["release", "perfect"])]
tun = cand[cand.dataset.isin(TUNE)].groupby("variant").agg(tune=("d_f1_exact", "mean"),
                                                           tune_changed=("n_changed_vs_release", "sum"))
for ds in TUNE:
    tun[f"d_{ds}"] = cand[cand.dataset == ds].set_index("variant").d_f1_exact
tun = tun.sort_values(["tune", "tune_changed"], ascending=[False, True])
tun.to_csv(f"{odir}/tuning_ranking.tsv", sep="\t", float_format="%.4f")

parts = tun.index.to_series().str.split("|", expand=True)
parts.columns = ["search", "evalue", "hsp", "depth", "pfam"]
sel = {"a5b24c3 (baseline)": BASE, "best_overall": tun.index[0]}
ga = tun[parts.pfam == "ga"]
sel["best_overall_ga"] = ga.index[0]
fam = {
    "evalue_x_search": (parts.hsp == "best") & (parts.depth == "top1") & (parts.pfam == "ga"),
    "depth_scheme": (parts.search == "fast") & (parts.evalue == "1e-05") & (parts.hsp == "best") & (parts.pfam == "ga"),
    "hsp": (parts.search == "fast") & (parts.evalue == "1e-05") & (parts.depth == "top1") & (parts.pfam == "ga"),
    "pfam": (parts.search == "k1"),
}
fam_rows = []
for f, mask in fam.items():
    members = tun[mask].index.tolist()
    sel[f"best_{f}"] = members[0]
    for v in members:
        for _, r in V[V.variant == v].iterrows():
            fam_rows.append(dict(family=f, variant=v, tune=tun.loc[v, "tune"], dataset=r.dataset,
                                 d_f1_exact=r.d_f1_exact))
pd.DataFrame(fam_rows).to_csv(f"{odir}/family_marginals.tsv", sep="\t", index=False, float_format="%.4f")

groups = {
    "sim.human.d1": lambda d: d == "sim.human.d1",
    "sim.human.d8": lambda d: d == "sim.human.d8",
    "sim.human.d32": lambda d: d == "sim.human.d32",
    "sim.athal.d1": lambda d: d == "sim.athal.d1",
    "bench.fragments": lambda d: d.startswith("bench.") and d.endswith(".fragments"),
    "bench.complete": lambda d: d.startswith("bench.") and d.endswith(".complete"),
}
rows = []
for name, v in list(sel.items()) + [("perfect (oracle on all candidates)", "perfect"), ("release (ab initio)", "release")]:
    x = V[V.variant == v]
    for g, fn in groups.items():
        y = x[x.dataset.map(fn)]
        if not len(y):
            continue
        r = dict(selection=name, variant=v, group=g, n_inputs=len(y), tune=tun.tune.get(v, np.nan))
        for mt in mets:
            d = y[f"d_{mt}"].dropna()
            r[f"median_d_{mt}"] = d.median() if len(d) else np.nan
            r[f"min_d_{mt}"] = d.min() if len(d) else np.nan
            r[f"max_d_{mt}"] = d.max() if len(d) else np.nan
        r["n_down_f1_exact"] = int((y.d_f1_exact < -1e-9).sum())
        r["n_up_f1_exact"] = int((y.d_f1_exact > 1e-9).sum())
        r["sum_d_tp_exact"] = int(y.tp_exact.sum() - V[(V.variant == BASE) & V.dataset.isin(y.dataset)].tp_exact.sum())
        rows.append(r)
S = pd.DataFrame(rows)
S.to_csv(f"{odir}/selected_variants.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
print(S[["selection", "variant", "group", "n_inputs", "median_d_f1_exact", "min_d_f1_exact", "max_d_f1_exact",
         "n_down_f1_exact", "sum_d_tp_exact"]].to_string(float_format=lambda z: f"{z:.3f}"))
print(tun.head(15).to_string(float_format=lambda z: f"{z:.4f}"))
