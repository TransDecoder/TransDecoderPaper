#!/usr/bin/env python3
"""E33 step 5: summarise UniRef50 vs Swiss-Prot (start-decision isolate); choose variants on human d1 + d8 only.

usage: python step5_07_summarize_v1.py <variants_f1.tsv.gz> <hit_stats.tsv> <out_dir>

Baseline = the a5b24c3 rules on Swiss-Prot -k 1 hits: 'sprot|k1|top1|none'.
Tuning score = mean change in exact F1 (points) vs the baseline over sim.human.d1 and sim.human.d8.
Selections (on the tuning sets only), each reported held out on human d32, athal d1 and d16, the 8 benchmark
fragments and the 8 complete inputs:
  uniref50 a5b24c3        'uniref50|k1|top1|none' (fixed, not selected)
  best_uniref50           best UniRef50 variant overall
  best_uniref50_<family>  best UniRef50 -k 25 variant within: topn (R1/R2 top-n schemes), cons (consensus
                          N-terminus), with each M filter (none / r2 / all)
  best_sprot              best Swiss-Prot variant overall (for reference: is any gain database-specific?)
Writes: main_by_dataset.tsv (release, both a5b24c3 arms, perfect: F1 overall and per type, rule firings,
wrong trims / keeps), tuning_ranking.tsv, selected_variants.tsv (held-out summaries), selected_by_dataset.tsv,
hit_stats.tsv (copy, ordered).
"""
import os
import sys

import numpy as np
import pandas as pd

vf, hf, odir = sys.argv[1:4]
os.makedirs(odir, exist_ok=True)
V = pd.read_csv(vf, sep="\t")
BASE = "sprot|k1|top1|none"
UBASE = "uniref50|k1|top1|none"
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
mets = ["f1_exact"] + [f"f1_{t}" for t in TYPES]
b = V[V.variant == BASE].set_index("dataset")
for mt in mets:
    V[f"d_{mt}"] = 100 * (V[mt] - V.dataset.map(b[mt]))
V["d_tp_exact"] = V.tp_exact - V.dataset.map(b.tp_exact)
order = list(dict.fromkeys(V.dataset))

# main comparison per dataset
keep_cols = ["dataset", "variant", "f1_exact"] + [f"f1_{t}" for t in TYPES] + [
    "tp_exact", "d_f1_exact", "n_changed_vs_release", "n_R1", "n_R2", "n_r3", "n_conflict", "n_blast_trim",
    "n_blast_keep", "wrong_trim", "wrong_keep", "right_trim", "right_keep"]
M = V[V.variant.isin(["release", BASE, UBASE, "perfect"])][keep_cols].copy()
M["variant"] = pd.Categorical(M.variant, ["release", BASE, UBASE, "perfect"])
M["dataset"] = pd.Categorical(M.dataset, order)
M.sort_values(["dataset", "variant"]).to_csv(f"{odir}/main_by_dataset.tsv", sep="\t", index=False,
                                              float_format="%.5f")

TUNE = ["sim.human.d1", "sim.human.d8"]
cand = V[~V.variant.isin(["release", "perfect"])]
tun = cand[cand.dataset.isin(TUNE)].groupby("variant").agg(tune=("d_f1_exact", "mean"),
                                                           tune_changed=("n_changed_vs_release", "sum"))
for ds in TUNE:
    tun[f"d_{ds}"] = cand[cand.dataset == ds].set_index("variant").d_f1_exact
tun = tun.sort_values(["tune", "tune_changed"], ascending=[False, True])
parts = tun.index.to_series().str.split("|", expand=True)
parts.columns = ["db", "search", "depth", "mfilt"]
tun = tun.join(parts)
tun.to_csv(f"{odir}/tuning_ranking.tsv", sep="\t", float_format="%.4f")

U = tun[tun.db == "uniref50"]
U25 = U[U.search == "k25"]
sel = {"sprot a5b24c3 (baseline)": BASE, "uniref50 a5b24c3 rules (k1 top1)": UBASE,
       "uniref50 k25 top1": "uniref50|k25|top1|none",
       "uniref50 k1 top1, R2 needs subject M": "uniref50|k1|top1|r2",
       "uniref50 k1 top1, non-M subjects dropped": "uniref50|k1|top1|all",
       "best_uniref50": U.index[0],
       "best_uniref50_topn": U25[(U25.depth.str.startswith("top")) & (U25.depth != "top1") & (U25.mfilt == "none")].index[0],
       "best_uniref50_cons": U25[U25.depth.str.startswith("cons") & (U25.mfilt == "none")].index[0]}
for mf in ("none", "r2", "all"):
    sel[f"best_uniref50_mfilt_{mf}"] = U[U.mfilt == mf].index[0]
sel["best_sprot"] = tun[tun.db == "sprot"].index[0]
sel["best_sprot_cons"] = tun[(tun.db == "sprot") & tun.depth.str.startswith("cons")].index[0]

groups = {
    "sim.human.d1 (tune)": lambda d: d == "sim.human.d1",
    "sim.human.d8 (tune)": lambda d: d == "sim.human.d8",
    "sim.human.d32": lambda d: d == "sim.human.d32",
    "sim.athal.d1": lambda d: d == "sim.athal.d1",
    "sim.athal.d16": lambda d: d == "sim.athal.d16",
    "bench.fragments": lambda d: d.startswith("bench.") and d.endswith(".fragments"),
    "bench.complete": lambda d: d.startswith("bench.") and d.endswith(".complete"),
}
rows, byds = [], []
for name, v in list(sel.items()) + [("perfect (oracle)", "perfect"), ("release (ab initio)", "release")]:
    x = V[V.variant == v]
    for _, r in x.iterrows():
        byds.append(dict(selection=name, variant=v, dataset=r.dataset, **{f"d_{m}": r[f"d_{m}"] for m in mets},
                         d_tp_exact=r.d_tp_exact, f1_exact=r.f1_exact, wrong_trim=r.wrong_trim, wrong_keep=r.wrong_keep,
                         n_changed_vs_release=r.n_changed_vs_release))
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
        r["sum_d_tp_exact"] = int(y.d_tp_exact.sum())
        r["sum_wrong_trim"] = y.wrong_trim.sum()
        r["sum_wrong_keep"] = y.wrong_keep.sum()
        rows.append(r)
S = pd.DataFrame(rows)
S.to_csv(f"{odir}/selected_variants.tsv", sep="\t", index=False, float_format="%.4f")
pd.DataFrame(byds).to_csv(f"{odir}/selected_by_dataset.tsv", sep="\t", index=False, float_format="%.4f")
H = pd.read_csv(hf, sep="\t")
H["dataset"] = pd.Categorical(H.dataset, order)
H.sort_values(["dataset", "db", "search"]).to_csv(f"{odir}/hit_stats.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
print(S[["selection", "variant", "group", "n_inputs", "median_d_f1_exact", "min_d_f1_exact", "max_d_f1_exact",
         "n_down_f1_exact", "sum_d_tp_exact", "sum_wrong_trim", "sum_wrong_keep"]].to_string(float_format=lambda z: f"{z:.3f}"))
print(tun.head(20).to_string(float_format=lambda z: f"{z:.4f}"))
