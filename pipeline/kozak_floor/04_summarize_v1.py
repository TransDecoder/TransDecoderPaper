#!/usr/bin/env python3
"""E35: summary tables of the Kozak-floor variants (deltas vs V0 in F1 points = 100 x difference).

usage: python 04_summarize_v1.py <eval_dir> <out_dir>

Inputs: <eval_dir>/{scores,changes,checks}.tsv from 02_evaluate_v1.py.
Writes to <out_dir>:
  scores.tsv              copy of the per-dataset x variant scores (absolute values)
  delta_overall.tsv       dataset x variant: exact F1, d exact F1, d stop F1, decisions changed
  delta_by_type.tsv       dataset x variant x truth type: d f1 (per type) and d f1w (per type with negatives)
  group_summary.tsv       arm x group x variant: mean / min / max d exact F1, inputs down, inputs down > 0.5
  e13_objective.tsv       arm x variant: mean exact F1 of the 8 complete and the 8 fragments benchmark inputs and
                          their equal-weighted mean (E13 objective), with deltas
  acceptance.tsv          arm x variant: pre-set criteria (mean d on complete inputs >= -0.1; no benchmark input
                          < -0.5), worst input, verdict
  changes_summary.tsv     arm x group x variant: changed decisions by class (w2r, r2w, w2w_later, w2w_upstream,
                          w2w_at_atg, nostop_truth, no_truth, no_truth_neg)
  changes_by_type.tsv     same, split by truth type
  kozak_veto.tsv          dataset: candidates, kept open, Kozak-veto candidates (V0)
"""
import os
import sys

import numpy as np
import pandas as pd

edir, odir = sys.argv[1:3]
os.makedirs(odir, exist_ok=True)
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
S = pd.read_csv(f"{edir}/scores.tsv", sep="\t")
Ch = pd.read_csv(f"{edir}/changes.tsv", sep="\t")
K = pd.read_csv(f"{edir}/checks.tsv", sep="\t")
VARS = list(dict.fromkeys(S.variant))
S.to_csv(f"{odir}/scores.tsv", sep="\t", index=False, float_format="%.6f")

base = S[S.variant == "V0"].set_index("dataset")
S = S.join(base[["f1_exact", "f1_stop"] + [f"f1_{t}" for t in TYPES] + [f"f1w_{t}" for t in TYPES]],
           on="dataset", rsuffix="_v0")
S["d_exact"] = 100 * (S.f1_exact - S.f1_exact_v0)
S["d_stop"] = 100 * (S.f1_stop - S.f1_stop_v0)
S[["dataset", "arm", "group", "org", "input", "variant", "f1_exact", "d_exact", "d_stop", "n_changed"]].to_csv(
    f"{odir}/delta_overall.tsv", sep="\t", index=False, float_format="%.4f")

rows = []
for _, r in S.iterrows():
    for t in TYPES:
        if r[f"n_{t}"] == 0:
            continue
        rows.append(dict(dataset=r.dataset, arm=r.arm, group=r.group, variant=r.variant, type=t, n=r[f"n_{t}"],
                         f1=r[f"f1_{t}"], d_f1=100 * (r[f"f1_{t}"] - r[f"f1_{t}_v0"]),
                         f1w=r[f"f1w_{t}"], d_f1w=100 * (r[f"f1w_{t}"] - r[f"f1w_{t}_v0"])))
pd.DataFrame(rows).to_csv(f"{odir}/delta_by_type.tsv", sep="\t", index=False, float_format="%.4f")

g = S.groupby(["arm", "group", "variant"], sort=False).agg(
    n_inputs=("d_exact", "size"), mean_d=("d_exact", "mean"), min_d=("d_exact", "min"), max_d=("d_exact", "max"),
    n_down=("d_exact", lambda x: int((x < -1e-9).sum())), n_down_gt_0_5=("d_exact", lambda x: int((x < -0.5).sum())),
    changed=("n_changed", "sum")).reset_index()
g.to_csv(f"{odir}/group_summary.tsv", sep="\t", index=False, float_format="%.4f")

e13, acc = [], []
for arm in ("abinitio", "homology"):
    A = S[S.arm == arm]
    for v in VARS:
        X = A[A.variant == v]
        c = X[X.group == "bench_complete"]
        cf = X[X.group.isin(["bench_complete", "bench_complete_fullutr"])]
        f = X[X.group == "bench_fragments"]
        e13.append(dict(arm=arm, variant=v, mean_complete=c.f1_exact.mean(), mean_fragments=f.f1_exact.mean(),
                        objective=(c.f1_exact.mean() + f.f1_exact.mean()) / 2,
                        d_complete=c.d_exact.mean(), d_fragments=f.d_exact.mean(),
                        d_objective=(c.d_exact.mean() + f.d_exact.mean()) / 2))
        B = X[X.group.str.startswith("bench")]
        w = B.loc[B.d_exact.idxmin()]
        crit1 = cf.d_exact.mean() >= -0.1
        crit2 = B.d_exact.min() >= -0.5
        acc.append(dict(arm=arm, variant=v, mean_d_complete=c.d_exact.mean(), mean_d_complete_and_fullutr=cf.d_exact.mean(),
                        n_bench_inputs=len(B), n_bench_down_gt_0_5=int((B.d_exact < -0.5).sum()),
                        worst_input=w.dataset, worst_d=w.d_exact,
                        crit_complete_ok=int(crit1), crit_no_input_below_0_5_ok=int(crit2),
                        e13_d_objective=e13[-1]["d_objective"], acceptable=int(crit1 and crit2)))
pd.DataFrame(e13).to_csv(f"{odir}/e13_objective.tsv", sep="\t", index=False, float_format="%.5f")
pd.DataFrame(acc).to_csv(f"{odir}/acceptance.tsv", sep="\t", index=False, float_format="%.4f")

CLS = ["w2r", "r2w", "w2w_later", "w2w_upstream", "w2w_at_atg", "nostop_truth", "no_truth", "no_truth_neg"]
Ch = Ch[Ch.variant != "V0"]
cs = Ch.pivot_table(index=["arm", "group", "variant"], columns="cls", values="n", aggfunc="sum", fill_value=0)
cs = cs.reindex(columns=CLS, fill_value=0)
cs.insert(0, "changed", cs.sum(axis=1))
cs["net_exact"] = cs.w2r - cs.r2w
cs.reset_index().to_csv(f"{odir}/changes_summary.tsv", sep="\t", index=False)
ct = Ch.pivot_table(index=["arm", "group", "variant", "truth_type"], columns="cls", values="n", aggfunc="sum",
                    fill_value=0).reindex(columns=CLS, fill_value=0)
ct.insert(0, "changed", ct.sum(axis=1))
ct.reset_index().to_csv(f"{odir}/changes_by_type.tsv", sep="\t", index=False)
K.to_csv(f"{odir}/kozak_veto.tsv", sep="\t", index=False)
print("done")
