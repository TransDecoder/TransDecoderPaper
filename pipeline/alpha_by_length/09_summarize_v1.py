#!/usr/bin/env python3
"""E36: summary tables from 07_evaluate_v1.py output.

usage: python 09_summarize_v1.py <eval_dir> <out_dir>

Writes (out_dir):
  scores.tsv              dataset x config: exact F1, per type f1 / f1w, deltas (points) vs V0 (copy of the eval table,
                          main columns only)
  bench_summary.tsv       per arm x config: E13 objective (8 complete + 8 fragments, equal weight; abs and delta),
                          mean deltas complete / complete_fullutr / fragments / mixed, inputs < -0.5, worst input,
                          E35 acceptance, both_up (complete and fragmented inputs both up on average)
  heldout_summary.tsv     per arm x config: delta exact F1 on each simulated assembly, mean / min over human and athal,
                          GM12878 E29 statistic deltas
  delta_by_type.tsv       simulated assemblies: per type delta f1w (with negatives)
  changes_summary.tsv     per arm x config x dataset group: changed decisions by class (w2r, r2w, w2w, masked, nostop,
                          no_truth, no_truth_neg) and net exact
"""
import sys

import numpy as np
import pandas as pd

edir, odir = sys.argv[1:3]
S = pd.read_csv(f"{edir}/scores.tsv", sep="\t")
C = pd.read_csv(f"{edir}/changes.tsv", sep="\t")
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
keep = ["dataset", "arm", "group", "org", "input", "config", "form", "mode", "tuned_arm", "params", "n_changed",
        "f1_exact", "d_f1_exact"] + [f"d_f1w_{t}" for t in TYPES] + [c for c in S.columns if c.startswith("e29_")]
S[keep].to_csv(f"{odir}/scores.tsv", sep="\t", index=False, float_format="%.5f")

# benchmark
B = S[S.group.str.startswith("bench")].copy()
B["cls"] = B.input.str.replace(r"_t[123]$", "", regex=True)
rows = []
for (arm, cfg), g in B.groupby(["arm", "config"], sort=False):
    t = g.pivot_table(index="org", columns="cls", values="d_f1_exact")
    a = g.pivot_table(index="org", columns="cls", values="f1_exact")
    w = g.loc[g.d_f1_exact.idxmin()]
    r = dict(arm=arm, config=cfg, params=g.params.iloc[0],
             E13_abs=float((0.5 * (a.complete + a.fragments)).mean()), E13_d=float((0.5 * (t.complete + t.fragments)).mean()),
             complete=t.complete.mean(), complete_fullutr=t.complete_fullutr.mean(), fragments=t.fragments.mean(),
             mixed=t.fragments_mixed.mean() if "fragments_mixed" in t else np.nan, n_inputs=len(g),
             n_below_05=int((g.d_f1_exact < -0.5).sum()), worst=f"{w.dataset} {w.d_f1_exact:+.2f}",
             n_complete_down=int((t.complete < 0).sum()), n_fragments_down=int((t.fragments < 0).sum()))
    r["acceptable"] = bool(r["complete"] >= -0.1 and r["n_below_05"] == 0)
    r["both_up"] = bool(r["complete"] > 0 and r["fragments"] > 0 and (np.isnan(r["mixed"]) or r["mixed"] > 0))
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{odir}/bench_summary.tsv", sep="\t", index=False, float_format="%.3f")

# held out
H = S[~S.group.str.startswith("bench")].copy()
H["col"] = H.dataset.str.replace("sim.", "", regex=False)
rows = []
for (arm, cfg), g in H.groupby(["arm", "config"], sort=False):
    x = g.set_index("col")
    r = dict(arm=arm, config=cfg, params=g.params.iloc[0])
    for c in x.index:
        if not c.startswith("gm12878"):
            r[c.replace(".H", "")] = x.loc[c, "d_f1_exact"]
    hu = g[g.group == "sim_human"].d_f1_exact
    at = g[g.group == "sim_athal"].d_f1_exact
    r.update(human_mean=hu.mean(), human_min=hu.min(), athal_mean=at.mean(), athal_min=at.min(),
             n_sims_down=int((g[g.group != "gm12878"].d_f1_exact < -0.005).sum()))
    gm = g[g.group == "gm12878"]
    if len(gm):
        v0 = S[(S.config == "V0") & (S.dataset == gm.dataset.iloc[0])].iloc[0]
        for t in ("5prime_partial", "complete"):
            r[f"gm_e29_{t}"] = 100 * gm[f"e29_exact_given_stop_{t}"].iloc[0]
            r[f"gm_e29_{t}_d"] = 100 * (gm[f"e29_exact_given_stop_{t}"].iloc[0] - v0[f"e29_exact_given_stop_{t}"])
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{odir}/heldout_summary.tsv", sep="\t", index=False, float_format="%.3f")

T = H[H.group != "gm12878"][["dataset", "arm", "config", "d_f1_exact"] + [f"d_f1w_{t}" for t in TYPES]]
T.to_csv(f"{odir}/delta_by_type.tsv", sep="\t", index=False, float_format="%.3f")

C["grp"] = np.where(C.group.str.startswith("bench"), C.group, C.group)
P = C.pivot_table(index=["arm", "config", "grp"], columns="cls", values="n", aggfunc="sum", fill_value=0)
for c in ["w2r", "r2w", "w2w", "masked", "nostop", "no_truth", "no_truth_neg"]:
    if c not in P:
        P[c] = 0
P["changed"] = P.sum(axis=1)
P["net_exact"] = P.w2r - P.r2w
P = P[["changed", "w2r", "r2w", "net_exact", "w2w", "masked", "nostop", "no_truth", "no_truth_neg"]].reset_index()
P.to_csv(f"{odir}/changes_summary.tsv", sep="\t", index=False)
