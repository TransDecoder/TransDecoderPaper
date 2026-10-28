#!/usr/bin/env python3
"""E33 step 1: which features carry the separation? Logistic-regression feature-set ablation.

usage: python 03_ablation.py <features_dir> <out.tsv>

Same candidates and labels as 02_separation.py. For each feature set, AUC and accuracy (p >= 0.5)
of a standardised logistic regression under: 5-fold GroupKFold on the pooled benchmark (cv_bench),
5-fold CV within each simulated assembly (cv_sim), and benchmark -> simulated assembly (bench2sim).
"""
import glob
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

fdir, out = sys.argv[1], sys.argv[2]
df = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in sorted(glob.glob(f"{fdir}/*.tsv.gz"))],
               ignore_index=True)
d = df[(df.label == "truncated") | ((df.label == "utr") & (df.true_atg_rank >= 1))].copy()
d["y"] = (d.label == "truncated").astype(int)
d["source"] = d.dataset.str.split(".").str[0]
d["org"] = d.dataset.str.split(".").str[1]
d["log_seg"] = np.log(d.seg)
d["coding_z"] = d.coding / np.sqrt(d.seg)
d["log_tx"] = np.log(d.tx_len)
d["log_orf"] = np.log(d.orf_len)
d["log_utr3"] = np.log1p(d.utr3)
d["log_orf_rel"] = np.log(d.orf_len_rel)
d["grp"] = d.org + ":" + d.transcript_id.astype(str)

RULE = ["seg", "coding", "kozak"]
SETS = {
    "rule3": RULE,
    "rule3+log_seg": RULE + ["log_seg"],
    "rule3+coding_z": RULE + ["coding_z"],
    "rule3+coding_z+log_seg": RULE + ["coding_z", "log_seg"],
    "coding_z+kozak": ["coding_z", "kozak"],
    "rule3+coding_z+log_seg+seg_gc": RULE + ["coding_z", "log_seg", "seg_gc"],
    "intr": RULE + ["coding_z", "log_seg", "coding_pb", "kozak_off", "seg_gc"],
}
INTR = SETS["intr"]
for f in ["log_tx", "log_orf", "log_orf_rel", "internal", "log_utr3"]:
    SETS[f"intr+{f}"] = INTR + [f]
SETS["full"] = INTR + ["log_tx", "log_orf", "log_orf_rel", "internal", "log_utr3"]


def lr():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))


rows = []
bench = d[d.source == "bench"].reset_index(drop=True)
sims = sorted(d[d.source == "sim"].dataset.unique())
for name, feats in SETS.items():
    p = np.zeros(len(bench))
    for tr, te in GroupKFold(5).split(bench, bench.y, bench.grp):
        p[te] = lr().fit(bench.loc[tr, feats], bench.y[tr]).predict_proba(bench.loc[te, feats])[:, 1]
    rows.append(dict(features=name, scheme="cv_bench", test="ALL", auc=roc_auc_score(bench.y, p),
                     acc=((p >= .5) == bench.y).mean()))
    m = lr().fit(bench[feats], bench.y)
    for ds in sims:
        s = d[d.dataset == ds].reset_index(drop=True)
        q = m.predict_proba(s[feats])[:, 1]
        rows.append(dict(features=name, scheme="bench2sim", test=ds, auc=roc_auc_score(s.y, q), acc=((q >= .5) == s.y).mean()))
        q = np.zeros(len(s))
        for tr, te in GroupKFold(5).split(s, s.y, s.transcript_id):
            q[te] = lr().fit(s.loc[tr, feats], s.y[tr]).predict_proba(s.loc[te, feats])[:, 1]
        rows.append(dict(features=name, scheme="cv_sim", test=ds, auc=roc_auc_score(s.y, q), acc=((q >= .5) == s.y).mean()))
r = pd.DataFrame(rows)
r.to_csv(out, sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
for m in ["auc", "acc"]:
    print(m)
    print(r[r.test.isin(["ALL", "sim.human.d1", "sim.human.d32", "sim.athal.d1"])]
          .pivot_table(index="features", columns=["scheme", "test"], values=m).round(3))
