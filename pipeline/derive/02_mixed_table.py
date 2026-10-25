#!/usr/bin/env python3
"""Medians for the balanced mixed inputs (E28), as used in the manuscript.

For each method, measure and ORF type (or the whole input): median over organisms within each
trial, then the mean over the three trials (the convention of notebook 02, "Balanced mixed
inputs"). Tools without the ciliate code are not applicable to T. thermophila (7 organisms).
F1 per type is computed from that type's sensitivity and precision.

Usage: 02_mixed_table.py  (run from the paper repo root)
Writes data/benchmark/final_v7/mixed_medians.tsv
"""
import pandas as pd

D = "data/benchmark/final_v7/"
NA = {"gmst", "codan", "cpc2", "rnasamba", "borf"}
MEAS = ["sensitivity_stop", "precision_stop", "f1_stop",
        "sensitivity_exact", "precision_exact", "f1_exact"]


def f1(s, p):
    return (2 * s * p / (s + p)).where(s + p > 0, 0.0)


def prep(d):
    d = d[d.variant.str.startswith("fragments_mixed")]
    return d[~((d.org == "tthermophila") & d.tool.isin(NA))].copy()


whole = prep(pd.read_csv(D + "mixed_scores.tsv", sep="\t"))
whole["type"] = "all"
types = prep(pd.read_csv(D + "partial_types.tsv", sep="\t"))
for lvl in ("stop", "exact"):
    types[f"f1_{lvl}"] = f1(types[f"sensitivity_{lvl}"], types[f"precision_{lvl}"])

d = pd.concat([whole, types], ignore_index=True)
per_trial = d.groupby(["tool", "type", "variant"])[MEAS].median()
out = per_trial.groupby(["tool", "type"])[MEAS].mean()
rng = per_trial.groupby(["tool", "type"])["f1_exact"].agg(lambda x: x.max() - x.min())
out["f1_exact_trial_range"] = rng
out["n_org"] = d.groupby(["tool", "type", "variant"]).org.nunique().groupby(["tool", "type"]).max()
out.round(4).reset_index().to_csv(D + "mixed_medians.tsv", sep="\t", index=False)
print(out.round(3).to_string())
