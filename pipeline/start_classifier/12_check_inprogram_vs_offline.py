#!/usr/bin/env python3
"""E33 adoption: in-program homology start decisions (PyTD start_refinement.homology.tsv) vs the offline
08_homology_rules.py policy H, on the ORFs both evaluate.

usage: python 12_check_inprogram_vs_offline.py <features_dir> <offline_decisions_dir> <out.tsv> <dataset>=<run_dir> ...

Offline: <offline_decisions_dir>/<dataset>.decisions.tsv.gz column H (1 = keep open), on the step-1
candidates (5'-open ORFs with an ATG selected by the ab initio pytd_v7 run), with the rule inputs
and ab initio decision from <features_dir>/<dataset>.tsv.gz (trim_run).
In-program: the homology run's table (5'-open ORFs with an ATG selected by that run; homology changes
the selection, so only ORF ids present in both are compared). Also compares the ab initio decision.
"""
import sys

import pandas as pd

fdir, odir, out = sys.argv[1:4]
rows = []
for arg in sys.argv[4:]:
    ds, run = arg.split("=", 1)
    off = pd.read_csv(f"{odir}/{ds}.decisions.tsv.gz", sep="\t", usecols=["orf_id", "H"])
    feat = pd.read_csv(f"{fdir}/{ds}.tsv.gz", sep="\t", usecols=["orf_id", "trim_run"])
    off = off.merge(feat, on="orf_id")
    inp = pd.read_csv(f"{run}/input.fa.transdecoder_dir/start_refinement.homology.tsv", sep="\t")
    m = off.merge(inp, on="orf_id")
    keep_in = (m.decision == "keep").astype(int)
    ab_in = (m.ab_initio_decision == "trim").astype(int)
    mis = m[keep_in != m.H.astype(int)]
    rows.append(dict(dataset=ds, offline_candidates=len(off), inprogram_candidates=len(inp), compared=len(m),
                     decision_mismatches=len(mis), ab_initio_mismatches=int((ab_in != m.trim_run).sum()),
                     homology_decided=int((m.source == "homology").sum()),
                     conflicts=int((m.source == "ab_initio(conflict)").sum())))
    if len(mis):
        mis.to_csv(f"{out}.{ds}.mismatches.tsv", sep="\t", index=False)
r = pd.DataFrame(rows)
r.to_csv(out, sep="\t", index=False)
print(r.to_string(index=False))
print("TOTAL compared", r.compared.sum(), "decision mismatches", r.decision_mismatches.sum(),
      "ab initio mismatches", r.ab_initio_mismatches.sum())
