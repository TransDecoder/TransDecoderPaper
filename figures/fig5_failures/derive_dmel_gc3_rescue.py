#!/usr/bin/env python3
"""D. melanogaster coding-test failures (ORFs >= 300 nt) by GC3 quintile, with the number rescued
by the GC3-matched Markov model, for Fig. 5c.

Reads analyses/long_orf_failures/results/per_cds/dmel.per_cds.tsv.gz (not in git; made by
long_orf_failures.py) and applies the definitions of summarize.py (fail_score: category in
fail_coding_*, fail_reason S1<=0 or shifted_frame>=S1; quintiles over matched ORFs >= 300 nt).
Check: n and fail_score must equal analyses/long_orf_failures/results/fail_rate_by_feature.tsv
(dmel, gc3).

usage: derive_dmel_gc3_rescue.py > data/dmel_gc3_rescue.tsv
"""
import sys
from pathlib import Path

import pandas as pd

RES = Path(__file__).resolve().parents[2] / "analyses/long_orf_failures/results"
d = pd.read_csv(RES / "per_cds/dmel.per_cds.tsv.gz", sep="\t", low_memory=False,
                keep_default_na=False, na_values=["NA", ""])
m = d[(d.category != "no_orf") & (d.orf_len >= 300)].copy()
m["fail_score"] = (m.category.isin(["fail_coding_long_rescued", "fail_coding_lost"])
                   & m.fail_reason.isin(["S1<=0", "shifted_frame>=S1"]))
m["q"] = pd.qcut(m.gc3, 5, labels=[1, 2, 3, 4, 5])
rows = []
for q, g in m.groupby("q", observed=True):
    rows.append(dict(group=f"gc3_q{q}", lo=round(g.gc3.min(), 4), hi=round(g.gc3.max(), 4), n=len(g),
                     fail_score=int(g.fail_score.sum()),
                     rescued_gc3_model=int((g.fail_score & (g.pass_gc3match500 == 1)).sum())))
for name, sel in [("chr4", m.seqid == "4"), ("chrY", m.seqid == "Y"),
                  ("other_chr", ~m.seqid.isin(["4", "Y"]))]:
    g = m[sel]
    rows.append(dict(group=name, lo="NA", hi="NA", n=len(g), fail_score=int(g.fail_score.sum()),
                     rescued_gc3_model=int((g.fail_score & (g.pass_gc3match500 == 1)).sum())))
pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)
