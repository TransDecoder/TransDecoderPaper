#!/usr/bin/env python3
"""E35: what the changed decisions look like (Kozak score, |u|, coding per base), by arm, group, variant and class.

usage: python 05_changed_profile_v1.py <eval_dir>/changed_candidates.tsv.gz <out.tsv>

Rows: arm x group x variant x class (w2r, r2w, w2w_later, other); n, median K (bits), median |u| (nt),
median coding(u) per base of u, share with |u| < 60 nt, and veto_precision = r2w / (w2r + r2w + w2w_later):
the share of changed decisions with a stop-matched truth CDS where the release keep-open was right.
"""
import sys

import numpy as np
import pandas as pd

C = pd.read_csv(sys.argv[1], sep="\t")
C["arm"] = np.where(C.dataset.str.endswith(".H"), "homology", "abinitio")
grp = C.dataset.str.replace(r"\.H$", "", regex=True)
C["group"] = np.select([grp.str.contains(r"\.complete_fullutr$"), grp.str.contains(r"\.complete$"),
                        grp.str.contains(r"\.fragments_mixed"), grp.str.contains(r"\.fragments$"),
                        grp.str.startswith("sim.human"), grp.str.startswith("sim.athal")],
                       ["bench_complete_fullutr", "bench_complete", "bench_fragments_mixed", "bench_fragments",
                        "sim_human", "sim_athal"], "other")
C["c"] = np.where(C.cls.isin(["w2r", "r2w", "w2w_later"]), C.cls, "other")
C["coding_pb"] = C.coding / C.seg
rows = []
for (arm, g, v), X in C.groupby(["arm", "group", "variant"], sort=False):
    n = X.c.value_counts()
    lab = n.get("w2r", 0) + n.get("r2w", 0) + n.get("w2w_later", 0)
    for c, Y in X.groupby("c"):
        rows.append(dict(arm=arm, group=g, variant=v, cls=c, n=len(Y), median_K=Y.kozak.median(),
                         median_seg=Y.seg.median(), median_coding_pb=Y.coding_pb.median(),
                         share_seg_lt60=(Y.seg < 60).mean(),
                         veto_precision=n.get("r2w", 0) / lab if lab else np.nan))
pd.DataFrame(rows).to_csv(sys.argv[2], sep="\t", index=False, float_format="%.3f")
print("done", len(rows))
