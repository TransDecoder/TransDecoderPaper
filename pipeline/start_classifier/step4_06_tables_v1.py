#!/usr/bin/env python3
"""E33 step 4: compact tables for the write-up.

usage: python step4_06_tables_v1.py <step4_big_dir> <out_dir>

Reads <step4_big_dir>/anatomy/anatomy_per_cds.tsv.gz, evidence.tsv, variants/variants_f1.tsv.gz and writes to
<out_dir>: anatomy_summary.tsv (dataset x truth type x class, with reasons folded into the class name),
evidence_summary.tsv, family_by_dataset.tsv (each family member's change in exact F1 per dataset group:
tuning mean, d32, athal, bench fragments median, bench complete median / min / n down).
"""
import os
import sys

import numpy as np
import pandas as pd

big, out = sys.argv[1:3]
os.makedirs(out, exist_ok=True)

D = pd.read_csv(f"{big}/anatomy/anatomy_per_cds.tsv.gz", sep="\t", keep_default_na=False)
D["klass"] = np.where(D.cls.isin(["selection", "hit_no_rule"]), D.cls + ":" + D.reason,
                      np.where(D.cls == "rule_wrong", "rule_wrong:" + D.rules,
                               np.where(D.cls == "start_elsewhere", "start_elsewhere", D.cls)))
A = D.groupby(["dataset", "type", "klass"]).size().unstack("klass", fill_value=0)
A.insert(0, "total_not_exact", A.sum(axis=1))
A.to_csv(f"{out}/anatomy_summary.tsv", sep="\t")
# complete CDSs wrongly left open: right stop, 5'-open candidate, trim needed, kept open
lo = D[(D.type == "complete") & (D.needed == "trim") & (D.decision == "keep")]
lo.groupby(["dataset", "cls"]).size().unstack("cls", fill_value=0).to_csv(f"{out}/anatomy_complete_left_open.tsv", sep="\t")

E = pd.read_csv(f"{big}/evidence.tsv", sep="\t")
E.to_csv(f"{out}/evidence_summary.tsv", sep="\t", index=False)

V = pd.read_csv(f"{big}/variants/variants_f1.tsv.gz", sep="\t")
BASE = "k1|1e-05|best|top1|ga"
b = V[V.variant == BASE].set_index("dataset").f1_exact
V["d"] = 100 * (V.f1_exact - V.dataset.map(b))
V = V[~V.variant.isin(["release", "perfect"])]
p = V.variant.str.split("|", expand=True)
p.columns = ["search", "evalue", "hsp", "depth", "pfam"]
V = pd.concat([V.drop(columns=[c for c in p.columns if c in V.columns]), p.set_index(V.index)], axis=1)
fam = {
    "evalue_x_search": (V.hsp == "best") & (V.depth == "top1") & (V.pfam == "ga"),
    "depth_scheme_fast_1e-5": (V.search == "fast") & (V.evalue == "1e-05") & (V.hsp == "best") & (V.pfam == "ga"),
    "depth_scheme_msens_10_all": (V.search == "msens") & (V.evalue == "10") & (V.hsp == "all") & (V.pfam == "ga"),
    "hsp": (V.evalue == "1e-05") & (V.depth == "top1") & (V.pfam == "ga"),
    "pfam_k1": (V.search == "k1"),
    "pfam_best_blast": (V.search == "msens") & (V.evalue == "10") & (V.hsp == "all") & (V.depth == "top3_bits"),
}
rows = []
for f, msk in fam.items():
    for v, x in V[msk].groupby("variant"):
        x = x.set_index("dataset").d
        bc = x[[i for i in x.index if i.startswith("bench.") and i.endswith(".complete")]]
        bf = x[[i for i in x.index if i.startswith("bench.") and i.endswith(".fragments")]]
        rows.append(dict(family=f, variant=v, tune=x.reindex(["sim.human.d1", "sim.human.d8"]).mean(),
                         d1=x.get("sim.human.d1"), d8=x.get("sim.human.d8"), d32=x.get("sim.human.d32"),
                         athal_d1=x.get("sim.athal.d1"), bench_frag_median=bf.median(), bench_frag_n=len(bf),
                         bench_frag_n_down=int((bf < -1e-9).sum()), bench_complete_median=bc.median(),
                         bench_complete_min=bc.min(), bench_complete_n=len(bc), bench_complete_n_down=int((bc < -1e-9).sum())))
F = pd.DataFrame(rows).sort_values(["family", "tune"], ascending=[True, False])
F.to_csv(f"{out}/family_by_dataset.tsv", sep="\t", index=False, float_format="%.3f")
print(F.to_string(float_format=lambda z: f"{z:.3f}"))
