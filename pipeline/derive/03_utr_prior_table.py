#!/usr/bin/env python3
"""Organism-independent UTR-length prior for start refinement, from the annotations only (E29).

usage: python 03_utr_prior_table.py <annot_dir> <out_dir>
  annot_dir: 03_utr_annotation.py outputs (<org>.tx.tsv.gz, <org>.trunc.tsv)

For d = the 5' UTR implied by trimming an open ORF to its first in-frame ATG:
  P(d | UTR)   lengths of annotated 5' UTRs (gene representatives) that are open (no in-frame
               stop) and whose annotated start is the first in-frame ATG -- the only UTRs for which
               the choice arises on a complete transcript
  P(d | trunc) distance from a uniformly placed cut inside an annotated CDS to the first in-frame
               ATG after the fragment's first codon
Each is binned, normalised per organism and averaged over the 8 organisms (equal weight). The prior
is log2 P(d|UTR)/P(d|trunc) in bits, given at the log-midpoint of each bin and clamped to [-6, 3].
Writes:
  annotation_utr5_summary.tsv  per organism: transcripts, fraction with a 5' UTR, quantiles,
                               fractions < 10/20/50 nt, for all UTRs and for open UTRs; for all
                               eligible transcripts ('all') and gene representatives ('rep')
  trunc_d_summary.tsv          per organism: fraction of cuts with d < 10/20/50 and no ATG
                               (uniform and N(0.5,0.1) cuts)
  utr_prior_lr.tsv             bin, lo, hi, mid, per-organism and pooled log2 LR (the prior)
"""
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

A, O = Path(sys.argv[1]), Path(sys.argv[2])
O.mkdir(parents=True, exist_ok=True)
x = pd.concat(pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(A / "*.tx.tsv.gz"))))
t = pd.concat(pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(A / "*.trunc.tsv"))))
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii"]

rows = []
for nm, sub in [("all", x), ("rep", x[x.truth_gene_rep == 1])]:
    for o in ORGS:
        g = sub[sub.org == o]
        for kind, u in [("any", g[g.utr5 > 0]), ("open", g[(g.utr5 > 0) & (g.open5 == 1)]),
                        ("open_first_atg", g[(g.utr5 > 0) & (g.open5 == 1) & (g.first_atg == 1)])]:
            rows.append(dict(set=nm, org=o, utr_kind=kind, n_tx=len(g), n=len(u),
                             frac_of_tx=round(len(u) / len(g), 4),
                             q10=u.utr5.quantile(.1), q25=u.utr5.quantile(.25), median=u.utr5.median(),
                             q75=u.utr5.quantile(.75), q90=u.utr5.quantile(.9),
                             lt10=round((u.utr5 < 10).mean(), 4), lt20=round((u.utr5 < 20).mean(), 4),
                             lt50=round((u.utr5 < 50).mean(), 4)))
pd.DataFrame(rows).to_csv(O / "annotation_utr5_summary.tsv", sep="\t", index=False)

rows = []
for o in ORGS:
    g = t[t.org == o]
    for w in ["uniform", "normal"]:
        tot = g[w].sum()
        ok = g[g.d >= 0]
        rows.append(dict(org=o, cut=w, no_atg=round(g[g.d < 0][w].sum() / tot, 4),
                         lt10=round(ok[ok.d < 10][w].sum() / tot, 4), lt20=round(ok[ok.d < 20][w].sum() / tot, 4),
                         lt50=round(ok[ok.d < 50][w].sum() / tot, 4),
                         median_d=float(np.interp(0.5, ok.sort_values("d")[w].cumsum() / ok[w].sum(),
                                                  ok.sort_values("d").d))))
pd.DataFrame(rows).to_csv(O / "trunc_d_summary.tsv", sep="\t", index=False)

B = [3, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 3000]
rep = x[x.truth_gene_rep == 1]
pu_all, pt_all, lr = [], [], {}
for o in ORGS:
    u = rep[(rep.org == o) & (rep.utr5 >= 3) & (rep.open5 == 1) & (rep.first_atg == 1)].utr5
    pu = np.histogram(u.clip(upper=2999), B)[0].astype(float)
    pu /= pu.sum()
    tt = t[(t.org == o) & (t.d >= 3)]
    pt = np.histogram(tt.d.clip(upper=2999), B, weights=tt.uniform.astype(float))[0].astype(float)
    pt /= pt.sum()
    pu_all.append(pu)
    pt_all.append(pt)
    with np.errstate(divide="ignore"):
        lr[o] = np.log2(pu / pt)
pu, pt = np.mean(pu_all, 0), np.mean(pt_all, 0)
pooled = np.clip(np.log2(np.maximum(pu, 1e-9) / pt), -6, 3)
tab = pd.DataFrame(dict(lo=B[:-1], hi=B[1:], mid=np.round(np.sqrt(np.array(B[:-1]) * np.array(B[1:])), 1),
                        p_utr=np.round(pu, 5), p_trunc=np.round(pt, 5), log2lr=np.round(pooled, 3)))
for o in ORGS:
    tab["lr_" + o] = np.round(np.clip(lr[o], -9, 9), 3)
tab.to_csv(O / "utr_prior_lr.tsv", sep="\t", index=False)
print(tab[["lo", "hi", "mid", "p_utr", "p_trunc", "log2lr"]].to_string(index=False))
