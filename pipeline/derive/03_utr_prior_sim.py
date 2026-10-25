#!/usr/bin/env python3
"""Offline effect of a UTR-length prior on the release start-refinement decisions (evidence E29).

usage: python 03_utr_prior_sim.py <features_dir> <prior_lr.tsv> <out.tsv> [tag=release] [inputs]

The release trims an open ORF to its first in-frame ATG when alpha*seg + w*K > coding. The prior adds
lambda * w * L(d) to the UTR side, L(d) = log2 P(d|UTR)/P(d|trunc) from utr_prior_lr.tsv (pooled,
annotation-derived; interpolated linearly in log d between bin midpoints, constant beyond them).
Because stops do not change, the effect on exact-CDS matches is computed directly from the decisions
on ORFs whose stop matches the reference: 5'-partial reference -> right iff kept open; start present
and the first ATG is the true start -> right iff trimmed; otherwise wrong either way.
Exact F1 = 2 TP / (predictions + references). Rows: org, input, ref_type, variant, tp_exact, delta,
gained, lost, exact F1 (input level, all reference types) for each lambda and variant:
  full      lambda * w * L(d)
  short     only for d < 20 (the 'very short UTR' hypothesis), L clamped at <= 0
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FD, LRF, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
TAG = sys.argv[4] if len(sys.argv) > 4 else "release"
INPUTS = (sys.argv[5] if len(sys.argv) > 5 else "complete,complete_fullutr,fragments").split(",")
W = 17.0
TRUNC = {"5prime_partial", "internal"}  # references missing the start codon
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii", "gm12878"]
PP = Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER")


def truth_path(org, inp):
    """benchmark truth, or for org 'gm12878' the assembly truth of 03_utr_gm12878_truth.py"""
    if org == "gm12878":
        return PP / f"work/utr_prior/gm12878/{inp}.truth.tsv"
    return PP / f"data/bench_v3/inputs/{org}/{inp}.truth.tsv"
P = Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER")
lr = pd.read_csv(LRF, sep="\t")
xm, ym = np.log(lr.mid.values), lr.log2lr.values


def L(d):
    return np.interp(np.log(np.maximum(d, 1)), xm, ym)


VARIANTS = [("release", 0, "full")] + [(f"full_l{l}", l, "full") for l in (0.5, 1, 2)] + [("short_l1", 1, "short")]
rows = []
for org in ORGS:
    for inp in INPUTS:
        f = FD / f"{TAG}.{org}.{inp}.tsv.gz"
        if not f.exists():
            continue
        x = pd.read_csv(f, sep="\t")
        truth = pd.read_csv(truth_path(org, inp), sep="\t")
        pf = P / f"work/bench_v3/runs/{org}/{inp}/pytd_v7/preds.tsv"
        if not pf.exists():  # mixed inputs: the lambda = 0 run of 03_utr_runs.sh (= release)
            pf = P / f"work/utr_prior/runs/{org}/{inp}/utr_l0/preds.tsv"
        npred = sum(1 for l in open(pf) if not l.startswith("transcript_id")) if pf.exists() else 0
        x = x[x.truth_type.notna()]
        s = x[x.stop_match == 1].drop_duplicates("transcript_id").copy()
        s["trim"] = s.trim.fillna(0).astype(int)
        pred5 = np.where(s.trim == 1, s.d + 1, s.orf5)
        s["exact0"] = (np.abs(pred5 - s.truth5) <= 2)
        dec = s.open5.eq(1) & s.d.notna()
        tp0 = int(s.exact0.sum())
        o = s[dec].copy()
        right_if_trim = (~o.truth_type.isin(TRUNC)) & ((o.d + 1 - o.truth5).abs() <= 2)
        right_if_keep = o.truth_type.isin(TRUNC)
        for name, lam, kind in VARIANTS:
            Ld = L(o.d.values)
            if kind == "short":
                Ld = np.where(o.d.values < 20, np.minimum(Ld, 0), 0)
            u = o.utr_score.values + lam * W * Ld
            trim = u > o.coding.values
            right = np.where(trim, right_if_trim, right_if_keep)
            right0 = np.where(o.trim.values == 1, right_if_trim, right_if_keep)
            for tt in ["all"] + sorted(s.truth_type.unique()):
                m = np.ones(len(o), bool) if tt == "all" else (o.truth_type == tt).values
                ns = len(s) if tt == "all" else int((s.truth_type == tt).sum())
                base = tp0 if tt == "all" else int(s[s.truth_type == tt].exact0.sum())
                g = int((right & ~right0 & m).sum())
                lo = int((~right & right0 & m).sum())
                tp = base + g - lo
                nref = len(truth) if tt == "all" else int((truth.type == tt).sum())
                row = dict(org=org, input=inp, ref_type=tt, variant=name, n_ref=nref, n_stop=ns,
                           tp_exact=tp, gained=g, lost=lo, delta=g - lo,
                           sens_exact=round(tp / nref, 4))
                if tt == "all":
                    row["f1_exact"] = round(2 * tp / (npred + nref), 4) if npred else float("nan")
                rows.append(row)
pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
