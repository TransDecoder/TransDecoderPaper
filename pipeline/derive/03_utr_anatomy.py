#!/usr/bin/env python3
"""5'-end error anatomy of start refinement on the stranded benchmark inputs (evidence E29).

usage: python 03_utr_anatomy.py <features_dir> <out_dir> [tag=release] [variants=complete,complete_fullutr,fragments]

Reads <features_dir>/<tag>.<org>.<variant>.tsv.gz (03_utr_features.py: the selected ORFs before
start refinement, with the release decision re-derived and checked against the run) and writes:

  <tag>.anatomy.tsv   per organism x input x reference type: references, references whose stop is
                      predicted (stop-level), and how the 5' end of that prediction relates to the
                      reference 5' end:
                        exact            within 2 nt (open at the fragment start for 5'-partials)
                        trim_trunc       5'-partial reference trimmed to an internal ATG
                        open_utr         start-present reference left open (extends through the
                                         true ATG to the sequence start)
                        trim_upstream    trimmed to an in-frame ATG upstream of the true start
                        atg_other        ORF begins (LongOrfs) at an ATG that is not the true start
                        other            anything else
                      and for trim_trunc the median d (nt between fragment start and chosen ATG) and
                      trimmed fraction of the reference CDS; for open_utr the median true UTR length.
  <tag>.decision_by_d.tsv  ORFs open at 5' with an in-frame ATG and a correct stop, by d (the 5'
                      UTR implied by trimming to that ATG): class 'truncated' (5'-partial
                      reference: keeping open is right), 'utr_first' (start-present reference whose
                      true start is that ATG: trimming is right), 'utr_later' (true start is a later
                      ATG: neither is right); n and number trimmed.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FD, OD = Path(sys.argv[1]), Path(sys.argv[2])
TAG = sys.argv[3] if len(sys.argv) > 3 else "release"
VARIANTS = (sys.argv[4] if len(sys.argv) > 4 else "complete,complete_fullutr,fragments").split(",")
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii", "gm12878"]
PP = Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER")


def truth_path(org, inp):
    """benchmark truth, or for org 'gm12878' the assembly truth of 03_utr_gm12878_truth.py"""
    if org == "gm12878":
        return PP / f"work/utr_prior/gm12878/{inp}.truth.tsv"
    return PP / f"data/bench_v3/inputs/{org}/{inp}.truth.tsv"
TRUNC = {"5prime_partial", "internal"}  # references missing the start codon
DBINS = [0, 10, 20, 50, 100, 200, 500, 10 ** 9]
DLAB = ["0-9", "10-19", "20-49", "50-99", "100-199", "200-499", ">=500"]
OD.mkdir(parents=True, exist_ok=True)


def classify(r):
    pred5 = r.d + 1 if r.trim == 1 else r.orf5
    if abs(pred5 - r.truth5) <= 2:
        return "exact"
    if r.truth_type in TRUNC:
        return "trim_trunc" if r.trim == 1 else "other"
    if r.open5 == 1 and r.trim != 1:
        return "open_utr"
    if r.trim == 1 and pred5 < r.truth5:
        return "trim_upstream"
    if r.open5 == 0:
        return "atg_other"
    return "other"


anat, dec = [], []
for org in ORGS:
    for v in VARIANTS:
        f = FD / f"{TAG}.{org}.{v}.tsv.gz"
        if not f.exists():
            continue
        x = pd.read_csv(f, sep="\t", dtype={"trim": "Int64", "d": "Int64"})
        x = x[x.truth_type.notna()]
        ref_n = x.drop_duplicates("transcript_id")  # transcripts with any prediction
        tfile = truth_path(org, v)
        truth = pd.read_csv(tfile, sep="\t")
        s = x[x.stop_match == 1].drop_duplicates("transcript_id").copy()
        s["trim"] = s.trim.fillna(0)
        s["cls"] = [classify(r) for r in s.itertuples()]
        s["true_utr"] = s.truth5 - 1
        s["reflen"] = s.truth3 - s.truth5 + 1
        for tt, g in s.groupby("truth_type"):
            row = dict(org=org, input=v, ref_type=tt, n_ref=int((truth.type == tt).sum()), n_stop=len(g))
            for c in ["exact", "trim_trunc", "open_utr", "trim_upstream", "atg_other", "other"]:
                row[c] = int((g.cls == c).sum())
            tr = g[g.cls == "trim_trunc"]
            row["trim_trunc_median_d"] = float(tr.d.median()) if len(tr) else np.nan
            row["trim_trunc_frac_d_lt20"] = round(float((tr.d < 20).mean()), 4) if len(tr) else np.nan
            row["trim_trunc_frac_d_lt50"] = round(float((tr.d < 50).mean()), 4) if len(tr) else np.nan
            row["trim_trunc_median_trimmed_frac"] = round(float(((tr.d + 1 - tr.truth5) / tr.reflen).median()), 4) if len(tr) else np.nan
            ou = g[g.cls == "open_utr"]
            row["open_utr_median_true_utr"] = float(ou.true_utr.median()) if len(ou) else np.nan
            row["open_utr_frac_true_utr_lt20"] = round(float((ou.true_utr < 20).mean()), 4) if len(ou) else np.nan
            anat.append(row)
        # decision set
        o = s[(s.open5 == 1) & s.d.notna()].copy()
        o["dclass"] = np.where(o.truth_type.isin(TRUNC), "truncated",
                               np.where((o.d + 1 - o.truth5).abs() <= 2, "utr_first", "utr_later"))
        o["dbin"] = pd.cut(o.d.astype(float), DBINS, right=False, labels=DLAB)
        for (dc, db), g in o.groupby(["dclass", "dbin"], observed=True):
            dec.append(dict(org=org, input=v, dclass=dc, dbin=db, n=len(g), n_trim=int(g.trim.sum())))

pd.DataFrame(anat).to_csv(OD / f"{TAG}.anatomy.tsv", sep="\t", index=False)
pd.DataFrame(dec).to_csv(OD / f"{TAG}.decision_by_d.tsv", sep="\t", index=False)
