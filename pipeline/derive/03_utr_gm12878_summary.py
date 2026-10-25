#!/usr/bin/env python3
"""5' truncation and 5' UTRs in GM12878 assemblies vs the benchmark fragments (evidence E29).

usage: python 03_utr_gm12878_summary.py <gm_dir> <out_dir>
  gm_dir: 03_utr_gm12878_truth.py outputs (<asm>.pairs.tsv.gz for lraa, stringtie, trinity)
Writes
  gm12878_five_prime_truncation.tsv  per assembly: transcripts that carry a MANE Select stop with the
      reference structure from the stop to their 5' end; how many begin inside the CDS (5' partial)
      and how much CDS they miss; of those that contain the start, how many have an open 5' UTR (no
      in-frame stop) and its length; pi = 5' partials / (5' partials + open UTRs of >= 3 nt), the share
      of truly truncated ORFs among open 5' ends. The benchmark rows give the same for the stranded
      fragments, complete and mixed inputs of the 8 organisms (from the truth tables).
  missing_frac_hist.tsv  distribution (10 bins) of the fraction of the CDS missing at the 5' end, for
      the GM12878 5' partials and the benchmark 5' partials (fragments, N(0.5,0.1) cut; mixed inputs)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

G, O = Path(sys.argv[1]), Path(sys.argv[2])
P = Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER")
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii"]
BINS = np.linspace(0, 1, 11)
rows, hist = [], []


def summarise(label, c, p, op):
    """c: complete (utr5), p: 5' partial (missing_frac), op: open-UTR subset of c"""
    n = len(c) + len(p)
    rows.append(dict(source=label, n=n, n_5p=len(p), frac_5p=round(len(p) / n, 4),
                     missing_frac_q25=round(p.missing_frac.quantile(.25), 3) if len(p) else "",
                     missing_frac_median=round(p.missing_frac.median(), 3) if len(p) else "",
                     missing_frac_q75=round(p.missing_frac.quantile(.75), 3) if len(p) else "",
                     frac_5p_missing_lt10pct=round((p.missing_frac < .1).mean(), 4) if len(p) else "",
                     n_complete=len(c), n_open_utr=len(op), frac_open=round(len(op) / len(c), 4) if len(c) else "",
                     open_utr_median=op.utr5.median() if len(op) else "",
                     open_utr_lt10=round((op.utr5 < 10).mean(), 4) if len(op) else "",
                     open_utr_lt20=round((op.utr5 < 20).mean(), 4) if len(op) else "",
                     open_utr_lt50=round((op.utr5 < 50).mean(), 4) if len(op) else "",
                     pi_trunc_among_open=round(len(p) / (len(p) + (op.utr5 >= 3).sum()), 4)
                     if len(p) + len(op) else ""))
    if len(p):
        h = np.histogram(p.missing_frac, BINS)[0]
        for lo, k in zip(BINS[:-1], h):
            hist.append(dict(source=label, bin_lo=round(lo, 1), n=int(k), frac=round(k / len(p), 4)))


for asm in ["lraa", "stringtie", "trinity"]:
    x = pd.read_csv(G / f"{asm}.pairs.tsv.gz", sep="\t")
    c = x[x.status == "complete"]
    summarise(f"gm12878_{asm}", c, x[x.status == "5prime_partial"], c[c.open5 == 1])

STOPS = {1: {"TAA", "TAG", "TGA"}, 6: {"TGA"}}


def read_fa(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


for v in ["complete", "fragments", "fragments_mixed_t1"]:
    cs, ps, ops = [], [], []
    for org in ORGS:
        full = pd.read_csv(P / f"data/bench_v3/inputs/{org}/complete.truth.tsv", sep="\t").set_index("transcript_id")
        t = pd.read_csv(P / f"data/bench_v3/inputs/{org}/{v}.truth.tsv", sep="\t")
        seqs = read_fa(P / f"data/bench_v3/inputs/{org}/{v}.fa")
        stops = STOPS[int(t.transl_table.iloc[0])]
        flen = (full.cds_end - full.cds_start + 1).reindex(t.transcript_id).values
        t["missing_frac"] = 1 - (t.cds_end - t.cds_start + 1) / flen
        t["utr5"] = t.cds_start - 1
        p = t[t.type.isin(["5prime_partial", "internal"])]
        c = t[t.type.isin(["complete", "3prime_partial"])].copy()
        c["open5"] = [int(not any(seqs[r.transcript_id][i:i + 3] in stops for i in range(r.utr5 - 3, -1, -3)))
                      for r in c.itertuples()]
        cs.append(c)
        ps.append(p)
        ops.append(c[c.open5 == 1])
    summarise(f"benchmark_{v}", pd.concat(cs), pd.concat(ps), pd.concat(ops))

pd.DataFrame(rows).to_csv(O / "gm12878_five_prime_truncation.tsv", sep="\t", index=False)
pd.DataFrame(hist).to_csv(O / "missing_frac_hist.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string())
