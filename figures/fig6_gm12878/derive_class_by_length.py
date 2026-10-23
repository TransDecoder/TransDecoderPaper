#!/usr/bin/env python3
"""Class composition of distinct predicted CDSs by predicted protein length (GM12878, final defaults).

Reads the per-prediction status tables written by eval_genome_cds.py (not in git):
  $P/work/gm12878_td/eval/status/<asm>_final_abinitio.pred_status.tsv
  columns: model_id, orf_type, aa, gencode_exact, gencode_compatible, no_gencode_cds_overlap, mane_exact
Classes as in gm12878_riboseq_support.tsv: exact; compatible (not exact); other_overlap (overlaps a
GENCODE CDS, not compatible); novel (no same-strand GENCODE CDS overlap).
Length bins: <100 aa (includes the few <50 aa left after start refinement), 100-149, 150-199,
200-299, >=300.

usage: derive_class_by_length.py [status_dir] > data/pred_class_by_length.tsv
"""
import csv
import sys
from collections import Counter

STATUS = sys.argv[1] if len(sys.argv) > 1 else \
    "/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/gm12878_td/eval/status"
BINS = [("lt100", 0, 100), ("100-149", 100, 150), ("150-199", 150, 200), ("200-299", 200, 300),
        ("ge300", 300, 10 ** 9)]
CLASSES = ["exact", "compatible", "other_overlap", "novel"]


def cls(r):
    if r["gencode_exact"] == "1":
        return "exact"
    if r["gencode_compatible"] == "1":
        return "compatible"
    if r["no_gencode_cds_overlap"] == "1":
        return "novel"
    return "other_overlap"


print("assembly\tlength_bin\tn\t" + "\t".join(CLASSES))
for asm in ("lraa", "stringtie", "trinity"):
    counts = {b[0]: Counter() for b in BINS}
    with open(f"{STATUS}/{asm}_final_abinitio.pred_status.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            aa = int(r["aa"])
            b = next(name for name, lo, hi in BINS if lo <= aa < hi)
            counts[b][cls(r)] += 1
    for name, _, _ in BINS:
        c = counts[name]
        print(f"{asm}\t{name}\t{sum(c.values())}\t" + "\t".join(str(c[k]) for k in CLASSES))
