#!/usr/bin/env python3
"""E32: Arabidopsis negatives split by source class (TAIR12 ncRNA_class lncRNA vs antisense_RNA).

usage: 14_negatives_by_class.py <outfile> <level> [<level> ...]
For each level (d1 = full depth), the 'noncoding' contigs of label_assembly.py are assigned to a class from
their compatible source transcripts (lncRNA, antisense_RNA, or 'both' when sources of both classes are
compatible), and for each tool the number of those contigs with >= 1 prediction is counted (preds.tsv of the
run, IDs cut at the first whitespace as in score_simtrinity.py; the longest-ORF baseline from pytd_v7).
Writes: level, class, n_negatives, tool, with_pred, frac_with_pred.
"""
import collections
import gzip
import os
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
TOOLS = ["pytd_v7", "perl_td_v600", "prodigal", "gmst", "cpat", "cpc2", "rnasamba", "borf", "longest_orf"]


def main():
    outf, levels = sys.argv[1], sys.argv[2:]
    ttype = {}
    with gzip.open(f"{P}/data/simtrinity_athal/ref/tair12.tx.tsv.gz", "rt") as fh:
        fh.readline()
        for line in fh:
            f = line.split("\t")
            ttype[f[0]] = f[4]
    rows = []
    for lev in levels:
        lab = f"{P}/data/simtrinity_athal/depth/{lev}/labels/sim_trinity.labels.tsv"
        runs = f"{P}/work/simtrinity_athal/depth/{lev}/runs"
        cls = {}
        with open(lab) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if r["label"] != "noncoding":
                    continue
                c = {ttype[s] for s in r["source_transcripts"].split(",")}
                cls[r["transcript_id"]] = c.pop() if len(c) == 1 else "both"
        n = collections.Counter(cls.values())
        for tool in TOOLS:
            pf = f"{runs}/pytd_v7/longest_orf.preds.tsv" if tool == "longest_orf" else f"{runs}/{tool}/preds.tsv"
            hit = set()
            with open(pf) as fh:
                for line in fh:
                    if line.strip() and not line.startswith(("#", "transcript_id\t")):
                        hit.add(line.split("\t", 1)[0].split()[0])
            h = collections.Counter(c for t, c in cls.items() if t in hit)
            for c in sorted(n):
                rows.append([lev, c, n[c], tool, h[c], f"{h[c] / n[c]:.4f}"])
    with open(outf, "w") as fh:
        fh.write("level\tclass\tn_negatives\ttool\twith_pred\tfrac_with_pred\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    sys.exit(main())
