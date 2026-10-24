#!/usr/bin/env python3
"""Oracle log prior odds that a candidate ORF of a given length is the annotated CDS (stringency study, E19).

From a LongOrfs run (longest_orfs.gff3) and the truth table: candidates are 'true' when on the annotated
strand with the annotated stop (within TOL nt), else 'false'. Per BIN-codon bin: log((true+0.5)/(false+0.5)),
made non-decreasing in length (running maximum) and written at bin midpoints, for Predict
--short-orf-rule curve --short-orf-curve. Uses the annotation, so it is an upper bound for a curve
estimated from the data alone.

usage: oracle_logodds_curve.py <longest_orfs.gff3> <truth.tsv> > curve.tsv   (columns: codons, logodds)
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, ends, read_truth  # noqa: E402

BIN, MAX = 10, 600


def main():
    truth = read_truth(sys.argv[2])
    nb = MAX // BIN
    t, f = [0] * (nb + 1), [0] * (nb + 1)
    with open(sys.argv[1]) as fh:
        for line in fh:
            x = line.split("\t")
            if len(x) < 9 or x[2] != "CDS":
                continue
            s, e, st = int(x[3]), int(x[4]), x[6]
            b = min((e - s + 1) // 3 // BIN, nb)
            tr = truth.get(x[0])
            if tr and st == tr[2] and abs(ends(s, e, st)[1] - ends(*tr[:3])[1]) <= TOL:
                t[b] += 1
            else:
                f[b] += 1
    print("codons\tlogodds")
    best = -math.inf
    for b in range(nb + 1):
        if t[b] + f[b] == 0:
            continue
        best = max(best, math.log((t[b] + 0.5) / (f[b] + 0.5)))
        print(f"{b * BIN + BIN / 2:g}\t{best:.4f}")


if __name__ == "__main__":
    main()
