#!/usr/bin/env python3
"""Stop-codon sensitivity by true CDS length, and negatives with a prediction by predicted ORF length.

Shows where along the length axis a short-ORF rule gains or loses (e.g. the step at 100 codons of the
'cliff' rule). Matching as score_negatives.py (stop within TOL nt, same strand).
  coding rows: true CDS length bin (codons, within the input sequence), n_truth, tp_stop
  negative rows: predicted ORF length bin, number of negatives whose prediction falls in the bin

usage: sens_by_length.py <inputs_dir> <variants> label=<path template with {org},{variant}> ... > out.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
BINS = [50, 75, 100, 125, 150, 200, 268, 400, 10 ** 9]


def bin_of(nt):
    c = nt // 3
    for lo, hi in zip([0] + BINS[:-1], BINS):
        if c < hi:
            return f"{lo}-{hi - 1}" if hi < 10 ** 9 else f">={lo}"


def main():
    inputs, variants = sys.argv[1], sys.argv[2].split(",")
    labels = [a.split("=", 1) for a in sys.argv[3:]]
    print("org\tvariant\tlabel\tset\tlength_bin\tn\thits")
    for v in variants:
        for org in ORGS:
            tf = os.path.join(inputs, org, v + ".truth.tsv")
            if not os.path.exists(tf):
                continue
            truth = read_truth(tf)
            for label, tmpl in labels:
                pf = tmpl.format(org=org, variant=v)
                if not os.path.exists(pf) or os.path.getsize(pf) == 0:
                    continue
                preds = read_preds(pf)
                n, hit, neg = {}, {}, {}
                for tid, t in truth.items():
                    b = bin_of(t[1] - t[0] + 1)
                    n[b] = n.get(b, 0) + 1
                    t5, t3 = ends(*t[:3])
                    if any(p[2] == t[2] and abs(ends(*p)[1] - t3) <= TOL for p in preds.get(tid, [])):
                        hit[b] = hit.get(b, 0) + 1
                for tid, plist in preds.items():
                    if tid not in truth:
                        b = bin_of(max(p[1] - p[0] + 1 for p in plist))
                        neg[b] = neg.get(b, 0) + 1
                for b in sorted(n, key=lambda x: int(x.lstrip(">=").split("-")[0])):
                    print(f"{org}\t{v}\t{label}\tcoding\t{b}\t{n[b]}\t{hit.get(b, 0)}")
                for b in sorted(neg, key=lambda x: int(x.lstrip(">=").split("-")[0])):
                    print(f"{org}\t{v}\t{label}\tnegative\t{b}\t\t{neg[b]}")


if __name__ == "__main__":
    main()
