#!/usr/bin/env python3
"""Sensitivity on coding transcripts and calls on negatives, per input and tool/setting.

For inputs that mix coding transcripts with negatives (random 'rand_' sequences or lncRNAs, i.e.
sequences without a truth row), reports per organism, variant and label:
  n_truth, tp_stop, tp_exact, n_pred, n_neg, neg_with_pred
  sens_stop, sens_exact, precision_stop (all predictions, as score_predictions.py), neg_rate
Matching and counts follow score_predictions.py (one truth credited per transcript).

usage: score_negatives.py <inputs_dir> <variants> label=<path template> ... > out.tsv
  the template is formatted with {org} and {variant}, e.g.
  s100_c200=work/bench_v3/stringency/{org}/{variant}/s100_c200/preds.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]


def seq_ids(fasta):
    with open(fasta) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def main():
    inputs, variants = sys.argv[1], sys.argv[2].split(",")
    labels = [a.split("=", 1) for a in sys.argv[3:]]
    print("org\tvariant\tlabel\tn_truth\ttp_stop\ttp_exact\tn_pred\tn_neg\tneg_with_pred\t"
          "sens_stop\tsens_exact\tprecision_stop\tneg_rate")
    for v in variants:
        for org in ORGS:
            fa = os.path.join(inputs, org, v + ".fa")
            if not os.path.exists(fa):
                continue
            truth = read_truth(os.path.join(inputs, org, v + ".truth.tsv"))
            n_neg = sum(1 for k in seq_ids(fa) if k not in truth)
            for label, tmpl in labels:
                pf = tmpl.format(org=org, variant=v)
                if not os.path.exists(pf) or os.path.getsize(pf) == 0:
                    continue
                preds = read_preds(pf)
                n_pred = tps = tpe = negp = 0
                for tid, plist in preds.items():
                    n_pred += len(plist)
                    t = truth.get(tid)
                    if not t:
                        negp += 1
                        continue
                    best = None
                    for p in plist:
                        if p[2] != t[2]:
                            continue
                        p5, p3 = ends(*p)
                        t5, t3 = ends(*t[:3])
                        if abs(p3 - t3) <= TOL:
                            best = "exact" if abs(p5 - t5) <= TOL or best == "exact" else "stop"
                    if best:
                        tps += 1
                        tpe += best == "exact"
                n_truth = len(truth)
                print(f"{org}\t{v}\t{label}\t{n_truth}\t{tps}\t{tpe}\t{n_pred}\t{n_neg}\t{negp}\t"
                      f"{tps / n_truth:.4f}\t{tpe / n_truth:.4f}\t{tps / n_pred if n_pred else 0:.4f}\t"
                      f"{negp / n_neg if n_neg else 0:.4f}")


if __name__ == "__main__":
    main()
