#!/usr/bin/env python3
"""Homology-assisted vs ab initio PyTD at the final defaults, per benchmark input.

Compares a homology run with its ab initio run (default pytd_dev_v2_homology vs pytd_dev_v2; the
release runs are pytd_v7_homology vs pytd_v7) in a summary.tsv from score_all.py. Prints the
number of predictions of each, the change (homology minus ab initio) in sensitivity, precision and
F1 at the stop and exact levels, and the predictions on negatives: rand_* transcripts for .wRand
inputs (n_pred_random) and, for .wLnc inputs, transcripts absent from the truth table (lncRNAs),
counted from preds.tsv.

usage: compare_homology.py <summary.tsv> <inputs_dir> <runs_dir> [base_tool homology_tool] > homology_vs_abinitio.tsv
"""
import csv
import os
import sys

BASE, HOM = "pytd_dev_v2", "pytd_dev_v2_homology"
METRICS = ["sensitivity_stop", "precision_stop", "f1_stop",
           "sensitivity_exact", "precision_exact", "f1_exact"]


def lnc_with_pred(inputs, runs, org, variant, tool):
    with open(os.path.join(inputs, org, variant + ".truth.tsv")) as fh:
        truth = {r["transcript_id"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(os.path.join(runs, org, variant, tool, "preds.tsv")) as fh:
        return len({r["transcript_id"] for r in csv.DictReader(fh, delimiter="\t")} - truth)


def main():
    global BASE, HOM
    summary, inputs, runs = sys.argv[1:4]
    if len(sys.argv) > 4:
        BASE, HOM = sys.argv[4:6]
    with open(summary) as fh:
        rows = {(r["org"], r["variant"], r["tool"]): r for r in csv.DictReader(fh, delimiter="\t")}
    keys = sorted((k for k in rows if k[2] == HOM), key=lambda k: (k[1], k[0]))
    print("org\tvariant\tn_pred_abinitio\tn_pred_homology\t" + "\t".join("d_" + m for m in METRICS)
          + "\tf1_exact_abinitio\tf1_exact_homology\tneg_pred_abinitio\tneg_pred_homology")
    for org, variant, _ in keys:
        b, h = rows[(org, variant, BASE)], rows[(org, variant, HOM)]
        if variant.endswith(".wLnc"):
            neg = [lnc_with_pred(inputs, runs, org, variant, t) for t in (BASE, HOM)]
        elif ".wRand" in variant:
            neg = [b["n_pred_random"], h["n_pred_random"]]
        else:
            neg = ["NA", "NA"]
        print(f"{org}\t{variant}\t{b['n_pred']}\t{h['n_pred']}\t"
              + "\t".join(f"{float(h[m]) - float(b[m]):+.4f}" for m in METRICS)
              + f"\t{b['f1_exact']}\t{h['f1_exact']}\t{neg[0]}\t{neg[1]}")


if __name__ == "__main__":
    sys.exit(main())
