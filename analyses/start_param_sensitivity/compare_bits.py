#!/usr/bin/env python3
"""Kozak score in bits (weight 17 per bit) vs the previous nats rule (weight 25 per nat).

The previous runs are the E13 grid cell a0.4_w25 (nats code, PyTD 48c870c). The bits code is run with
weight 17 (the new default) and 17.32868 (= 25 x ln 2, the exact equivalent; differences then come only
from the Kozak model file being written to 3 decimals). Reports, per input, the predictions whose
coordinates differ and the exact F1 of each.

usage: compare_bits.py [grid_dir_nats] [grid_dir_bits] > results/bits_vs_nats.tsv
"""
import csv
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/start_param_grid"
NATS = sys.argv[1] if len(sys.argv) > 1 else f"{P}/runs"
BITS = sys.argv[2] if len(sys.argv) > 2 else f"{P}/runs_bits"
ORGS = "athal creinhardtii dmel human mouse pfalciparum spombe tthermophila".split()


def preds(path):
    with open(path) as fh:
        return {tuple(r[k] for k in ("transcript_id", "start", "end", "strand"))
                for r in csv.DictReader(fh, delimiter="\t")}


def f1(path):
    with open(path) as fh:
        return float(next(csv.DictReader(fh, delimiter="\t"))["f1_exact"])


print("org\tinput\tn_pred_nats25\tdiffer_bits17\tdiffer_bits17.33\tf1_exact_nats25\tf1_exact_bits17\tf1_exact_bits17.33")
for org in ORGS:
    for inp in ("complete", "fragments"):
        old = f"{NATS}/{org}/{inp}/a0.4_w25"
        new17 = f"{BITS}/{org}/{inp}/a0.4_w17"
        newx = f"{BITS}/{org}/{inp}/a0.4_w17.3287"
        a, b, c = preds(f"{old}/preds.tsv"), preds(f"{new17}/preds.tsv"), preds(f"{newx}/preds.tsv")
        print(f"{org}\t{inp}\t{len(a)}\t{len(a - b)}\t{len(a - c)}\t"
              f"{f1(f'{old}/score.tsv'):.4f}\t{f1(f'{new17}/score.tsv'):.4f}\t{f1(f'{newx}/score.tsv'):.4f}")
