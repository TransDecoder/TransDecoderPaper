#!/usr/bin/env python3
"""Predictions on lncRNA negatives that match translated ORFs of the GENCODE Ribo-seq ORF catalog.

Annotated lncRNAs are not all noncoding: some carry translated small ORFs. For the lncRNA negatives of
an input (sequences without a truth row), catalog ORFs (TransCODE/GENCODE phase 2, spliced .fna) are
located in each lncRNA by their last 30 nt (stop codon included) and confirmed by the full ORF sequence,
or, when the catalog start lies upstream of the lncRNA, by the overlapping part. A prediction 'matches'
when it is on the + strand and ends at a catalog ORF's stop (within TOL nt).
Reported per organism, variant and label:
  lnc_pred        lncRNAs with a prediction
  lnc_pred_cat    of those, predictions ending at a catalog ORF stop
  cat_orfs        catalog ORFs of >= MIN_CODONS codons located in the lncRNAs (one per lncRNA and stop)
  cat_found       of those, with a prediction ending at their stop

usage: score_lnc_riboseq_orfs.py <catalog.fna> <inputs_dir> <org> <variants> label=<template> ... > out.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, read_preds, read_truth  # noqa: E402

K = 30
MIN_CODONS = 50


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf).upper()


def main():
    cat_fna, inputs, org, variants = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(",")
    labels = [a.split("=", 1) for a in sys.argv[5:]]
    tail = {}
    for name, seq in read_fasta(cat_fna):
        if len(seq) >= K:
            tail.setdefault(seq[-K:], []).append(seq)
    print("org\tvariant\tlabel\tn_lnc\tlnc_pred\tlnc_pred_cat\tcat_orfs\tcat_found")
    for v in variants:
        truth = read_truth(os.path.join(inputs, org, v + ".truth.tsv"))
        stops, long_stops = {}, {}  # lncRNA -> catalog stop ends (1-based, last nt of stop)
        n_lnc = 0
        for tid, seq in read_fasta(os.path.join(inputs, org, v + ".fa")):
            if tid in truth:
                continue
            n_lnc += 1
            for i in range(len(seq) - K + 1):
                hits = tail.get(seq[i:i + K])
                if not hits:
                    continue
                end = i + K
                for orf in hits:
                    up = seq[max(0, end - len(orf)):end]
                    if orf.endswith(up):
                        stops.setdefault(tid, set()).add(end)
                        if len(orf) // 3 >= MIN_CODONS and len(up) == len(orf):
                            long_stops.setdefault(tid, set()).add(end)
        n_cat = sum(len(s) for s in long_stops.values())
        for label, tmpl in labels:
            pf = tmpl.format(org=org, variant=v)
            if not os.path.exists(pf):
                continue
            preds = read_preds(pf)
            lp = lpc = found = 0
            for tid, plist in preds.items():
                if tid in truth:
                    continue
                lp += 1
                ends_ = [p[1] for p in plist if p[2] == "+"]
                if any(abs(e - s) <= TOL for e in ends_ for s in stops.get(tid, ())):
                    lpc += 1
                found += sum(1 for s in long_stops.get(tid, ()) if any(abs(e - s) <= TOL for e in ends_))
            print(f"{org}\t{v}\t{label}\t{n_lnc}\t{lp}\t{lpc}\t{n_cat}\t{found}")


if __name__ == "__main__":
    main()
