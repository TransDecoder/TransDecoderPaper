#!/usr/bin/env python3
"""Benchmark accuracy on fragments, split by partial type of the reference CDS.

Fragments (methods_benchmark.tex) are cut once inside the CDS: the 3' piece keeps the stop but not
the start (5prime_partial), the 5' piece keeps the start but not the stop (3prime_partial). For a
3' partial the reference 3' end is the last whole codon at the fragment end, so a 'stop-level'
match means the predicted CDS runs to the fragment end. Matching is that of
benchmark_v2/scripts/score_predictions.py (ends within 2 nt, same strand, one credit per reference).

Per (org, variant, tool, type): references, stop- and exact-level true positives, predictions on
transcripts of that type and how many of them are credited, giving sensitivity and precision
(precision excludes predictions on random negatives, which have no type).

usage: 02_partial_types.py <inputs_dir> <runs_dir> > data/benchmark/final_v7/partial_types.tsv
"""
import collections
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark_v2", "scripts"))
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

VARIANTS = ["fragments", "fragments.unstranded", "fragments.wRand", "fragments.unstranded.wRand"]
TOOLS = ["pytd_v7", "perl_td_v7", "perl_td_v600", "prodigal", "gmst", "codan", "cpat", "cpc2", "rnasamba", "borf",
         "longest_orf"]


def score(truth, preds):
    c = collections.defaultdict(collections.Counter)
    for tid, t in truth.items():
        typ = t[3]
        c[typ]["n_truth"] += 1
        plist = preds.get(tid, [])
        c[typ]["n_pred"] += len(plist)
        best = None
        for p in plist:
            if p[2] != t[2]:
                continue
            p5, p3 = ends(*p)
            t5, t3 = ends(*t[:3])
            if abs(p3 - t3) <= TOL:
                lvl = "exact" if abs(p5 - t5) <= TOL else "stop"
                if best != "exact":
                    best = lvl
        if best:
            c[typ]["tp_stop"] += 1
            c[typ]["tp_exact"] += best == "exact"
    return c


def main():
    inputs, runs = sys.argv[1], sys.argv[2]
    cols = ["org", "variant", "tool", "type", "n_truth", "n_pred", "tp_stop", "tp_exact",
            "sensitivity_stop", "sensitivity_exact", "precision_stop", "precision_exact"]
    print("\t".join(cols))
    for org in sorted(os.listdir(runs)):
        for variant in VARIANTS:
            tf = os.path.join(inputs, org, variant + ".truth.tsv")
            if not os.path.exists(tf):
                continue
            truth = read_truth(tf)
            for tool in TOOLS:
                if tool == "longest_orf":
                    pf = os.path.join(runs, org, variant, "pytd", "longest_orf.preds.tsv")
                else:
                    pf = os.path.join(runs, org, variant, tool, "preds.tsv")
                if not os.path.exists(pf):
                    continue
                for typ, k in sorted(score(truth, read_preds(pf)).items()):
                    n, npred = k["n_truth"], k["n_pred"]
                    row = [org, variant, tool, typ, n, npred, k["tp_stop"], k["tp_exact"],
                           k["tp_stop"] / n, k["tp_exact"] / n,
                           k["tp_stop"] / npred if npred else 0, k["tp_exact"] / npred if npred else 0]
                    print("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in row))


if __name__ == "__main__":
    main()
