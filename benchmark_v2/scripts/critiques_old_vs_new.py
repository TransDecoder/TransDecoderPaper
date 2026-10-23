#!/usr/bin/env python3
"""Measure, for an old and a new TransDecoder, the failure modes reported in earlier benchmarks.

Per organism and input, for each tool label given:
  n_truth, n_pred
  multi_orf_tx     transcripts with more than one prediction (GeneMarkS-T 2015 Table S2; inGAP-CDG 2016)
  non_atg_start    predictions on coding transcripts whose first codon is not ATG, i.e. 5'-open
                   (TranSuite 2020; Borf 2021)
  start_upstream   stop-level TPs whose start lies upstream of the annotated start
                   (GeneMarkS-T 2015 Table 2 'longer'; ORFanage 2023)
  start_upstream_open  ... of which the prediction runs to the 5' end of the transcript (no ATG)
  wrong_strand     predictions on the opposite strand of a coding transcript (Borf 2021; unstranded inputs)
  no_pred          coding transcripts with no prediction (ORFanage 2023)
  neg_with_pred    negatives (random 'rand_' or lncRNA, i.e. sequences without truth) with a prediction
  n_neg            negatives in the input
  exact, stop      exact and stop-level TPs (as score_predictions.py)

usage: critiques_old_vs_new.py <inputs_dir> <runs_dir> <variants> <label=subdir/file> ... > out.tsv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
RC = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta(path):
    seqs, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def first_codon(seq, s, e, strand):
    cds = seq[s - 1:e] if strand == "+" else seq[s - 1:e].translate(RC)[::-1]
    return cds[:3]


def main():
    inputs, runs, variants = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    labels = [a.split("=", 1) for a in sys.argv[4:]]
    cols = ["n_truth", "n_pred", "multi_orf_tx", "non_atg_start", "start_upstream", "start_upstream_open",
            "wrong_strand", "no_pred", "neg_with_pred", "n_neg", "exact", "stop"]
    print("org\tvariant\ttool\t" + "\t".join(cols))
    for v in variants:
        for org in ORGS:
            fa = os.path.join(inputs, org, v + ".fa")
            tf = os.path.join(inputs, org, v + ".truth.tsv")
            if not os.path.exists(fa):
                continue
            seqs, truth = read_fasta(fa), read_truth(tf)
            n_neg = sum(1 for k in seqs if k not in truth)
            for label, rel in labels:
                pf = os.path.join(runs, org, v, rel)
                if not os.path.exists(pf):
                    continue
                preds = read_preds(pf)
                c = dict.fromkeys(cols, 0)
                c["n_truth"], c["n_neg"] = len(truth), n_neg
                for tid, plist in preds.items():
                    t = truth.get(tid)
                    c["n_pred"] += len(plist)
                    if len(plist) > 1:
                        c["multi_orf_tx"] += 1
                    if not t:
                        c["neg_with_pred"] += 1
                        continue
                    best = None
                    for p in plist:
                        if first_codon(seqs[tid], *p) != "ATG":
                            c["non_atg_start"] += 1
                        if p[2] != t[2]:
                            c["wrong_strand"] += 1
                            continue
                        p5, p3 = ends(*p)
                        t5, t3 = ends(*t[:3])
                        if abs(p3 - t3) <= TOL:
                            if abs(p5 - t5) <= TOL:
                                best = ("exact", p)
                            elif best is None or best[0] != "exact":
                                best = ("stop", p)
                    if best:
                        c["stop"] += 1
                        if best[0] == "exact":
                            c["exact"] += 1
                        else:
                            p = best[1]
                            p5 = ends(*p)[0]
                            t5 = ends(*t[:3])[0]
                            up = p5 < t5 if t[2] == "+" else p5 > t5
                            if up:
                                c["start_upstream"] += 1
                                if first_codon(seqs[tid], *p) != "ATG":
                                    c["start_upstream_open"] += 1
                for tid in truth:
                    if tid not in preds:
                        c["no_pred"] += 1
                print(f"{org}\t{v}\t{label}\t" + "\t".join(str(c[k]) for k in cols))


if __name__ == "__main__":
    main()
