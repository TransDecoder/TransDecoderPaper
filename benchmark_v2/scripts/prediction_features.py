#!/usr/bin/env python3
"""Per-prediction coding evidence for PyTD runs, labelled against the truth, to design a precise mode.

For every final ORF of a PyTD run (single best per transcript): ORF length (codons, LongOrfs ORF),
in-frame score, frame contrast (in-frame score minus the better shifted frame), the training ORFs'
median per-base contrast mu, and the call: TP_exact, TP_stop, wrong (coding transcript, other ORF),
neg (random or lncRNA sequence).

usage: prediction_features.py <inputs_dir> <run_root> <setting> <variants> > out.tsv
  run_root/<org>/<variant>/<setting>/ holds input.fa.transdecoder.gff3, predict.log and
  input.fa.transdecoder_dir/longest_orfs.cds.scores
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, ends, read_truth  # noqa: E402

ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
MU = re.compile(r"median frame contrast ([0-9.eE+-]+) per base")


def main():
    inputs, root, setting, variants = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(",")
    print("org\tvariant\ttranscript_id\torf_id\tlength_codons\tscore\tcontrast\tmu\tcall")
    for v in variants:
        for org in ORGS:
            d = os.path.join(root, org, v, setting)
            gff = os.path.join(d, "input.fa.transdecoder.gff3")
            if not os.path.exists(gff):
                continue
            truth = read_truth(os.path.join(inputs, org, v + ".truth.tsv"))
            mu = None
            with open(os.path.join(d, "predict.log")) as fh:
                for line in fh:
                    m = MU.search(line)
                    if m:
                        mu = float(m.group(1))
            scores = {}
            with open(os.path.join(d, "input.fa.transdecoder_dir", "longest_orfs.cds.scores")) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split("\t")
                    s = [float(x) for x in f[3:6]]
                    scores[f[0]] = (int(f[2]) // 3, s[0], s[0] - max(s[1], s[2]))
            with open(gff) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 9 or f[2] != "CDS":
                        continue
                    orf = re.search(r"Parent=([^;]+)", f[8]).group(1)
                    tid, s, e, strand = f[0], int(f[3]), int(f[4]), f[6]
                    length, sc, con = scores[orf]
                    t = truth.get(tid)
                    if not t:
                        call = "neg"
                    elif strand != t[2]:
                        call = "wrong"
                    else:
                        p5, p3 = ends(s, e, strand)
                        t5, t3 = ends(*t[:3])
                        call = ("TP_exact" if abs(p5 - t5) <= TOL else "TP_stop") if abs(p3 - t3) <= TOL else "wrong"
                    print(f"{org}\t{v}\t{tid}\t{orf}\t{length}\t{sc:.2f}\t{con:.2f}\t{mu}\t{call}")


if __name__ == "__main__":
    main()
