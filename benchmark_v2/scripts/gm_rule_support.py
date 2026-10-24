#!/usr/bin/env python3
"""Ribo-seq and MS support of GM12878 ORFs kept or dropped by a short-ORF rule, relative to a baseline run.

The rule runs use the same transcripts as the final_abinitio runs whose ORFs have Ribo-seq support
(riboseq_orf_support.py) and MS peptide tables (E01, E02). Each predicted ORF is matched to a
final_abinitio ORF by transcript, strand and stop codon (3' end), and inherits its support:
  ribo_n      ORFs with >= MIN_PSITES P-sites;  ribo_sup  of those, frame-supported (supported == 1)
  ms_sup      ORFs whose final_abinitio model has >= 1 MS peptide
Rows per assembly and label: all predictions, and kept / dropped / added relative to the baseline label
(dropped = baseline ORF without a prediction with the same stop in the label's run).

usage: gm_rule_support.py <gm12878_td dir> <baseline label> asm:label=preds.tsv ... > out.tsv
"""
import csv
import os
import re
import sys

MIN_PSITES = 10


def orf_stops(gff3):
    """final_abinitio ORF id -> (transcript, strand, 3' end) from CDS lines."""
    out = {}
    with open(gff3) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            orf = re.search(r"Parent=([^;]+)", f[8]).group(1)
            out[orf] = (f[0], f[6], int(f[4]) if f[6] == "+" else int(f[3]))
    return out


def pred_stops(path):
    s = set()
    with open(path) as fh:
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if x[0] == "transcript_id" or len(x) < 4:
                continue
            s.add((x[0], x[3], int(x[2]) if x[3] == "+" else int(x[1])))
    return s


def main():
    gm, base = sys.argv[1], sys.argv[2]
    runs = {}
    for a in sys.argv[3:]:
        key, path = a.split("=", 1)
        asm, label = key.split(":")
        runs.setdefault(asm, {})[label] = pred_stops(path)
    print("asm\tlabel\tset\tn\tmatched\tribo_n\tribo_sup\tribo_rate\tms_sup\tms_rate")
    for asm, labels in runs.items():
        stop_of = orf_stops(os.path.join(gm, asm, "final_abinitio", "transcripts.fa.transdecoder.gff3"))
        sup = {}
        with open(os.path.join(gm, "ribo_gm12878", f"gm12878.{asm}.orf_support.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["model_id"] in stop_of:
                    sup.setdefault(stop_of[r["model_id"]], {})["ribo"] = (int(r["psites"]), r["supported"] == "1")
        with open(os.path.join(gm, "proteomics", f"ms_peptide_support.{asm}.model_hits.tsv")) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["model_id"] in stop_of:
                    sup.setdefault(stop_of[r["model_id"]], {})["ms"] = int(r["n_peptides"]) > 0
        known = set(stop_of.values())

        def row(label, name, keys):
            m = [k for k in keys if k in known]
            rb = [sup[k]["ribo"][1] for k in m if "ribo" in sup.get(k, {}) and sup[k]["ribo"][0] >= MIN_PSITES]
            ms = sum(1 for k in m if sup.get(k, {}).get("ms"))
            print(f"{asm}\t{label}\t{name}\t{len(keys)}\t{len(m)}\t{len(rb)}\t{sum(rb)}\t"
                  f"{sum(rb) / len(rb) if rb else 0:.3f}\t{ms}\t{ms / len(m) if m else 0:.3f}")

        b = labels[base]
        for label, s in labels.items():
            row(label, "all", s)
            if label != base:
                row(label, "kept", s & b)
                row(label, "dropped", b - s)
                row(label, "added", s - b)


if __name__ == "__main__":
    main()
