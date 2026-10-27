#!/usr/bin/env python3
"""Small summary tables of the simulated GM12878 Trinity benchmark for the notebook (evidence E30).

usage: 09_summarize.py <sim_labels_prefix> <real_labels_prefix> <counts.tsv> <sim_e29.pairs.tsv.gz>
                       <real_e29.pairs.tsv.gz> <outdir>
  *_labels_prefix   label_assembly.py outputs (<prefix>.labels.tsv) for the simulated and the real
                    GM12878 Trinity GG assemblies (the real one labelled against the same transcript set)
  counts.tsv        03_make_counts.py output (simulated read pairs and salmon TPM per kept transcript)
  *_e29.pairs.tsv.gz  pipeline/derive/03_utr_gm12878_truth.py outputs (E29 statistic)
Writes to <outdir>:
  label_counts.tsv        assembly x label x reason: transcripts
  type_by_expr.tsv        simulated assembly, scored coding contigs: TPM bin x type (contigs, fraction)
  complete_by_source_expr.tsv  simulated coding transcripts with >= 1 read pair, by TPM bin: how many
                          have >= 1 assembled contig labelled complete with them as a source, and how
                          many have any labelled (scored) contig
  e29_five_prime.tsv      assembly: contigs carrying a MANE Select stop with reference structure to their
                          5' end, how many begin inside the CDS (E29 statistic)
"""
import collections
import gzip
import os
import sys

EXPR = [(0, 1, "<1"), (1, 5, "1-5"), (5, 20, "5-20"), (20, 100, "20-100"), (100, float("inf"), ">100")]
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]


def ebin(x):
    return next(lab for lo, hi, lab in EXPR if lo <= x < hi)


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]


def main():
    simp, realp, counts_f, sim_e29, real_e29, outdir = sys.argv[1:7]
    os.makedirs(outdir, exist_ok=True)
    labs = {"simulated": read_tsv(simp + ".labels.tsv"), "real": read_tsv(realp + ".labels.tsv")}
    with open(os.path.join(outdir, "label_counts.tsv"), "w") as fh:
        fh.write("assembly\tlabel\treason\tn\n")
        for asm, rows in labs.items():
            c = collections.Counter((r["label"], r["reason"] if r["label"] in ("unresolved", "ambiguous") else "")
                                    for r in rows)
            for (l, r), n in sorted(c.items()):
                fh.write(f"{asm}\t{l}\t{r or '.'}\t{n}\n")

    sim = labs["simulated"]
    c = collections.Counter()
    for r in sim:
        if r["label"] in TYPES:
            c[(ebin(float(r["tpm"])), r["label"])] += 1
    with open(os.path.join(outdir, "type_by_expr.tsv"), "w") as fh:
        fh.write("tpm_bin\ttype\tn\tfrac\n")
        for _, _, b in EXPR:
            tot = sum(c[(b, t)] for t in TYPES)
            for t in TYPES:
                fh.write(f"{b}\t{t}\t{c[(b, t)]}\t{c[(b, t)] / tot if tot else 0:.4f}\n")

    comp_src, any_src = set(), set()
    for r in sim:
        if r["label"] in TYPES:
            for s in r["source_transcripts"].split(","):
                any_src.add(s)
                if r["label"] == "complete":
                    comp_src.add(s)
    cc = collections.Counter()
    for r in read_tsv(counts_f):
        if r["class"] != "coding" or int(r["count"]) == 0:
            continue
        b = ebin(float(r["tpm"]))
        cc[(b, "n")] += 1
        cc[(b, "complete")] += r["transcript_id"] in comp_src
        cc[(b, "any")] += r["transcript_id"] in any_src
    with open(os.path.join(outdir, "complete_by_source_expr.tsv"), "w") as fh:
        fh.write("tpm_bin\tn_coding_simulated\twith_complete_contig\tfrac_complete\twith_scored_contig\tfrac_scored\n")
        for _, _, b in EXPR:
            n = cc[(b, "n")]
            fh.write(f"{b}\t{n}\t{cc[(b, 'complete')]}\t{cc[(b, 'complete')] / n:.4f}\t{cc[(b, 'any')]}\t{cc[(b, 'any')] / n:.4f}\n")

    with open(os.path.join(outdir, "e29_five_prime.tsv"), "w") as fh:
        fh.write("assembly\tn_with_mane_stop\tn_5prime_partial\tfrac_5prime_partial\tmissing_frac_median\n")
        for asm, f in (("simulated", sim_e29), ("real", real_e29)):
            rows = read_tsv(f)
            p = sorted(float(r["missing_frac"]) for r in rows if r["status"] == "5prime_partial")
            med = p[len(p) // 2] if len(p) % 2 else (p[len(p) // 2 - 1] + p[len(p) // 2]) / 2
            fh.write(f"{asm}\t{len(rows)}\t{len(p)}\t{len(p) / len(rows):.4f}\t{med:.4f}\n")


if __name__ == "__main__":
    main()
