#!/usr/bin/env python3
"""Reference sensitivity restricted to genes expressed in the sample.

usage: expressed_sensitivity.py <stringtie -e ref_quant.gtf> <label=ref_status.tsv> [...]
       [--tpm 1,5,10]

Expression: StringTie -e (reference transcripts only) on the sample's RNA-seq alignments. Gene TPM is
the sum of its transcripts' TPMs. A reference transcript counts as expressed when its gene's TPM is
at least the threshold. Sensitivities are computed from eval_genome_cds.py's ref_status.tsv.

Output (stdout, TSV): label, set, tpm_min, n_expressed, assembled_frac, sens_exact, sens_compatible,
sens_exact_given_assembled
"""
import argparse
import re
import sys
from collections import defaultdict

ATTR = re.compile(r'(\S+) "([^"]*)"')


def read_quant(path):
    tx_gene, gene_tpm = {}, defaultdict(float)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        a = dict(ATTR.findall(f[8]))
        tid, gid = a.get("reference_id", a["transcript_id"]), a.get("ref_gene_id", a["gene_id"])
        tx_gene[tid.split(".")[0]] = gid
        gene_tpm[gid] += float(a.get("TPM", 0))
    return tx_gene, gene_tpm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("quant_gtf")
    ap.add_argument("ref_status", nargs="+", help="label=path/to/ref_status.tsv")
    ap.add_argument("--tpm", default="1,5,10")
    args = ap.parse_args()
    tx_gene, gene_tpm = read_quant(args.quant_gtf)
    thresholds = [float(t) for t in args.tpm.split(",")]
    print("label\tset\ttpm_min\tn\tassembled_frac\tsens_exact\tsens_compatible\tsens_exact_given_assembled")
    for spec in args.ref_status:
        label, path = spec.split("=", 1)
        rows, missing = [], 0
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if f[0] == "set":
                continue
            gid = tx_gene.get(f[1].split(".")[0])
            if gid is None:
                missing += 1
                continue
            rows.append((f[0], gene_tpm[gid], int(f[3]), int(f[4]), int(f[5])))
        if missing:
            print(f"[{label}] {missing} reference transcripts not in the quantification", file=sys.stderr)
        for setname in ("MANE", "GENCODE_complete"):
            for t in [0.0] + thresholds:
                sel = [r for r in rows if r[0] == setname and r[1] >= t]
                n = len(sel)
                asm = sum(r[2] for r in sel)
                ex = sum(r[3] for r in sel)
                comp = sum(r[4] for r in sel)
                ex_asm = sum(r[2] and r[3] for r in sel)
                print(f"{label}\t{setname}\t{t:g}\t{n}\t{asm / max(n, 1):.4f}\t{ex / max(n, 1):.4f}\t"
                      f"{comp / max(n, 1):.4f}\t{ex_asm / max(asm, 1):.4f}")


if __name__ == "__main__":
    main()
