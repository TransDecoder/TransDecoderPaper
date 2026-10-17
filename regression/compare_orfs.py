#!/usr/bin/env python3
"""Compare ORF predictions between two TransDecoder-style GFF3 files.

ORFs are keyed by transcript, strand and CDS span. Exact matches need identical
CDS coordinates; stop-matches share the strand and the stop-codon end (i.e. they
differ only in the chosen start).

usage: compare_orfs.py A.gff3 B.gff3 [--label-a perl --label-b py] [--tsv out.tsv]
"""
import argparse
import sys
from collections import defaultdict


def read_cds(path):
    """Return {(transcript, strand, lend, rend)} for CDS features, merging per-parent."""
    spans = defaultdict(lambda: [None, None, None, None])
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            parent = attrs.get("Parent", attrs.get("ID"))
            lend, rend = int(f[3]), int(f[4])
            s = spans[parent]
            s[0], s[1] = f[0], f[6]
            s[2] = lend if s[2] is None else min(s[2], lend)
            s[3] = rend if s[3] is None else max(s[3], rend)
    return {tuple(v) for v in spans.values()}


def stop_key(orf):
    tx, strand, lend, rend = orf
    return (tx, strand, rend if strand == "+" else lend)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--tsv", help="write per-ORF discordance table")
    args = ap.parse_args()

    a, b = read_cds(args.a), read_cds(args.b)
    exact = a & b
    a_only, b_only = a - b, b - a
    a_stop = {stop_key(o): o for o in a_only}
    b_stop = {stop_key(o): o for o in b_only}
    start_diff = set(a_stop) & set(b_stop)

    la, lb = args.label_a, args.label_b
    print(f"{la}_orfs\t{len(a)}")
    print(f"{lb}_orfs\t{len(b)}")
    print(f"exact_match\t{len(exact)}")
    print(f"same_stop_diff_start\t{len(start_diff)}")
    print(f"{la}_only\t{len(a_only) - len(start_diff)}")
    print(f"{lb}_only\t{len(b_only) - len(start_diff)}")
    union = len(a | b) - len(start_diff)
    print(f"jaccard_exact\t{len(exact) / max(1, len(a | b)):.4f}")
    print(f"jaccard_stop\t{(len(exact) + len(start_diff)) / max(1, union):.4f}")

    if args.tsv:
        with open(args.tsv, "w") as out:
            out.write("category\ttranscript\tstrand\tlend_a\trend_a\tlend_b\trend_b\n")
            for k in sorted(start_diff):
                oa, ob = a_stop[k], b_stop[k]
                out.write(f"diff_start\t{oa[0]}\t{oa[1]}\t{oa[2]}\t{oa[3]}\t{ob[2]}\t{ob[3]}\n")
            for k, o in sorted(a_stop.items()):
                if k not in start_diff:
                    out.write(f"{la}_only\t{o[0]}\t{o[1]}\t{o[2]}\t{o[3]}\t.\t.\n")
            for k, o in sorted(b_stop.items()):
                if k not in start_diff:
                    out.write(f"{lb}_only\t{o[0]}\t{o[1]}\t.\t.\t{o[2]}\t{o[3]}\n")


if __name__ == "__main__":
    sys.exit(main())
