#!/usr/bin/env python3
"""Simulate the length distribution of ORFs that LongOrfs extracts from random sequence.

TransDecoder.Predict automatically retains any candidate ORF whose CDS length meets a
GC-dependent threshold. That threshold is the 0.999 quantile of the lengths of ORFs
found in random sequence of the same GC content. This script regenerates the table.

Random sequences are i.i.d. with P(G)=P(C)=gc/2 and P(A)=P(T)=(1-gc)/2. ORFs follow
LongOrfs: for each stop codon, take the most upstream in-frame ATG that follows the
previous in-frame stop. Both strands are scanned. With --partials, an ORF may also
start at the sequence 5' end (no upstream stop in frame) and/or run off the 3' end.
Lengths are CDS nucleotides, including the stop codon, and only ORFs >= --min-len-nt
are kept (LongOrfs default -m 100 codons = 300 nt).

Output: a TSV of quantiles per GC level, plus optionally the raw lengths.
"""
import argparse
import gzip
import random
import sys

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGT", "TGCA")


def orf_lengths(seq, partials, min_len):
    lens = []
    for s in (seq, seq.translate(COMP)[::-1]):
        L = len(s)
        for frame in range(3):
            start = None  # first ATG (or 5' partial start) since the last stop
            if partials:
                start = frame  # 5' partial: open frame at the sequence start
            last_codon = frame
            for i in range(frame, L - 2, 3):
                codon = s[i:i + 3]
                last_codon = i
                if codon in STOPS:
                    if start is not None:
                        n = i + 3 - start
                        if n >= min_len:
                            lens.append(n)
                    start = None
                elif start is None and codon == "ATG":
                    start = i
            if partials and start is not None:
                n = last_codon + 3 - start
                if n >= min_len:
                    lens.append(n)
    return lens


def quantile(sorted_vals, q):
    """Linear-interpolation quantile (R type 7, numpy default)."""
    if not sorted_vals:
        return float("nan")
    h = (len(sorted_vals) - 1) * q
    lo = int(h)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (h - lo) * (sorted_vals[hi] - sorted_vals[lo])


def read_lengths(path):
    op = gzip.open if path.endswith(".gz") else open
    lens, n = [], 0
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if n:
                    lens.append(n)
                n = 0
            else:
                n += len(line.strip())
    if n:
        lens.append(n)
    return lens


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gc", default="25,30,35,40,45,50,55,60,65,70,75,80")
    ap.add_argument("--n-seqs", type=int, default=10000)
    ap.add_argument("--seq-len", type=int, default=2000, help="fixed sequence length (ignored with --length-fasta)")
    ap.add_argument("--length-fasta", help="draw sequence lengths from this FASTA's length distribution")
    ap.add_argument("--partials", action="store_true", help="include 5'/3' partial ORFs (LongOrfs default)")
    ap.add_argument("--min-len-nt", type=int, default=300)
    ap.add_argument("--quantiles", default="0.95,0.99,0.999,0.9999")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--raw-out", help="write gc<TAB>length for every ORF")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    qs = [float(q) for q in args.quantiles.split(",")]
    pool = read_lengths(args.length_fasta) if args.length_fasta else None
    raw = open(args.raw_out, "w") if args.raw_out else None

    print("gc\tn_orfs\t" + "\t".join(f"q{q}" for q in qs))
    for gc in (int(g) for g in args.gc.split(",")):
        p = gc / 200.0
        weights = [0.5 - p, p, p, 0.5 - p]  # A C G T
        lens = []
        for _ in range(args.n_seqs):
            L = rng.choice(pool) if pool else args.seq_len
            seq = "".join(rng.choices("ACGT", weights=weights, k=L))
            lens.extend(orf_lengths(seq, args.partials, args.min_len_nt))
        lens.sort()
        if raw:
            raw.writelines(f"{gc}\t{x}\n" for x in lens)
        print(f"{gc}\t{len(lens)}\t" + "\t".join(f"{quantile(lens, q):.1f}" for q in qs), flush=True)


if __name__ == "__main__":
    sys.exit(main())
