#!/usr/bin/env python3
"""Trim predicted CDS coordinates that run past either end of the transcript back inside it, in
whole codons (r2012 reports 3'-partial ORFs as if a stop codon followed the transcript end).

usage: clip_preds_to_transcript.py <transcripts.fa> <preds.tsv> > clipped.preds.tsv
"""
import sys


def lengths(fasta):
    out, name = {}, None
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            else:
                out[name] += len(line.strip())
    return out


def main():
    tlen = lengths(sys.argv[1])
    n = 0
    with open(sys.argv[2]) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "transcript_id":
                print(line, end="")
                continue
            s, e, L = int(f[1]), int(f[2]), tlen[f[0]]
            if e > L or s < 1:
                n += 1
            while e > L:
                e -= 3
            while s < 1:
                s += 3
            f[1], f[2] = str(s), str(e)
            print("\t".join(f))
    print(f"clipped {n} predictions", file=sys.stderr)


if __name__ == "__main__":
    main()
