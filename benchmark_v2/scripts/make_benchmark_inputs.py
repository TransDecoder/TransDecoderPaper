#!/usr/bin/env python3
"""Derive benchmark input sets from a truth set (build_truth_from_gbff.py output).

For each truth set this writes, under --outdir:

  complete.fa / complete.truth.tsv
      the reference transcripts as given (CDS on the + strand)
  fragments.fa / fragments.truth.tsv
      each transcript cut once inside its CDS at a fraction drawn from N(0.5, 0.1)
      (clipped to [0.05, 0.95]), as in the 2018 benchmark. If the cut falls 3' of the
      CDS midpoint the 3' piece is kept (stop present, start lost; type 5prime_partial),
      otherwise the 5' piece (start present, stop lost; type 3prime_partial). Partial
      CDS boundaries are trimmed to whole codons in the original frame.
  <set>.unstranded.fa / .truth.tsv
      the same set with a random half of the transcripts reverse-complemented
      (truth strand '-'), for evaluating strand-agnostic prediction
  <set>[.unstranded].wRand.fa / .truth.tsv
      the set plus one random negative per transcript ('rand_<id>'): i.i.d. sequence
      of the same length and the same G+C fraction as that transcript. Negatives
      have no truth rows; any prediction on them is a false positive.

All randomness comes from --seed.
"""
import argparse
import os
import random
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, acc, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def read_truth(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return hdr, rows


def write_set(outdir, name, hdr, seqs, rows):
    with open(os.path.join(outdir, name + ".fa"), "w") as fa:
        for acc, s in seqs.items():
            fa.write(f">{acc}\n{s}\n")
    with open(os.path.join(outdir, name + ".truth.tsv"), "w") as tt:
        tt.write("\t".join(hdr) + "\n")
        for r in rows:
            tt.write("\t".join(str(r[h]) for h in hdr) + "\n")


def fragment(seqs, rows, rng):
    out_seqs, out_rows = {}, []
    for r in rows:
        acc = r["transcript_id"]
        seq = seqs[acc]
        cs, ce = int(r["cds_start"]), int(r["cds_end"])  # 1-based inclusive
        frac = min(0.95, max(0.05, rng.gauss(0.5, 0.1)))
        cut = cs - 1 + int((ce - cs + 1) * frac)  # 0-based cut point inside the CDS
        r2 = dict(r)
        if cut - (cs - 1) > (ce - cs + 1) / 2:
            # keep the 3' piece: seq[cut:]; first whole codon in the original frame
            piece = seq[cut:]
            first = cut + (-(cut - (cs - 1))) % 3  # 0-based codon start >= cut
            r2.update(cds_start=first - cut + 1, cds_end=ce - cut, type="5prime_partial")
        else:
            piece = seq[:cut]
            ncod = (cut - (cs - 1)) // 3
            r2.update(cds_start=cs, cds_end=cs - 1 + 3 * ncod, type="3prime_partial")
        r2["transcript_length"] = len(piece)
        out_seqs[acc] = piece
        out_rows.append(r2)
    return out_seqs, out_rows


def unstrand(seqs, rows, rng):
    out_seqs, out_rows = {}, []
    for r in rows:
        acc = r["transcript_id"]
        s = seqs[acc]
        r2 = dict(r)
        if rng.random() < 0.5:
            L = len(s)
            cs, ce = int(r["cds_start"]), int(r["cds_end"])
            s = revcomp(s)
            r2.update(cds_start=L - ce + 1, cds_end=L - cs + 1, strand="-")
        out_seqs[acc] = s
        out_rows.append(r2)
    return out_seqs, out_rows


def with_random(seqs, rng):
    out = dict(seqs)
    for acc, s in seqs.items():
        gc = (s.count("G") + s.count("C")) / max(1, len(s))
        w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        out["rand_" + acc] = "".join(rng.choices("ACGT", weights=w, k=len(s)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    seqs = read_fasta(args.fasta)
    hdr, rows = read_truth(args.truth)
    sets = {"complete": (seqs, rows)}
    sets["fragments"] = fragment(seqs, rows, random.Random(args.seed))
    for name in list(sets):
        sets[name + ".unstranded"] = unstrand(*sets[name], random.Random(args.seed + 1))
    for i, name in enumerate(sorted(sets)):
        s, r = sets[name]
        write_set(args.outdir, name, hdr, s, r)
        write_set(args.outdir, name + ".wRand", hdr, with_random(s, random.Random(args.seed + 10 + i)), r)


if __name__ == "__main__":
    sys.exit(main())
