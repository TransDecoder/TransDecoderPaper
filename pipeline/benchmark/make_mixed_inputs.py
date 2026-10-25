#!/usr/bin/env python3
"""Make the balanced mixed benchmark inputs fragments_mixed_t1, _t2, _t3 (E28).

The file names start with 'fragments' on purpose: run_tool.sh selects CodAn's *_partial models and
RNAsamba's partial-length weights for inputs whose name starts with 'fragments', which is the
appropriate setting for an input that contains partial CDSs (75% here).

Each input mixes complete transcripts and the three kinds of partial CDS in equal proportions, so
self-training methods (TransDecoder's coding and start models, GeneMarkS-T, Prodigal) see a
realistic mixture rather than one partial type only. Three independent trials (seeds 1, 2, 3)
give the variability due to the random assignment and cuts.

Source: the `complete` input of each organism (reference transcripts with UTRs <= 250 nt; the
transcripts `fragments` is cut from), <inputs>/<org>/complete.{fa,truth.tsv}. Stranded (all on
the + strand), no negatives.

Per trial t (seed t) and organism, with rng = random.Random(t):
  1. Category assignment: the transcripts (truth-file order) are shuffled with rng and assigned in
     turn to complete, 5prime_partial, 3prime_partial, internal (position mod 4), i.e. 25% each
     (counts differ by at most one).
  2. Cuts are then drawn with the same rng, in truth-file order. For a CDS at 1-based cs..ce
     (incl. ATG and stop), N = (ce - cs + 1) / 3 codons (codon 1 = ATG, codon N = stop) and
     MIN_CDS = 150 nt (50 codons) of CDS must remain:
     - complete: the transcript unchanged (no draw, never dropped).
     - 5prime_partial: cut codon k ~ U{2, ..., K}, K = min(floor(0.8 N), N - 50); phase
       p ~ U{0, 1, 2}; the fragment starts at base p of codon k (everything 5' of it, incl. the
       5' UTR and the ATG, removed) and keeps the rest of the CDS and the 3' UTR. Dropped if K < 2.
     - 3prime_partial: last whole kept codon m ~ U{max(ceil(0.2 N), 50), ..., N - 1}; extra bases
       q ~ U{0, 1, 2} of codon m + 1 also kept; the 5' UTR and codons 1..m are kept, everything 3'
       (incl. the stop codon and 3' UTR) removed. Dropped if max(ceil(0.2 N), 50) > N - 1.
     - internal: the allowed region is the CDS without its start and stop
       codons (M = 3N - 6 nt); window length W = round(f * 3N), f ~ N(0.5, 0.1) clipped to
       [0.05, 0.95], then clipped to [152, M]; offset ~ U{0, ..., M - W}. Both ends lie inside
       the CDS, no UTR remains; 152 nt in any phase contain >= 150 nt of whole codons. Dropped
       if M < 152.
  3. Truth, in fragment coordinates: cds_start = first whole codon of the original frame,
     cds_end = end of the stop codon (complete, 5' partial) or of the last whole codon (3'
     partial, internal); type as above; transcript_length; utr5_source / utr3_source 'none'
     where that UTR was removed. Extra columns:
       cut_codon_5p  codons of the CDS removed at the 5' end, wholly or in part (0 if none);
       cut_codon_3p  the same at the 3' end, counting the stop codon (0 if none);
       dist_first_inframe_atg  for 5' partials and internals: nt from the fragment start to the
                     first ATG in the CDS frame within the reference CDS of the fragment (0 = the
                     fragment starts with an in-frame ATG); NA if none, and NA for the other types.
Writes <outdir>/<org>/fragments_mixed_t<t>.{fa,truth.tsv} and prints per-organism counts to stdout.

usage: make_mixed_inputs.py <inputs_dir> [--outdir DIR] [--trials 1,2,3] [--orgs a,b]
"""
import argparse
import math
import os
import random
import sys

MIN_CDS = 150
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]


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
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]
    return hdr, rows


def cut(seq, cs, ce, typ, rng):
    """Return (fragment, cds_start, cds_end, cut5, cut3, dist_atg) or None (category impossible)."""
    N = (ce - cs + 1) // 3
    c0 = cs - 1                                        # 0-based first CDS base
    if typ == "complete":
        return seq, cs, ce, 0, 0, "NA"
    if typ == "5prime_partial":
        K = min(int(0.8 * N), N - MIN_CDS // 3)
        if K < 2:
            return None
        k, p = rng.randint(2, K), rng.randint(0, 2)
        s, e = c0 + 3 * (k - 1) + p, len(seq)
    elif typ == "3prime_partial":
        lo = max(math.ceil(0.2 * N), MIN_CDS // 3)
        if lo > N - 1:
            return None
        m, q = rng.randint(lo, N - 1), rng.randint(0, 2)
        s, e = 0, c0 + 3 * m + q
    else:  # internal
        M = 3 * N - 6
        f = min(0.95, max(0.05, rng.gauss(0.5, 0.1)))
        if M < MIN_CDS + 2:
            return None
        W = min(M, max(MIN_CDS + 2, int(round(f * 3 * N))))
        s = c0 + 3 + rng.randint(0, M - W)
        e = s + W
    first = c0 if s <= c0 else s + (c0 - s) % 3        # first whole CDS codon (0-based) >= s
    last_end = min(ce, first + 3 * ((e - first) // 3))  # 0-based exclusive end of the last whole codon
    cut5 = (first - c0) // 3
    cut3 = N - (last_end - c0) // 3
    dist = "NA"
    if typ in ("5prime_partial", "internal"):
        for i in range(first, min(last_end, ce - 3), 3):
            if seq[i:i + 3] == "ATG":
                dist = i - s
                break
    return seq[s:e], first - s + 1, last_end - s, cut5, cut3, dist


def make_trial(seqs, hdr, rows, seed):
    rng = random.Random(seed)
    order = list(range(len(rows)))
    rng.shuffle(order)
    cat = {}
    for pos, i in enumerate(order):
        cat[i] = TYPES[pos % 4]
    out_seqs, out_rows, dropped = {}, [], {t: 0 for t in TYPES}
    for i, r in enumerate(rows):
        assert r["strand"] == "+", r
        typ = cat[i]
        res = cut(seqs[r["transcript_id"]], int(r["cds_start"]), int(r["cds_end"]), typ, rng)
        if res is None:
            dropped[typ] += 1
            continue
        frag, a, b, c5, c3, dist = res
        r2 = dict(r)
        r2.update(cds_start=a, cds_end=b, transcript_length=len(frag), type=typ,
                  cut_codon_5p=c5, cut_codon_3p=c3, dist_first_inframe_atg=dist)
        if typ in ("5prime_partial", "internal"):
            r2["utr5_source"] = "none"
        if typ in ("3prime_partial", "internal"):
            r2["utr3_source"] = "none"
        out_seqs[r["transcript_id"]] = frag
        out_rows.append(r2)
    return out_seqs, out_rows, dropped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs")
    ap.add_argument("--outdir", help="default: the inputs directory")
    ap.add_argument("--trials", default="1,2,3", help="trial numbers; trial t uses seed t")
    ap.add_argument("--orgs", help="comma-separated organisms (default: all with complete.fa)")
    args = ap.parse_args()
    outdir = args.outdir or args.inputs
    orgs = args.orgs.split(",") if args.orgs else sorted(
        o for o in os.listdir(args.inputs) if os.path.exists(os.path.join(args.inputs, o, "complete.fa")))
    print("org\ttrial\tn_source\t" + "\t".join(f"n_{t}" for t in TYPES) + "\t"
          + "\t".join(f"dropped_{t}" for t in TYPES))
    for org in orgs:
        seqs = read_fasta(os.path.join(args.inputs, org, "complete.fa"))
        hdr, rows = read_truth(os.path.join(args.inputs, org, "complete.truth.tsv"))
        hdr2 = hdr + ["cut_codon_5p", "cut_codon_3p", "dist_first_inframe_atg"]
        for t in map(int, args.trials.split(",")):
            out_seqs, out_rows, dropped = make_trial(seqs, hdr, rows, t)
            name = f"fragments_mixed_t{t}"
            od = os.path.join(outdir, org)
            os.makedirs(od, exist_ok=True)
            with open(os.path.join(od, name + ".fa"), "w") as fa:
                for acc, s in out_seqs.items():
                    fa.write(f">{acc}\n{s}\n")
            with open(os.path.join(od, name + ".truth.tsv"), "w") as tt:
                tt.write("\t".join(hdr2) + "\n")
                for r in out_rows:
                    tt.write("\t".join(str(r[h]) for h in hdr2) + "\n")
            n = {ty: sum(r["type"] == ty for r in out_rows) for ty in TYPES}
            print(f"{org}\t{t}\t{len(rows)}\t" + "\t".join(str(n[ty]) for ty in TYPES) + "\t"
                  + "\t".join(str(dropped[ty]) for ty in TYPES))


if __name__ == "__main__":
    sys.exit(main())
