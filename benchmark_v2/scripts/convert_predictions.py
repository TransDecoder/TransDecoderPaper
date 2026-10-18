#!/usr/bin/env python3
"""Convert predictor outputs to the standard prediction TSV used by score_predictions.py.

usage: convert_predictions.py <format> [options] <file>

formats:
  transdecoder_gff3   TransDecoder / PyTransDecoder .gff3 (CDS features; also longest_orfs.gff3)
  prodigal_gff        Prodigal -f gff (coordinates include the stop codon)
  gmst_gff            GeneMarkS-T --format GFF (same layout as Prodigal; coordinates include the stop)
  codan_gtf           CodAn annotation.gtf (CDS excludes the stop; stop_codon rows are merged in)

Output columns: transcript_id, start, end, strand, score (1-based inclusive, transcript
coordinates; CDS including the stop codon when the tool reports one).
"""
import argparse
import collections
import re
import sys


def transdecoder(path):
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if len(x) < 9 or x[2] != "CDS":
            continue
        m = re.search(r"Parent=([^;]+)", x[8])
        yield x[0], int(x[3]), int(x[4]), x[6], ".", (m.group(1) if m else None)


def prodigal(path, plus_only):
    for line in open(path):
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        if len(x) < 9 or x[2] != "CDS" or (plus_only and x[6] != "+"):
            continue
        yield x[0], int(x[3]), int(x[4]), x[6], x[5], None


def codan(path):
    feats = collections.OrderedDict()
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if len(x) < 9 or x[2] not in ("CDS", "stop_codon"):
            continue
        tid = re.search(r'transcript_id "([^"]+)"', x[8]).group(1)
        key = (x[0], tid)
        s, e = int(x[3]), int(x[4])
        if key in feats:
            ps, pe, strand, score = feats[key]
            feats[key] = (min(ps, s), max(pe, e), strand, score)
        else:
            feats[key] = (s, e, x[6], x[5].strip())
    for (seqid, _), (s, e, strand, score) in feats.items():
        yield seqid, s, e, strand, score, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("format", choices=["transdecoder_gff3", "prodigal_gff", "gmst_gff", "codan_gtf"])
    ap.add_argument("file")
    ap.add_argument("--plus-only", action="store_true", help="prodigal: keep + strand predictions only")
    ap.add_argument("--longest-per-transcript", action="store_true",
                    help="keep only the longest prediction per transcript (ties: first listed)")
    args = ap.parse_args()

    if args.format == "transdecoder_gff3":
        rows = transdecoder(args.file)
    elif args.format in ("prodigal_gff", "gmst_gff"):
        rows = prodigal(args.file, args.plus_only)
    else:
        rows = codan(args.file)

    rows = list(rows)
    if args.longest_per_transcript:
        best = collections.OrderedDict()
        for r in rows:
            if r[0] not in best or (r[2] - r[1]) > (best[r[0]][2] - best[r[0]][1]):
                best[r[0]] = r
        rows = list(best.values())

    out = sys.stdout
    out.write("transcript_id\tstart\tend\tstrand\tscore\n")
    for r in rows:
        out.write("\t".join(map(str, r[:5])) + "\n")


if __name__ == "__main__":
    sys.exit(main())
