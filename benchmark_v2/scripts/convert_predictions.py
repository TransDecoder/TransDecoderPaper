#!/usr/bin/env python3
"""Convert predictor outputs to the standard prediction TSV used by score_predictions.py.

usage: convert_predictions.py <format> [options] <file>

formats:
  transdecoder_gff3   TransDecoder / PyTransDecoder .gff3 (CDS features; also longest_orfs.gff3)
  prodigal_gff        Prodigal -f gff (coordinates include the stop codon)
  gmst_gff            GeneMarkS-T --format GFF (same layout as Prodigal; coordinates include the stop)
  codan_gtf           CodAn annotation.gtf (CDS excludes the stop; stop_codon rows are merged in)
  cpat_tsv            CPAT *.ORF_prob.best.tsv; best ORF kept if Coding_prob >= --cutoff
  cpc2_txt            CPC2 output run with --ORF; transcripts labelled coding. CPC2 reports no
                      strand, so with --fasta the ORF is placed on the strand where it translates
  borf_txt            Borf .txt (minus-strand ORFs are given on the reverse complement)
  rnasamba            RNAsamba --protein_fasta output (sequences classified coding); each protein is
                      located in its transcript (forward strand, standard code), needs --fasta
  angel_cds           ANGEL *.ANGEL.cds headers (pos: forward-strand coordinates including the stop)

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


COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CODE = {a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
        for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def translate(seq):
    return "".join(CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def id_mapper(fasta):
    """Map tool-reported IDs back to input IDs (CPAT upper-cases them)."""
    if not fasta:
        return lambda x: x
    upper = {k.upper(): k for k in read_fasta(fasta)}
    return lambda x: upper.get(x.upper(), x)


def cpat(path, cutoff, fasta):
    fix = id_mapper(fasta)
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if x[0] == "seq_ID" or float(x[10]) < cutoff:
            continue
        a, b = int(x[5]), int(x[6])
        yield fix(x[0]), min(a, b), max(a, b), x[3], x[10], None


def cpc2(path, fasta):
    seqs = read_fasta(fasta) if fasta else None
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if line.startswith("#") or x[8] != "coding":
            continue
        tid, tlen, plen, start = x[0], int(x[1]), int(x[2]), int(x[6])
        end = start + 3 * plen - 1          # peptide_length counts the stop when there is one
        strand = "+"
        if seqs is not None:
            s = seqs[tid]
            fwd = translate(s[start - 1:end]).rstrip("*")
            if "*" in fwd or not s[start - 1:start + 2] == "ATG":
                strand = "-"
        if strand == "-":
            start, end = tlen - end + 1, tlen - start + 1
        yield tid, max(start, 1), min(end, tlen), strand, x[7], None


def borf(path):
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if x[0] == "orf_id":
            continue
        tid, strand, L, a, b = x[1], x[3], int(x[4]), int(x[5]), int(x[6])
        if strand == "-":
            a, b = L - b + 1, L - a + 1
        yield tid, a, b, strand, ".", None


def rnasamba(path, fasta):
    seqs = read_fasta(fasta)
    for tid, prot in read_fasta(path).items():
        s, hit = seqs[tid], None
        for f in range(3):
            i = translate(s[f:]).find(prot)
            if i >= 0:
                a = f + 3 * i + 1
                b = a + 3 * len(prot) - 1
                if s[b:b + 3] and CODE.get(s[b:b + 3]) == "*":
                    b += 3
                hit = (a, b)
                break
        if hit is None:
            print("rnasamba: protein not found in %s" % tid, file=sys.stderr)
            continue
        yield tid, hit[0], hit[1], "+", ".", None


def angel(path):
    for line in open(path):
        if not line.startswith(">"):
            continue
        m = re.match(r">(\S+)\|\S+ type:(\S*) .*strand:([+-]) pos:(\d+)-(\d+)", line)
        yield m.group(1), int(m.group(4)), int(m.group(5)), m.group(3), m.group(2), None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("format", choices=["transdecoder_gff3", "prodigal_gff", "gmst_gff", "codan_gtf",
                                       "cpat_tsv", "cpc2_txt", "borf_txt", "rnasamba", "angel_cds"])
    ap.add_argument("file")
    ap.add_argument("--plus-only", action="store_true", help="prodigal: keep + strand predictions only")
    ap.add_argument("--fasta", help="input transcripts (cpat_tsv: restore ID case; cpc2_txt: strand; rnasamba: required)")
    ap.add_argument("--cutoff", type=float, default=0.0, help="cpat_tsv: minimum coding probability")
    ap.add_argument("--longest-per-transcript", action="store_true",
                    help="keep only the longest prediction per transcript (ties: first listed)")
    args = ap.parse_args()

    if args.format == "transdecoder_gff3":
        rows = transdecoder(args.file)
    elif args.format in ("prodigal_gff", "gmst_gff"):
        rows = prodigal(args.file, args.plus_only)
    elif args.format == "codan_gtf":
        rows = codan(args.file)
    elif args.format == "cpat_tsv":
        rows = cpat(args.file, args.cutoff, args.fasta)
    elif args.format == "cpc2_txt":
        rows = cpc2(args.file, args.fasta)
    elif args.format == "borf_txt":
        rows = borf(args.file)
    elif args.format == "angel_cds":
        rows = angel(args.file)
    else:
        rows = rnasamba(args.file, args.fasta)

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
