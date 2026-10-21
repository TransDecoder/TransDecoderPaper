#!/usr/bin/env python3
"""Ribo-seq support for genome-mapped ORF predictions (sample-matched footprints).

usage: riboseq_orf_support.py --bam footprints.genome.bam --mane MANE.gtf[.gz]
                              --pred <name>=<genome.gff3>[,<pred_status.tsv>] [...] -o out_prefix
                              [--min-psites 10] [--alpha 0.01]

1. Calibration on MANE Select CDSs (spliced, stop codon excluded):
   - footprint strand (sense or antisense to the CDS);
   - for each read length, the 5'-end-to-P-site offset (0..20) maximising the frame-0 fraction;
   - a read length is used only if that fraction is >= --min-length-frame (default 0.5) and it has
     >= 1000 reads.
2. For each predicted CDS: P-sites of usable reads are mapped into spliced CDS coordinates, and counted
   by frame (frame 0 = codon position 1 of the predicted reading frame). The first and last codons
   are excluded (initiation/termination peaks).
   supported = P-sites >= --min-psites, frame 0 the largest, and binomial P(frame0 | n, 1/3) < --alpha.

Outputs:
- <out_prefix>.calibration.tsv
- <out_prefix>.<name>.orf_support.tsv (one row per predicted model; with pred_status, the
  GENCODE-overlap class of the distinct prediction is added)
- <out_prefix>.summary.tsv
"""
import argparse
import gzip
import math
import re
import sys
from collections import Counter, defaultdict

import pysam

ATTR = re.compile(r'(\S+) "([^"]*)"')


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_mane_cds(path):
    cds = defaultdict(list)
    strand, chrom = {}, {}
    for line in opener(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "CDS" or "MANE_Select" not in f[8]:
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        cds[tid].append((int(f[3]), int(f[4])))
        strand[tid], chrom[tid] = f[6], f[0]
    return {t: (chrom[t], strand[t], sorted(s)) for t, s in cds.items()}


def read_pred_gff3(path):
    models = {}
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] == "mRNA":
            mid = re.search(r"ID=([^;]+)", f[8]).group(1)
            t = re.search(r"ORF type:(\w+)", f[8])
            models[mid] = [f[0], f[6], [], t.group(1) if t else "NA"]
        elif f[2] == "CDS":
            mid = re.search(r"Parent=([^;]+)", f[8]).group(1)
            models[mid][2].append((int(f[3]), int(f[4])))
    for m in models.values():
        m[2].sort()
    return models


def cds_index(segs, strand):
    """genomic position -> 0-based offset in the spliced CDS (5'->3')."""
    pos = []
    for a, b in segs:
        pos.extend(range(a, b + 1))
    if strand == "-":
        pos.reverse()
    return {p: i for i, p in enumerate(pos)}, len(pos)


def read_5p_positions(read):
    """1-based reference positions of the aligned bases, ordered from the read's 5' end."""
    pos = [p + 1 for p in read.get_reference_positions()]
    return pos[::-1] if read.is_reverse else pos


def iter_reads(bam, chrom, segs):
    seen = set()
    for a, b in segs:
        for r in bam.fetch(chrom, a - 1, b):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 10:
                continue
            key = (r.query_name, r.reference_start, r.is_reverse)
            if key in seen:
                continue
            seen.add(key)
            yield r


def calibrate(bam, mane, max_tx=4000):
    sense = anti = 0
    lens = Counter()
    frame_by_len_off = defaultdict(lambda: [[0, 0, 0] for _ in range(21)])
    for i, (tid, (chrom, strand, segs)) in enumerate(mane.items()):
        if i >= max_tx:
            break
        if chrom not in bam.references:
            continue
        idx, n = cds_index(segs, strand)
        for r in iter_reads(bam, chrom, segs):
            rstrand = "-" if r.is_reverse else "+"
            if rstrand == strand:
                sense += 1
            else:
                anti += 1
                continue
            L = r.query_alignment_length
            lens[L] += 1
            p5 = read_5p_positions(r)
            for off in range(min(21, len(p5))):
                j = idx.get(p5[off])
                if j is not None and 3 <= j < n - 3:
                    frame_by_len_off[L][off][j % 3] += 1
    calib = {}
    rows = []
    for L, offs in sorted(frame_by_len_off.items()):
        best = max(range(21), key=lambda o: offs[o][0] / max(sum(offs[o]), 1) if sum(offs[o]) >= 200 else 0)
        tot = sum(offs[best])
        frac = offs[best][0] / max(tot, 1)
        rows.append((L, lens[L], best, tot, *offs[best], round(frac, 3)))
        calib[L] = (best, frac, tot)
    return sense, anti, rows, calib


def binom_sf(k, n, p=1 / 3):
    """P(X >= k), X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    logp, logq = math.log(p), math.log(1 - p)
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * logp + (n - i) * logq
             for i in range(k, n + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--mane", required=True)
    ap.add_argument("--pred", action="append", required=True)
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--min-psites", type=int, default=10)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--min-length-frame", type=float, default=0.5)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    mane = read_mane_cds(args.mane)
    sense, anti, rows, calib = calibrate(bam, mane)
    use = {L: off for L, (off, frac, tot) in calib.items() if frac >= args.min_length_frame and tot >= 1000}
    with open(f"{args.out_prefix}.calibration.tsv", "w") as fh:
        fh.write(f"# MANE CDS reads: sense {sense}, antisense {anti}\n")
        fh.write("read_len\treads\tbest_offset\tpsites_in_cds\tframe0\tframe1\tframe2\tframe0_frac\tused\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + f"\t{int(r[0] in use)}\n")
    print(f"strand: sense {sense} antisense {anti}; usable lengths/offsets {use}", file=sys.stderr)
    if not use:
        sys.exit("no read length shows frame periodicity; stopping")
    fp_sense = sense >= anti

    summary = []
    for spec in args.pred:
        name, files = spec.split("=", 1)
        gff3, status = (files.split(",") + [None])[:2]
        klass = {}
        if status:
            for line in open(status):
                f = line.rstrip("\n").split("\t")
                if f[0] != "model_id":
                    klass[f[0]] = ("exact" if f[3] == "1" else "compatible" if f[4] == "1"
                                   else "novel" if f[5] == "1" else "other_overlap", int(f[2]))
        models = read_pred_gff3(gff3)
        out = open(f"{args.out_prefix}.{name}.orf_support.tsv", "w")
        out.write("model_id\torf_type\tclass\taa\tpsites\tframe0\tframe1\tframe2\tframe0_frac\tp_frame0\tsupported\n")
        agg = defaultdict(lambda: [0, 0, 0])  # class -> [n, with>=min psites, supported]
        for mid, (chrom, strand, segs, otype) in models.items():
            if status and mid not in klass:
                continue  # not the distinct representative
            if chrom not in bam.references or not segs:
                continue
            idx, n = cds_index(segs, strand)
            fr = [0, 0, 0]
            for r in iter_reads(bam, chrom, segs):
                L = r.query_alignment_length
                if L not in use or (("-" if r.is_reverse else "+") == strand) != fp_sense:
                    continue
                p5 = read_5p_positions(r)
                if use[L] >= len(p5):
                    continue
                j = idx.get(p5[use[L]])
                if j is not None and 3 <= j < n - 3:
                    fr[j % 3] += 1
            tot = sum(fr)
            p = binom_sf(fr[0], tot) if tot else 1.0
            sup = int(tot >= args.min_psites and fr[0] == max(fr) and p < args.alpha)
            c, aa = klass.get(mid, ("NA", n // 3))
            out.write(f"{mid}\t{otype}\t{c}\t{aa}\t{tot}\t{fr[0]}\t{fr[1]}\t{fr[2]}\t"
                      f"{fr[0] / tot if tot else 0:.3f}\t{p:.3g}\t{sup}\n")
            for key in (c, f"{c}:{'short' if aa < 100 else 'long'}", "all"):
                agg[key][0] += 1
                agg[key][1] += tot >= args.min_psites
                agg[key][2] += sup
        out.close()
        for key, (n, cov, sup) in sorted(agg.items()):
            summary.append((name, key, n, cov, sup, f"{sup / max(n, 1):.4f}", f"{sup / max(cov, 1):.4f}"))
    with open(f"{args.out_prefix}.summary.tsv", "w") as fh:
        fh.write("assembly\tclass\tn\twith_min_psites\tsupported\tsupported_frac\tsupported_frac_of_covered\n")
        for r in summary:
            fh.write("\t".join(map(str, r)) + "\n")
            print("\t".join(map(str, r)))


if __name__ == "__main__":
    main()
