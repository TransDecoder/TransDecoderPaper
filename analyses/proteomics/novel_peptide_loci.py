#!/usr/bin/env python3
"""Where unannotated ORFs that share MS peptides with known proteins lie, relative to the parent gene.

usage: novel_peptide_loci.py <ms_peptide_support.novel_with_peptides.tsv> <gencode.gtf>
                             <name>=<transdecoder.genome.gff3> [...] -o out_prefix

For each novel-class ORF with peptides (ms_peptide_support.py output), the parent gene is the gene of its
best-matching reference protein (best_ref_genes). The ORF's genome span comes from the genome-projected
GFF3. Location class (relative to the nearest locus of any parent gene of that name):
  parent_overlap   ORF span overlaps the parent gene span (intron, UTR or another frame)
  near_le2Mb       same chromosome, within 2 Mb (tandem / segmental duplication)
  same_chr_far     same chromosome, farther away
  other_chr        different chromosome (e.g. processed pseudogene)
  parent_not_found no GENCODE gene with that name
Independently, the ORF is flagged if its span overlaps a GENCODE gene whose gene_type ends in
'pseudogene' (the overlapping pseudogene names are listed).

Outputs: <out_prefix>.tsv (one row per ORF), <out_prefix>.summary.tsv
"""
import argparse
import bisect
import re
from collections import Counter, defaultdict

ATTR = re.compile(r'(\S+) "([^"]*)"')


def read_gencode(path):
    by_name = defaultdict(list)          # gene_name -> [(chrom, start, end)]
    pseudo = defaultdict(list)           # chrom -> [(start, end, name)]
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        a = dict(ATTR.findall(f[8]))
        s, e = int(f[3]), int(f[4])
        by_name[a.get("gene_name", "")].append((f[0], s, e))
        if a.get("gene_type", "").endswith("pseudogene"):
            pseudo[f[0]].append((s, e, a.get("gene_name", a["gene_id"])))
    idx = {}
    for c, iv in pseudo.items():
        iv.sort()
        idx[c] = (iv, [x[0] for x in iv], max(e - s for s, e, _ in iv))
    return by_name, idx


def pseudo_overlaps(idx, chrom, s, e):
    if chrom not in idx:
        return []
    iv, starts, maxlen = idx[chrom]
    i = bisect.bisect_left(starts, s - maxlen)
    out = []
    while i < len(iv) and iv[i][0] <= e:
        if iv[i][1] >= s:
            out.append(iv[i][2])
        i += 1
    return out


def read_spans(path):
    span = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        mid = re.search(r"Parent=([^;]+)", f[8]).group(1)
        s, e = int(f[3]), int(f[4])
        c, a, b = span.get(mid, (f[0], s, e))
        span[mid] = (c, min(a, s), max(b, e))
    return span


def classify(orf, loci):
    if not loci:
        return "parent_not_found", ""
    c, s, e = orf
    best = None
    for lc, ls, le in loci:
        if lc != c:
            cand = (3, 0)
        elif ls <= e and le >= s:
            cand = (0, 0)
        else:
            d = max(ls - e, s - le)
            cand = (1 if d <= 2_000_000 else 2, d)
        best = cand if best is None or cand < best else best
    return ["parent_overlap", "near_le2Mb", "same_chr_far", "other_chr"][best[0]], str(best[1] or "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("novel_tsv")
    ap.add_argument("gencode_gtf")
    ap.add_argument("pred", nargs="+", help="name=genome.gff3")
    ap.add_argument("-o", "--out-prefix", required=True)
    args = ap.parse_args()
    by_name, pidx = read_gencode(args.gencode_gtf)
    spans = {n: read_spans(p) for n, p in (x.split("=", 1) for x in args.pred)}
    summ = Counter()
    with open(f"{args.out_prefix}.tsv", "w") as out:
        out.write("assembly\tmodel_id\taa\tn_peptides\tparent_genes\torf_locus\tlocation_class\t"
                  "distance_bp\toverlapping_pseudogenes\n")
        for line in open(args.novel_tsv):
            f = line.rstrip("\n").split("\t")
            if f[0] == "assembly" or f[0] not in spans:
                continue
            asm, mid = f[0], f[1]
            orf = spans[asm].get(mid)
            if orf is None:
                continue
            genes = [g for g in re.split(r"[;,]", f[8]) if g]
            loci = [l for g in genes for l in by_name.get(g, [])]
            cls, dist = classify(orf, loci)
            ps = sorted(set(pseudo_overlaps(pidx, *orf)))
            summ[(asm, cls, bool(ps))] += 1
            summ[(asm, "all", bool(ps))] += 1
            out.write(f"{asm}\t{mid}\t{f[2]}\t{f[3]}\t{','.join(genes)}\t{orf[0]}:{orf[1]}-{orf[2]}\t{cls}\t"
                      f"{dist}\t{','.join(ps)}\n")
    with open(f"{args.out_prefix}.summary.tsv", "w") as fh:
        fh.write("assembly\tlocation_class\tn\toverlaps_pseudogene\n")
        for asm in spans:
            for cls in ("all", "parent_overlap", "near_le2Mb", "same_chr_far", "other_chr", "parent_not_found"):
                n = summ[(asm, cls, True)] + summ[(asm, cls, False)]
                row = f"{asm}\t{cls}\t{n}\t{summ[(asm, cls, True)]}"
                fh.write(row + "\n")
                print(row)


if __name__ == "__main__":
    main()
