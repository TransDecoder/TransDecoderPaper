#!/usr/bin/env python3
"""Match TransDecoder genome CDS models to GENCODE Ribo-seq ORFs.

usage: match_riboseq.py --pred <genome.gff3> --phase1 <Ribo-seq_ORFs.all.bed>
                        --phase2 <Ribo-seq_ORFs.comprehensive.bed> [--ids <file of model ids>] > out.tsv
       (PYTHONPATH must include TransDecoderPaper/genome_annotation/scripts)

A prediction matches a Ribo-seq ORF at the 'stop' level when both are on the same chromosome and
strand and end at the same stop codon (genomic 3' end of the CDS including the stop), and 'exact'
when the CDS segments are identical. Phase 1 itemRgb gives the peptide-support tier.
Output per prediction: model_id, orf_type, aa, stop_match, exact_match, catalog, ribo_orf, ribo_aa,
peptide_support.
"""
import argparse
from eval_genome_cds import read_pred_gff3, merge

TIERS = {"0,128,0": "peptidein_nonHLA", "200,100,0": "peptidein_HLA", "30,100,200": "nonHLA_peptide",
         "210,150,30": "HLA_tier1B", "210,190,80": "HLA_tier2B", "150,150,150": "none"}


def read_bed(path, catalog):
    orfs = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        start, ts, te = int(f[1]), int(f[6]), int(f[7])
        sizes = [int(x) for x in f[10].strip(",").split(",")]
        starts = [int(x) for x in f[11].strip(",").split(",")]
        segs = []
        for s, z in zip(starts, sizes):
            a, b = start + s, start + s + z   # 0-based half-open
            a, b = max(a, ts), min(b, te)
            if a < b:
                segs.append((a + 1, b))
        segs = merge(segs)
        orfs.append(dict(name=f[3], chrom=f[0], strand=f[5], segs=segs,
                         aa=sum(b - a + 1 for a, b in segs) // 3 - 1,
                         tier=TIERS.get(f[8], "NA") if catalog == "phase1" else "NA", catalog=catalog))
    return orfs


def stop_key(chrom, strand, segs):
    return (chrom, strand, segs[-1][1] if strand == "+" else segs[0][0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--phase1", required=True)
    ap.add_argument("--phase2", required=True)
    ap.add_argument("--ids")
    args = ap.parse_args()
    by_stop, by_exact = {}, {}
    for o in read_bed(args.phase1, "phase1") + read_bed(args.phase2, "phase2"):
        k = stop_key(o["chrom"], o["strand"], o["segs"])
        # prefer phase 1 (peptide tiers) when an ORF is in both
        if k not in by_stop or (o["catalog"] == "phase1" and by_stop[k]["catalog"] != "phase1"):
            by_stop[k] = o
        by_exact[(o["chrom"], o["strand"], o["segs"])] = o
    preds = read_pred_gff3(args.pred)
    keep = set(open(args.ids).read().split()) if args.ids else None
    print("model_id\torf_type\taa\tstop_match\texact_match\tcatalog\tribo_orf\tribo_aa\tpeptide_support")
    for mid, p in preds.items():
        if keep is not None and mid not in keep:
            continue
        aa = sum(b - a + 1 for a, b in p["segs"]) // 3 - 1
        o = by_stop.get(stop_key(p["chrom"], p["strand"], p["segs"]))
        ex = (p["chrom"], p["strand"], p["segs"]) in by_exact
        if o:
            print(f"{mid}\t{p['orf_type']}\t{aa}\t1\t{int(ex)}\t{o['catalog']}\t{o['name']}\t{o['aa']}\t{o['tier']}")
        else:
            print(f"{mid}\t{p['orf_type']}\t{aa}\t0\t0\tNA\tNA\tNA\tNA")


if __name__ == "__main__":
    main()
