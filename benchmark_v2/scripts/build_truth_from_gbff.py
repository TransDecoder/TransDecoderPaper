#!/usr/bin/env python3
"""Build a coding-region truth set (one transcript per gene) from a RefSeq rna.gbff file.

This follows the 2018 benchmark (TransDecoder_benchmarking_data/utilities/
assign_longest_refcds_to_refmrna_per_gene.pl) but reads CDS coordinates directly
from the RefSeq transcript records instead of re-locating CDS strings in mRNAs.

A transcript is eligible when it:
  - is an mRNA with an accepted accession prefix (default NM_; e.g. NM_,XM_ for
    genomes with model-only annotation),
  - optionally carries the 'MANE Select' keyword (--mane),
  - has exactly one CDS, a single contiguous span on the + strand of the record,
    not flagged partial, with no /exception, /ribosomal_slippage or /transl_except,
  - starts with ATG, ends with a stop codon of its translation table, has length
    divisible by 3 and no internal stops, and its translation matches /translation,
  - contains only ACGT.

Per gene (GeneID db_xref, else /gene or /locus_tag) the transcript with the longest
CDS is kept (ties: longer mRNA, then accession). UTRs longer than --max-utr are
trimmed (the 2018 benchmark used 250; use --max-utr -1 to keep full UTRs).

Outputs (prefix given by --out):
  <out>.fasta      transcripts, header '>accession'
  <out>.truth.tsv  transcript_id, cds_start, cds_end (1-based, incl. stop), strand,
                   transcript_length, gene_id, transl_table, type
  <out>.log        counts of exclusions by reason
"""
import argparse
import collections
import gzip
import sys

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition


def gene_key(feat_quals, rec_id):
    for x in feat_quals.get("db_xref", []):
        if x.startswith("GeneID:"):
            return x
    for k in ("gene", "locus_tag"):
        if k in feat_quals:
            return f"{k}:{feat_quals[k][0]}"
    return f"acc:{rec_id}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gbff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefixes", default="NM_")
    ap.add_argument("--mane", action="store_true", help="require the 'MANE Select' keyword")
    ap.add_argument("--max-utr", type=int, default=250)
    ap.add_argument("--min-cds", type=int, default=0, help="minimum CDS length in nt incl. stop")
    args = ap.parse_args()

    prefixes = tuple(args.prefixes.split(","))
    reasons = collections.Counter()
    best = {}  # gene -> candidate dict

    opener = gzip.open if args.gbff.endswith(".gz") else open
    with opener(args.gbff, "rt") as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            reasons["records"] += 1
            if not rec.id.startswith(prefixes):
                reasons["prefix"] += 1
                continue
            if rec.annotations.get("molecule_type", "mRNA") != "mRNA":
                reasons["not_mRNA"] += 1
                continue
            if args.mane and "MANE Select" not in rec.annotations.get("keywords", []):
                reasons["not_MANE"] += 1
                continue
            cdss = [f for f in rec.features if f.type == "CDS"]
            if len(cdss) != 1:
                reasons["n_cds_ne_1"] += 1
                continue
            cds = cdss[0]
            q = cds.qualifiers
            loc = cds.location
            if len(loc.parts) != 1 or loc.strand != 1:
                reasons["cds_not_single_plus_span"] += 1
                continue
            if not (isinstance(loc.start, ExactPosition) and isinstance(loc.end, ExactPosition)):
                reasons["cds_partial"] += 1
                continue
            if any(k in q for k in ("exception", "ribosomal_slippage", "transl_except", "artificial_location")):
                reasons["cds_exception"] += 1
                continue
            if q.get("codon_start", ["1"])[0] != "1":
                reasons["codon_start_ne_1"] += 1
                continue
            table = int(q.get("transl_table", ["1"])[0])
            seq = str(rec.seq).upper()
            if set(seq) - set("ACGT"):
                reasons["non_ACGT"] += 1
                continue
            s, e = int(loc.start), int(loc.end)  # 0-based half open
            cds_seq = seq[s:e]
            if len(cds_seq) % 3 or len(cds_seq) < max(6, args.min_cds):
                reasons["cds_len"] += 1
                continue
            stops = set(CodonTable.unambiguous_dna_by_id[table].stop_codons)
            if not cds_seq.startswith("ATG"):
                reasons["non_ATG_start"] += 1
                continue
            if cds_seq[-3:] not in stops:
                reasons["no_terminal_stop"] += 1
                continue
            if any(cds_seq[i:i + 3] in stops for i in range(0, len(cds_seq) - 3, 3)):
                reasons["internal_stop"] += 1
                continue
            prot = str(Seq(cds_seq).translate(table=table, cds=False)).rstrip("*")
            given = q.get("translation", [None])[0]
            # first residue: alternative initiation codons translate as M in /translation
            if given is not None and prot[1:] != given[1:]:
                reasons["translation_mismatch"] += 1
                continue

            g = gene_key(q, rec.id)
            cand = dict(acc=rec.id, seq=seq, s=s, e=e, table=table, gene=g)
            cur = best.get(g)
            key = (e - s, len(seq), rec.id)
            if cur is None or key > (cur["e"] - cur["s"], len(cur["seq"]), cur["acc"]):
                if cur is not None:
                    reasons["not_longest_for_gene"] += 1
                best[g] = cand
            else:
                reasons["not_longest_for_gene"] += 1

    with open(args.out + ".fasta", "w") as fa, open(args.out + ".truth.tsv", "w") as tt:
        tt.write("transcript_id\tcds_start\tcds_end\tstrand\ttranscript_length\tgene_id\ttransl_table\ttype\n")
        for g in sorted(best, key=lambda k: best[k]["acc"]):
            c = best[g]
            seq, s, e = c["seq"], c["s"], c["e"]
            if args.max_utr >= 0:
                lo = max(0, s - args.max_utr)
                hi = min(len(seq), e + args.max_utr)
                seq, s, e = seq[lo:hi], s - lo, e - lo
            fa.write(f">{c['acc']}\n{seq}\n")
            tt.write(f"{c['acc']}\t{s + 1}\t{e}\t+\t{len(seq)}\t{g}\t{c['table']}\tcomplete\n")
    reasons["kept_genes"] = len(best)
    with open(args.out + ".log", "w") as lg:
        for k, v in reasons.most_common():
            lg.write(f"{k}\t{v}\n")
    sys.stderr.write("".join(f"{k}\t{v}\n" for k, v in reasons.most_common()))


if __name__ == "__main__":
    sys.exit(main())
