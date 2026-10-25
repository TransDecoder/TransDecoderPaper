#!/usr/bin/env python3
"""Real 5' UTR lengths and first-ATG distances from the reference annotations (evidence E29).

usage: python 03_utr_annotation.py <org> <out_prefix>
   (P = /broad/hptmp/bhaas/TRANSDECODER_PAPER; reads $P/data/annot via
    TransDecoderPaper/benchmark_v2/truth_sources.tsv and parses with build_truth_from_annotation.py)

All eligible protein-coding transcripts (same checks as the truth builder: nuclear, CDS in exons
and contiguous, ATG...stop, no internal stop, ACGT only; for mouse transcript_type protein_coding
and no *_NF tags), not one per gene and without added or clipped UTRs. Human: all MANE transcripts
(the MANE GTF). Writes
  <out_prefix>.tx.tsv.gz   per transcript: utr5 (nt, 0 = none annotated), cds_len, open5 (no in-frame
                           stop in the 5' UTR, so an ORF would run open to the transcript start),
                           first_atg (no in-frame ATG in the UTR downstream of the last in-frame
                           stop: the annotated start is the first in-frame ATG), truth_gene_rep
                           (1 if it is the transcript the benchmark keeps for its gene)
  <out_prefix>.trunc.tsv   distance d (nt) from a cut at a codon boundary inside the CDS to the first
                           in-frame ATG downstream of the cut's first codon (the UTR start
                           refinement would imply), counted over every codon position of every
                           CDS (uniform cut, weight 1) and weighted by the benchmark's N(0.5,0.1)
                           cut density ('normal'); cuts with no ATG before the stop are counted
                           as d = -1.
"""
import collections
import gzip
import math
import sys
from pathlib import Path

P = Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER")
sys.path.insert(0, str(P / "TransDecoderPaper/benchmark_v2/scripts"))
import build_truth_from_annotation as B  # noqa: E402
from Bio.Data import CodonTable  # noqa: E402

org, outp = sys.argv[1], sys.argv[2]
src = {}
for line in open(P / "TransDecoderPaper/benchmark_v2/truth_sources.tsv"):
    x = line.rstrip("\n").split("\t")
    src[x[0]] = x
_, gpath, apath, extra = src[org][:4] if len(src[org]) >= 4 else src[org] + [""]
table = 6 if "--table 6" in extra else 1
stops = set(CodonTable.unambiguous_dna_by_id[table].stop_codons)
fmt = "gtf" if apath.endswith((".gtf", ".gtf.gz")) else "gff3"
genome = B.read_genome(str(P / "data/annot" / gpath))
txs = B.load_annotation(str(P / "data/annot" / apath), fmt, {"mRNA", "transcript"})
excl_seq = __import__("re").compile(B.DEFAULT_EXCLUDE_SEQIDS)
excl_tags = {"cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"} if org == "mouse" else set()

rows, best = [], {}
for t in txs.values():
    if excl_seq.match(t.seqid or ""):
        continue
    if org == "mouse" and t.attrs.get("transcript_type", ["protein_coding"])[0] != "protein_coding":
        continue
    if set(t.attrs.get("tag", [])) & excl_tags:
        continue
    if t.cds_table is not None and t.cds_table != table:
        continue
    rec = B.tx_to_record(t, genome)
    if isinstance(rec, str):
        continue
    seq, s, e = rec
    if set(seq) - set("ACGT"):
        continue
    if seq[e - 3:e] not in stops and e + 3 <= len(seq) and seq[e:e + 3] in stops and (e - s) % 3 == 0:
        e += 3
    cds = seq[s:e]
    if len(cds) % 3 or len(cds) < 6 or not cds.startswith("ATG") or cds[-3:] not in stops:
        continue
    if any(cds[i:i + 3] in stops for i in range(0, len(cds) - 3, 3)):
        continue
    # scan the 5' UTR upstream in frame
    open5, first_atg = 1, 1
    i = s - 3
    while i >= 0:
        c = seq[i:i + 3]
        if c in stops:
            open5 = 0
            break
        if c == "ATG":
            first_atg = 0
        i -= 3
    rows.append([t.id, t.gene or t.id, s, len(cds), open5, first_atg, cds])
    g = t.gene or t.id
    key = (e - s, len(seq), t.id)
    if g not in best or key > best[g]:
        best[g] = key
rep = {k[2] for k in best.values()}

trunc_u, trunc_n = collections.Counter(), collections.Counter()
with gzip.open(outp + ".tx.tsv.gz", "wt") as o:
    o.write("org\ttranscript_id\tgene_id\tutr5\tcds_len\topen5\tfirst_atg\ttruth_gene_rep\n")
    for tid, gid, u, cl, op, fa, cds in rows:
        o.write(f"{org}\t{tid}\t{gid}\t{u}\t{cl}\t{op}\t{fa}\t{int(tid in rep)}\n")
    # truncation distances over the gene representatives (one CDS per gene)
    for tid, gid, u, cl, op, fa, cds in rows:
        if tid not in rep:
            continue
        n = cl // 3
        mets = [k for k in range(n) if cds[3 * k:3 * k + 3] == "ATG"]
        j = 0
        for c in range(1, n - 1):  # fragment's first codon is c; ATG searched from codon c+1
            while j < len(mets) and mets[j] <= c:
                j += 1
            d = 3 * (mets[j] - c) if j < len(mets) and mets[j] < n - 1 else -1
            f = c / n
            trunc_u[d] += 1
            if 0.05 <= f <= 0.95:
                trunc_n[d] += math.exp(-0.5 * ((f - 0.5) / 0.1) ** 2)
with open(outp + ".trunc.tsv", "w") as o:
    o.write("org\td\tuniform\tnormal\n")
    for d in sorted(set(trunc_u) | set(trunc_n)):
        o.write(f"{org}\t{d}\t{trunc_u[d]}\t{trunc_n[d]:.3f}\n")
print(org, len(rows), "transcripts", len(rep), "genes", file=sys.stderr)
