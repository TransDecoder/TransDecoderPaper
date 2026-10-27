#!/usr/bin/env python3
"""Parse GENCODE v50 (primary assembly GTF) into a per-transcript table and select the
transcripts to simulate.

usage: 01_gencode_tables.py <gencode.gtf[.gz]> <all_transcripts.fa> <outdir>

<all_transcripts.fa> is `gffread -w` of the same GTF (all transcripts; headers = transcript_id).

Writes to <outdir>:
  gencode_v50.tx.tsv.gz   every transcript: transcript_id gene_id gene_name gene_type
                          transcript_type tags chrom strand exons cds
                          (exons / cds = genomic segments 'start-end,...' ascending, 1-based;
                          cds = GTF CDS + stop_codon merged, i.e. including the stop codon)
  kept.tsv                the simulated set: transcript_id class(coding|lncRNA) gene_id gene_name
                          tx_len cds_start cds_end (1-based transcript coords incl. stop; NA for lncRNA)
  kept.fa                 their sequences (from <all_transcripts.fa>)
  selection_counts.tsv    counts at each filter

Coding = transcript_type protein_coding with an annotated CDS, none of the tags cds_start_NF,
cds_end_NF, mRNA_start_NF, mRNA_end_NF, and a complete CDS (ATG ... stop, length % 3 == 0, no
internal stop, standard code). Negatives = every transcript_type lncRNA transcript.
"""
import collections
import gzip
import os
import re
import sys

STOPS = {"TAA", "TAG", "TGA"}
NF_TAGS = {"cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"}
ATTR = re.compile(r'(\S+) "([^"]*)"')


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def merge(segs):
    segs = sorted(segs)
    out = []
    for s, e in segs:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def read_fasta(path):
    seqs = {}
    name, buf = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def genomic_to_tx(exons, strand, pos):
    """1-based transcript coordinate of genomic position pos (must be exonic)."""
    off = 0
    ex = exons if strand == "+" else exons[::-1]
    for s, e in ex:
        if s <= pos <= e:
            return off + (pos - s + 1 if strand == "+" else e - pos + 1)
        off += e - s + 1
    return None


def main():
    gtf, fa, outdir = sys.argv[1:4]
    os.makedirs(outdir, exist_ok=True)
    info, exons, cds = {}, collections.defaultdict(list), collections.defaultdict(list)
    with opener(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            feat = f[2]
            if feat not in ("transcript", "exon", "CDS", "stop_codon"):
                continue
            a = f[8]
            tid = a.split('transcript_id "', 1)[1].split('"', 1)[0]
            s, e = int(f[3]), int(f[4])
            if feat == "transcript":
                kv = ATTR.findall(a)
                d = {}
                tags = []
                for k, v in kv:
                    if k == "tag":
                        tags.append(v)
                    else:
                        d.setdefault(k, v)
                info[tid] = (d["gene_id"], d.get("gene_name", "."), d.get("gene_type", "."),
                             d.get("transcript_type", "."), ",".join(tags) or ".", f[0], f[6])
            elif feat == "exon":
                exons[tid].append((s, e))
            else:
                cds[tid].append((s, e))

    seqs = read_fasta(fa)
    cnt = collections.OrderedDict()
    def inc(k, n=1):
        cnt[k] = cnt.get(k, 0) + n

    with gzip.open(os.path.join(outdir, "gencode_v50.tx.tsv.gz"), "wt") as out, \
            open(os.path.join(outdir, "kept.tsv"), "w") as kept, \
            open(os.path.join(outdir, "kept.fa"), "w") as kfa:
        out.write("transcript_id\tgene_id\tgene_name\tgene_type\ttranscript_type\ttags\tchrom\tstrand\texons\tcds\n")
        kept.write("transcript_id\tclass\tgene_id\tgene_name\ttx_len\tcds_start\tcds_end\n")
        for tid, (gid, gname, gtype, ttype, tags, chrom, strand) in info.items():
            ex = merge(exons[tid])
            cd = merge(cds[tid]) if cds[tid] else []
            out.write("\t".join([tid, gid, gname, gtype, ttype, tags, chrom, strand,
                                 ",".join(f"{s}-{e}" for s, e in ex),
                                 ",".join(f"{s}-{e}" for s, e in cd) or "."]) + "\n")
            inc("all_transcripts")
            if ttype == "lncRNA":
                inc("lncRNA")
                seq = seqs[tid]
                kept.write(f"{tid}\tlncRNA\t{gid}\t{gname}\t{len(seq)}\tNA\tNA\n")
                kfa.write(f">{tid}\n{seq}\n")
                continue
            if ttype != "protein_coding":
                continue
            inc("protein_coding")
            if not cd:
                inc("pc_drop_no_cds")
                continue
            tagset = set(tags.split(","))
            if tagset & NF_TAGS:
                inc("pc_drop_NF_tag")
                for t in sorted(tagset & NF_TAGS):
                    inc(f"pc_NF_tag:{t}")
                continue
            seq = seqs[tid]
            if strand == "+":
                a, b = genomic_to_tx(ex, strand, cd[0][0]), genomic_to_tx(ex, strand, cd[-1][1])
            else:
                a, b = genomic_to_tx(ex, strand, cd[-1][1]), genomic_to_tx(ex, strand, cd[0][0])
            c = seq[a - 1:b]
            cdslen = sum(e - s + 1 for s, e in cd)
            if len(c) != cdslen:
                inc("pc_drop_cds_not_contiguous_in_tx")
                continue
            reason = None
            if len(c) % 3:
                reason = "pc_drop_len_not_mult3"
            elif not c.startswith("ATG"):
                reason = "pc_drop_no_ATG"
            elif c[-3:] not in STOPS:
                reason = "pc_drop_no_stop"
            elif any(c[i:i + 3] in STOPS for i in range(0, len(c) - 3, 3)):
                reason = "pc_drop_internal_stop"
            if reason:
                inc(reason)
                if any("seleno" in t for t in tagset):
                    inc(reason + ":selenoprotein")
                continue
            inc("coding_kept")
            kept.write(f"{tid}\tcoding\t{gid}\t{gname}\t{len(seq)}\t{a}\t{b}\n")
            kfa.write(f">{tid}\n{seq}\n")
    with open(os.path.join(outdir, "selection_counts.tsv"), "w") as fh:
        for k, v in cnt.items():
            fh.write(f"{k}\t{v}\n")
    for k, v in cnt.items():
        print(k, v, sep="\t")


if __name__ == "__main__":
    main()
