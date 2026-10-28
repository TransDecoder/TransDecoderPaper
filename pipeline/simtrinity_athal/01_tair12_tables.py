#!/usr/bin/env python3
"""TAIR12 (Arabidopsis) per-transcript tables and the transcript set to simulate (evidence E32).

The Arabidopsis counterpart of pipeline/simtrinity/01_gencode_tables.py (E30). It writes the same table
format, so label_assembly.py and the rest of the E30 pipeline can be reused unchanged.

usage: 01_tair12_tables.py <TAIR12.gff3.gz> <genome.fa.gz> <outdir>

Transcripts are the GFF3 features that own exon features (mRNA, ncRNA, transcript [pseudogene], tRNA,
misc_RNA, miRNA_primary_transcript, ...). Features with no exon children but an RNA-like type (rRNA,
precursor_RNA) become single-exon transcripts, so that the salmon index covers them as well.
  transcript_type = 'protein_coding' for mRNA with CDS; ncRNA_class for ncRNA (lncRNA, antisense_RNA, ...);
  otherwise the GFF3 type ('pseudogene' for type 'transcript' with product=pseudogene).
CDS = GFF3 CDS segments. If the CDS does not end in a stop codon but the next 3 genomic bases are a stop,
those bases are added (records whether the annotation includes the stop codon).

Selection (same eligibility rules as E30):
  coding    nuclear (Chr1-5) protein_coding with ATG start, stop codon end, length % 3 == 0 and no internal
            stop (standard code). Organellar (ChrM, ChrC) transcripts are not simulated (RNA editing).
  negatives nuclear ncRNA of class lncRNA or antisense_RNA whose exons overlap no CDS (any TAIR12
            transcript) on the same strand; written to kept.tsv with class 'lncRNA'.
Writes to <outdir>: gencode_v50.tx.tsv.gz-format table 'tair12.tx.tsv.gz', all_transcripts.fa (salmon index),
tair12.exons.gtf (exon GTF for the STAR index), kept.tsv, kept.fa, selection_counts.tsv.
"""
import collections
import gzip
import os
import re
import sys

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
NEG_CLASSES = {"lncRNA", "antisense_RNA"}
NUCLEAR = {"Chr1", "Chr2", "Chr3", "Chr4", "Chr5"}


def attrs(s):
    return dict(kv.split("=", 1) for kv in s.strip().strip(";").split(";") if "=" in kv)


def read_genome(path):
    g, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    g[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        g[name] = "".join(buf).upper()
    return g


def merge(segs):
    out = []
    for s, e in sorted(segs):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def spliced(g, chrom, strand, segs):
    s = "".join(g[chrom][a - 1:b] for a, b in segs)
    return s if strand == "+" else s.translate(COMP)[::-1]


def main():
    gff, fa, outdir = sys.argv[1:4]
    os.makedirs(outdir, exist_ok=True)
    feats, exons, cds, gene_name = {}, collections.defaultdict(list), collections.defaultdict(list), {}
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            s, e = int(f[3]), int(f[4])
            if f[2] == "exon":
                for p in a.get("Parent", "").split(","):
                    exons[p].append((s, e))
            elif f[2] == "CDS":
                for p in a.get("Parent", "").split(","):
                    cds[p].append((s, e))
            elif f[2] in ("gene", "pseudogene"):
                gene_name[a.get("ID")] = a.get("Name", a.get("ID"))
            elif "ID" in a:
                feats[a["ID"]] = (f[2], f[0], f[6], s, e, a)
    g = read_genome(fa)
    cnt = collections.OrderedDict()

    def inc(k, n=1):
        cnt[k] = cnt.get(k, 0) + n

    rows = []
    for tid, (ftype, chrom, strand, s, e, a) in feats.items():
        ex = merge(exons[tid]) if tid in exons else None
        if ex is None:
            if ftype not in ("rRNA", "precursor_RNA"):
                continue
            ex = [(s, e)]
        if ftype == "mRNA":
            ttype = "protein_coding"
        elif ftype == "ncRNA":
            ttype = a.get("ncRNA_class", "ncRNA")
        elif ftype == "transcript":
            ttype = "pseudogene" if "pseudogene" in a.get("product", "") else "transcript"
        else:
            ttype = ftype
        gid = a.get("Parent", tid).split(",")[0]
        rows.append([tid, gid, gene_name.get(gid, gid), ttype, chrom, strand, ex, merge(cds[tid]) if tid in cds else []])

    # same-strand CDS index for the lncRNA overlap filter
    cds_iv = collections.defaultdict(list)
    for r in rows:
        for s, e in r[7]:
            cds_iv[(r[4], r[5])].append((s, e))
    for k in cds_iv:
        cds_iv[k].sort()
    starts = {k: [s for s, e in v] for k, v in cds_iv.items()}
    import bisect

    def overlaps_cds(chrom, strand, ex):
        v, st = cds_iv.get((chrom, strand), []), starts.get((chrom, strand), [])
        for s, e in ex:
            i = bisect.bisect_right(st, e)
            for cs, ce in v[max(0, i - 200):i]:
                if cs <= e and ce >= s:
                    return True
        return False

    stop_added = 0
    with gzip.open(os.path.join(outdir, "tair12.tx.tsv.gz"), "wt") as out, \
            open(os.path.join(outdir, "all_transcripts.fa"), "w") as allfa, \
            open(os.path.join(outdir, "tair12.exons.gtf"), "w") as gtf, \
            open(os.path.join(outdir, "kept.tsv"), "w") as kept, open(os.path.join(outdir, "kept.fa"), "w") as kfa:
        out.write("transcript_id\tgene_id\tgene_name\tgene_type\ttranscript_type\ttags\tchrom\tstrand\texons\tcds\n")
        kept.write("transcript_id\tclass\tgene_id\tgene_name\ttx_len\tcds_start\tcds_end\n")
        for tid, gid, gname, ttype, chrom, strand, ex, cd in rows:
            inc("all_transcripts")
            inc(f"type:{ttype}")
            if cd and ttype == "protein_coding":
                c = spliced(g, chrom, strand, cd)
                if c[-3:] not in STOPS:  # annotation without the stop codon? add it if the next codon is a stop
                    nxt = (g[chrom][cd[-1][1]:cd[-1][1] + 3] if strand == "+"
                           else g[chrom][cd[0][0] - 4:cd[0][0] - 1].translate(COMP)[::-1])
                    if nxt in STOPS:
                        if strand == "+":
                            cd = cd[:-1] + [(cd[-1][0], cd[-1][1] + 3)]
                        else:
                            cd = [(cd[0][0] - 3, cd[0][1])] + cd[1:]
                        stop_added += 1
            seq = spliced(g, chrom, strand, ex)
            out.write("\t".join([tid, gid, gname, ".", ttype, ".", chrom, strand,
                                 ",".join(f"{s}-{e}" for s, e in ex), ",".join(f"{s}-{e}" for s, e in cd) or "."]) + "\n")
            allfa.write(f">{tid}\n{seq}\n")
            for s, e in ex:
                gtf.write(f'{chrom}\tTAIR12\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}";\n')
            nuclear = chrom in NUCLEAR
            if ttype in NEG_CLASSES:
                inc(f"neg_candidates:{ttype}")
                if not nuclear:
                    inc("neg_drop_organellar")
                elif overlaps_cds(chrom, strand, ex):
                    inc(f"neg_drop_overlaps_CDS:{ttype}")
                else:
                    inc(f"neg_kept:{ttype}")
                    kept.write(f"{tid}\tlncRNA\t{gid}\t{gname}\t{len(seq)}\tNA\tNA\n")
                    kfa.write(f">{tid}\n{seq}\n")
                continue
            if ttype != "protein_coding":
                continue
            inc("protein_coding")
            if not cd:
                inc("pc_drop_no_cds")
                continue
            if not nuclear:
                inc("pc_drop_organellar")
                continue
            # CDS in transcript coordinates
            exl = ex if strand == "+" else ex[::-1]

            def toff(p):
                o = 0
                for s, e in exl:
                    if s <= p <= e:
                        return o + (p - s + 1 if strand == "+" else e - p + 1)
                    o += e - s + 1
                return None
            a_, b_ = (toff(cd[0][0]), toff(cd[-1][1])) if strand == "+" else (toff(cd[-1][1]), toff(cd[0][0]))
            if a_ is None or b_ is None:
                inc("pc_drop_cds_outside_exons")
                continue
            c = seq[a_ - 1:b_]
            if len(c) != sum(e - s + 1 for s, e in cd):
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
                continue
            inc("coding_kept")
            kept.write(f"{tid}\tcoding\t{gid}\t{gname}\t{len(seq)}\t{a_}\t{b_}\n")
            kfa.write(f">{tid}\n{seq}\n")
    cnt["cds_stop_codon_added_from_genome"] = stop_added
    with open(os.path.join(outdir, "selection_counts.tsv"), "w") as fh:
        for k, v in cnt.items():
            fh.write(f"{k}\t{v}\n")
            print(k, v, sep="\t")


if __name__ == "__main__":
    main()
