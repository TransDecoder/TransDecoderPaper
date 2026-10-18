#!/usr/bin/env python3
"""Build a coding-region truth set (one transcript per gene) from a genome + GFF3/GTF.

Replaces build_truth_from_gbff.py so that each organism can use its most current
reference annotation (GENCODE/MANE, FlyBase, TAIR, PomBase, TGD, RefSeq), most of which
are distributed as genome + GFF3/GTF rather than RefSeq transcript records.

Parsing
  GFF3: transcripts are features (type --transcript-types) that are the Parent of CDS
        features; gene = the transcript's Parent. Exons are the exon features, or, where
        a source has none (PomBase, TGD), the merged CDS + UTR pieces.
  GTF:  transcripts are grouped by transcript_id, gene by gene_id.
  stop_codon features are merged into the CDS (GTF and GENCODE CDS exclude the stop).
  If the CDS still does not end in a stop but the next codon in the transcript is one,
  the CDS is extended by that codon (counted as 'stop_appended').

A transcript is eligible when it:
  - is on a nuclear sequence (not matching --exclude-seqids) and, if --biotype-attr is
    given, has that attribute equal to protein_coding,
  - carries the --require-tag tag if given (e.g. MANE_Select, Ensembl_canonical), and
    none of the --exclude-tags (e.g. GENCODE cds_start_NF),
  - has a CDS contained in its exons, contiguous in transcript coordinates,
  - has transl_table (when annotated on the CDS) equal to --table,
  - starts with ATG, ends with a stop of --table, has length divisible by 3, no internal
    stops and length >= --min-cds; the transcript is ACGT only.
Per gene the transcript with the longest CDS is kept (ties: longer transcript, then id).

Missing UTRs (--add-utr)
  Transcripts with no 5' (or 3') UTR get one: a length drawn from this organism's
  empirical 5' (3') UTR lengths among the kept transcripts that have one, filled with the
  genomic sequence flanking the transcript (unspliced; truncated at a sequence end or at
  the first non-ACGT base). The utr5_source / utr3_source columns record
  annotated / added / none. Randomness comes from --seed and does not depend on --max-utr.

UTRs longer than --max-utr are then trimmed (the 2018 benchmark used 250; -1 keeps all).

Outputs (prefix --out):
  <out>.fasta      transcripts, header '>transcript_id'
  <out>.truth.tsv  transcript_id, cds_start, cds_end (1-based, incl. stop), strand,
                   transcript_length, gene_id, transl_table, type, utr5_source, utr3_source
  <out>.log        counts of exclusions by reason
"""
import argparse
import collections
import gzip
import random
import re
import sys
import urllib.parse

from Bio.Data import CodonTable

COMP = str.maketrans("ACGTN", "TGCAN")
DEFAULT_EXCLUDE_SEQIDS = (r"(?i)^(chr)?(M|MT|Mito|mitochondrion.*|Pt|C|chloroplast|plastid|apicoplast|"
                          r"NC_037526\.1|NC_036769\.1|NC_001638\.1|NC_005353\.1)$")


def revcomp(s):
    return s.translate(COMP)[::-1]


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_genome(path):
    seqs, name, buf = {}, None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def parse_attrs(field, fmt):
    attrs = collections.defaultdict(list)
    if fmt == "gtf":
        for m in re.finditer(r'(\S+)\s+"([^"]*)"', field):
            attrs[m.group(1)].append(m.group(2))
    else:
        for kv in field.strip().strip(";").split(";"):
            if "=" in kv:
                k, v = kv.strip().split("=", 1)
                attrs[k].extend(urllib.parse.unquote(x) for x in v.split(","))
    return attrs


class Tx:
    __slots__ = ("id", "gene", "seqid", "strand", "attrs", "exons", "cds", "utr", "stops", "cds_table", "declared")

    def __init__(self, tid):
        self.id, self.gene, self.seqid, self.strand = tid, None, None, None
        self.attrs = collections.defaultdict(list)
        self.exons, self.cds, self.utr, self.stops = [], [], [], []
        self.cds_table = None
        self.declared = False  # GFF3: seen as a feature of a --transcript-types type


def load_annotation(path, fmt, tx_types):
    txs, parent_of, rows = {}, {}, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            x = line.rstrip("\n").split("\t")
            if len(x) != 9:
                continue
            a = parse_attrs(x[8], fmt)
            rows.append((x[0], x[2], int(x[3]), int(x[4]), x[6], a))
            if fmt == "gff3" and "ID" in a and "Parent" in a:
                parent_of[a["ID"][0]] = a["Parent"][0]

    def get(tid):
        t = txs.get(tid)
        if t is None:
            t = txs[tid] = Tx(tid)
        return t

    for seqid, ftype, s, e, strand, a in rows:
        if fmt == "gtf":
            if "transcript_id" not in a:
                continue
            parents = a["transcript_id"][:1]
        elif ftype in tx_types and "ID" in a:
            t = get(a["ID"][0])
            t.gene = a.get("Parent", [a["ID"][0]])[0]
            t.declared = True
            t.seqid, t.strand = seqid, strand
            for k, v in a.items():
                t.attrs[k].extend(v)
            continue
        else:
            parents = a.get("Parent", [])
        for p in parents:
            t = get(p)
            t.seqid, t.strand = seqid, strand
            if fmt == "gtf":
                t.gene = a["gene_id"][0]
                if ftype == "transcript":
                    for k, v in a.items():
                        t.attrs[k].extend(v)
            if ftype == "exon":
                t.exons.append((s, e))
            elif ftype == "CDS":
                t.cds.append((s, e))
                if "transl_table" in a:
                    t.cds_table = int(a["transl_table"][0])
            elif ftype == "stop_codon":
                t.stops.append((s, e))
            elif ftype in ("five_prime_UTR", "three_prime_UTR", "UTR", "5UTR", "3UTR"):
                t.utr.append((s, e))
            if fmt == "gtf" and ftype in ("exon", "CDS") and not t.attrs:
                for k, v in a.items():
                    if k not in ("exon_number", "exon_id"):
                        t.attrs[k].extend(v)
    # GTF transcripts carry tags on every line; keep transcript-level attrs when present
    # GFF3: CDS parents of other types (e.g. PomBase pseudogenic_transcript) are dropped
    return {k: t for k, t in txs.items() if t.cds and (fmt == "gtf" or t.declared)}


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def tx_to_record(t, genome):
    """Return (seq, cds_start0, cds_end0_exclusive) in transcript orientation, or a reason string."""
    chrom = genome.get(t.seqid)
    if chrom is None:
        return "seqid_not_in_genome"
    cds = merge(t.cds + t.stops)
    exons = merge(t.exons) if t.exons else merge(t.cds + t.utr + t.stops)
    exons = merge(exons + cds)  # stop codons are sometimes annotated past the last exon end
    for cs, ce in cds:
        if not any(es <= cs and ce <= ee for es, ee in exons):
            return "cds_outside_exons"
    seq = "".join(chrom[s - 1:e] for s, e in exons)
    # transcript offset of genomic position p (1-based)
    def offset(p):
        o = 0
        for s, e in exons:
            if s <= p <= e:
                return o + p - s
            o += e - s + 1
        return None
    gs, ge = cds[0][0], cds[-1][1]
    a, b = offset(gs), offset(ge)
    if t.strand == "-":
        L = len(seq)
        seq = revcomp(seq)
        a, b = L - 1 - b, L - 1 - a
    if b - a + 1 != sum(e - s + 1 for s, e in cds):
        return "cds_not_contiguous_in_transcript"
    return seq, a, b + 1


def flank(chrom, t, exons_span, side, length):
    """Genomic flank 5' ('up') or 3' ('down') of the transcript, in transcript orientation."""
    lo, hi = exons_span
    if (side == "up") == (t.strand == "+"):
        s = chrom[max(0, lo - 1 - length):lo - 1]
        seq = s if t.strand == "+" else revcomp(s)
        # keep the part adjacent to the transcript up to the first non-ACGT base
        if t.strand == "+":
            m = re.search(r"[ACGT]*$", seq)
            return m.group(0)
        m = re.match(r"[ACGT]*", seq)
        return m.group(0)
    s = chrom[hi:hi + length]
    seq = s if t.strand == "+" else revcomp(s)
    if t.strand == "+":
        return re.match(r"[ACGT]*", seq).group(0)
    return re.search(r"[ACGT]*$", seq).group(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--format", choices=["gff3", "gtf"], help="default: from the file name")
    ap.add_argument("--out", required=True)
    ap.add_argument("--table", type=int, default=1, help="NCBI translation table")
    ap.add_argument("--transcript-types", default="mRNA,transcript")
    ap.add_argument("--biotype-attr", help="attribute that must equal protein_coding (e.g. transcript_type)")
    ap.add_argument("--require-tag", help="keep only transcripts with this value in their tag attribute")
    ap.add_argument("--exclude-tags", default="", help="comma list of tag values to exclude")
    ap.add_argument("--exclude-seqids", default=DEFAULT_EXCLUDE_SEQIDS, help="regex of organelle sequence ids")
    ap.add_argument("--min-cds", type=int, default=0)
    ap.add_argument("--add-utr", action="store_true")
    ap.add_argument("--max-utr", type=int, default=250)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    fmt = args.format or ("gtf" if re.search(r"\.gtf(\.gz)?$", args.annotation) else "gff3")
    excl_tags = set(filter(None, args.exclude_tags.split(",")))
    excl_seq = re.compile(args.exclude_seqids)
    stops = set(CodonTable.unambiguous_dna_by_id[args.table].stop_codons)
    reasons = collections.Counter()

    genome = read_genome(args.genome)
    txs = load_annotation(args.annotation, fmt, set(args.transcript_types.split(",")))
    reasons["transcripts_with_cds"] = len(txs)

    best = {}
    for t in txs.values():
        tags = set(t.attrs.get("tag", []))
        if excl_seq.match(t.seqid or ""):
            reasons["organelle_seqid"] += 1
            continue
        if args.biotype_attr and t.attrs.get(args.biotype_attr, ["protein_coding"])[0] != "protein_coding":
            reasons["not_protein_coding"] += 1
            continue
        if args.require_tag and args.require_tag not in tags:
            reasons["missing_required_tag"] += 1
            continue
        if tags & excl_tags:
            reasons["excluded_tag"] += 1
            continue
        if t.cds_table is not None and t.cds_table != args.table:
            reasons["other_transl_table"] += 1
            continue
        rec = tx_to_record(t, genome)
        if isinstance(rec, str):
            reasons[rec] += 1
            continue
        seq, s, e = rec
        if set(seq) - set("ACGT"):
            reasons["non_ACGT"] += 1
            continue
        if seq[e - 3:e] not in stops and e + 3 <= len(seq) and seq[e:e + 3] in stops and (e - s) % 3 == 0:
            e += 3
            reasons["stop_appended"] += 1
        cds = seq[s:e]
        if len(cds) % 3 or len(cds) < max(6, args.min_cds):
            reasons["cds_len"] += 1
            continue
        if not cds.startswith("ATG"):
            reasons["non_ATG_start"] += 1
            continue
        if cds[-3:] not in stops:
            reasons["no_terminal_stop"] += 1
            continue
        if any(cds[i:i + 3] in stops for i in range(0, len(cds) - 3, 3)):
            reasons["internal_stop"] += 1
            continue
        g = t.gene or t.id
        key = (e - s, len(seq), t.id)
        cur = best.get(g)
        if cur is None or key > cur["key"]:
            if cur is not None:
                reasons["not_longest_for_gene"] += 1
            best[g] = dict(t=t, seq=seq, s=s, e=e, key=key)
        else:
            reasons["not_longest_for_gene"] += 1

    # empirical UTR lengths from transcripts that have them
    u5 = sorted(c["s"] for c in best.values() if c["s"] > 0)
    u3 = sorted(len(c["seq"]) - c["e"] for c in best.values() if c["e"] < len(c["seq"]))
    rng = random.Random(args.seed)
    with open(args.out + ".fasta", "w") as fa, open(args.out + ".truth.tsv", "w") as tt:
        tt.write("transcript_id\tcds_start\tcds_end\tstrand\ttranscript_length\tgene_id\ttransl_table\ttype"
                 "\tutr5_source\tutr3_source\n")
        for g in sorted(best, key=lambda k: best[k]["t"].id):
            c = best[g]
            t, seq, s, e = c["t"], c["seq"], c["s"], c["e"]
            src5 = "annotated" if s > 0 else "none"
            src3 = "annotated" if e < len(seq) else "none"
            if args.add_utr:
                # draw both lengths for every transcript so the stream does not depend on which are missing
                l5, l3 = rng.choice(u5), rng.choice(u3)
                chrom = genome[t.seqid]
                exons = merge(t.exons) if t.exons else merge(t.cds + t.utr + t.stops)
                span = (exons[0][0], exons[-1][1])
                if src5 == "none":
                    f5 = flank(chrom, t, span, "up", l5)
                    if f5:
                        seq, s, e, src5 = f5 + seq, s + len(f5), e + len(f5), "added"
                        reasons["utr5_added"] += 1
                if src3 == "none":
                    f3 = flank(chrom, t, span, "down", l3)
                    if f3:
                        seq, src3 = seq + f3, "added"
                        reasons["utr3_added"] += 1
            if args.max_utr >= 0:
                lo, hi = max(0, s - args.max_utr), min(len(seq), e + args.max_utr)
                seq, s, e = seq[lo:hi], s - lo, e - lo
            fa.write(f">{t.id}\n{seq}\n")
            tt.write(f"{t.id}\t{s + 1}\t{e}\t+\t{len(seq)}\t{g}\t{args.table}\tcomplete\t{src5}\t{src3}\n")
    reasons["kept_genes"] = len(best)
    reasons["kept_with_annotated_utr5"] = len(u5)
    reasons["kept_with_annotated_utr3"] = len(u3)
    with open(args.out + ".log", "w") as lg:
        for k, v in reasons.most_common():
            lg.write(f"{k}\t{v}\n")
    sys.stderr.write("".join(f"{k}\t{v}\n" for k, v in reasons.most_common()))


if __name__ == "__main__":
    sys.exit(main())
