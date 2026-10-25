#!/usr/bin/env python3
"""Where do assembled GM12878 transcripts begin relative to the MANE CDS? (evidence E29)

usage: python 03_utr_gm12878_truth.py <structures.gtf|.bam> <transcripts.fa> <mane.gtf.gz> <out_prefix>

For every assembled transcript (exons from the assembler's GTF, or from the primary minimap2
alignment for Trinity; strand = alignment strand, the libraries are stranded) and every MANE Select
CDS on the same strand whose stop codon lies in one of its exons, the CDS is followed from the stop
toward the start in transcript coordinates. The pair is kept when the CDS positions map to
consecutive transcript positions from the stop up to either the start codon or the transcript's 5'
end (so the transcript carries the CDS 3' part with the reference splicing):
  complete        the start codon is in the transcript; utr5 = transcript bases upstream of it;
                  open5 = no in-frame stop in that UTR; first_atg = no in-frame ATG in it
  5prime_partial  the transcript begins inside the CDS; missing_nt = CDS bases upstream of the
                  transcript's first base; missing_frac = missing_nt / CDS length
Pairs whose structure departs from the reference between stop and 5' end are 'discordant'.
Writes <out_prefix>.pairs.tsv.gz (one row per transcript and MANE CDS; first MANE CDS per transcript)
and <out_prefix>.truth.tsv in the benchmark truth format (transcript coordinates, '+' = transcript
orientation; 5' partials start at the first whole in-frame codon), for 03_utr_features.py.
"""
import bisect
import collections
import gzip
import re
import subprocess
import sys

structs, fasta, mane, outp = sys.argv[1:5]
COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def gtf_exons(path):
    ex = collections.defaultdict(list)
    info = {}
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        if len(x) < 9 or x[2] != "exon":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', x[8]).group(1)
        ex[tid].append((int(x[3]), int(x[4])))
        info[tid] = (x[0], x[6])
    out = {}
    for t, v in ex.items():
        q, blocks = 0, []
        for a, b in sorted(v):
            blocks.append((a, b, q))
            q += b - a + 1
        out[t] = (info[t][0], info[t][1], blocks)
    return out


def bam_exons(path):
    """Aligned blocks (genome start, genome end, query offset) of each primary alignment; query
    offsets in alignment orientation, counting soft clips and insertions."""
    out = {}
    p = subprocess.Popen(["samtools", "view", "-F", "0x904", path], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        x = line.split("\t", 6)
        flag, chrom, pos, cig = int(x[1]), x[2], int(x[3]), x[5]
        blocks, g, q = [], pos, 0
        for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            n = int(n)
            if o in "M=X":
                blocks.append((g, g + n - 1, q))
                g += n
                q += n
            elif o in "DN":
                g += n
            elif o in "IS":
                q += n
        out[x[0]] = (chrom, "-" if flag & 16 else "+", blocks)
    return out


def mane_cds():
    cds = collections.defaultdict(list)
    info = {}
    for line in gzip.open(mane, "rt"):
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        if x[2] not in ("CDS", "stop_codon") or "MANE_Select" not in x[8]:
            continue
        tid = re.search(r'transcript_id "([^"]+)"', x[8]).group(1)
        cds[tid].append((int(x[3]), int(x[4])))
        info[tid] = (x[0], x[6], re.search(r'gene_name "([^"]+)"', x[8]).group(1))
    out = {}
    for t, iv in cds.items():
        chrom, strand, gname = info[t]
        pos = sorted({p for s, e in iv for p in range(s, e + 1)})
        if strand == "-":
            pos = pos[::-1]
        out[t] = (chrom, strand, gname, pos)  # CDS positions 5'->3' incl. stop
    return out


tx = bam_exons(structs) if structs.endswith(".bam") else gtf_exons(structs)
seqs = read_fasta(fasta)
ref = mane_cds()
# index MANE stops by chrom/strand
stops_idx = collections.defaultdict(list)
for t, (c, s, g, pos) in ref.items():
    stops_idx[(c, s)].append((min(pos[-3:]), max(pos[-3:]), t))
for k in stops_idx:
    stops_idx[k].sort()
starts = {k: [a for a, b, t in v] for k, v in stops_idx.items()}

counts = collections.Counter()
with gzip.open(outp + ".pairs.tsv.gz", "wt") as o, open(outp + ".truth.tsv", "w") as tt:
    o.write("transcript_id\tmane_id\tgene\tstatus\ttx_len\tcds_len\tutr5\topen5\tfirst_atg\tmissing_nt\tmissing_frac\n")
    tt.write("transcript_id\tcds_start\tcds_end\tstrand\ttranscript_length\tgene_id\ttransl_table\ttype\tutr5_source\tutr3_source\n")
    for tid, (chrom, strand, blocks) in tx.items():
        seq = seqs.get(tid)
        if seq is None:
            continue
        L = len(seq)
        if blocks[-1][2] + blocks[-1][1] - blocks[-1][0] + 1 > L:
            counts["length_mismatch"] += 1
            continue
        lo, hi = blocks[0][0], blocks[-1][1]
        key = (chrom, strand)
        if key not in stops_idx:
            continue
        i0 = bisect.bisect_left(starts[key], lo - 3)
        cands = [t for a, b, t in stops_idx[key][i0:] if a <= hi and a >= lo][:50]
        if not cands:
            continue
        gs = [b[0] for b in blocks]

        def toff(p):
            j = bisect.bisect_right(gs, p) - 1
            if j < 0 or p > blocks[j][1]:
                return None
            g = blocks[j][2] + p - blocks[j][0]
            return g if strand == "+" else L - 1 - g

        for m in cands:
            pos = ref[m][3]
            n = len(pos)
            stop_off = [toff(p) for p in pos[-3:]]
            if None in stop_off or stop_off != list(range(stop_off[0], stop_off[0] + 3)):
                continue
            # walk upstream from the stop
            k, prev = n - 1, stop_off[2]
            while k > 0:
                q = toff(pos[k - 1])
                if q is None or q != prev - 1:
                    break
                prev, k = q, k - 1
            # prev = transcript offset of CDS position k (0-based in CDS)
            if k == 0:
                status = "complete"
                u = prev
                ups = [seq[i:i + 3] for i in range(u - 3, -1, -3)]
                op5 = int(not any(c in STOPS for c in ups))
                fa = int("ATG" not in [c for c in ups[:next((j for j, c in enumerate(ups) if c in STOPS), len(ups))]])
                miss = 0
            elif prev == 0:
                status, u, op5, fa, miss = "5prime_partial", 0, "", "", k
            else:
                counts["discordant"] += 1
                continue
            counts[status] += 1
            o.write(f"{tid}\t{m}\t{ref[m][2]}\t{status}\t{L}\t{n}\t{u}\t{op5}\t{fa}\t{miss}\t{miss / n:.4f}\n")
            if status == "complete":
                cs, ce = prev + 1, prev + n
            else:
                ph = (3 - k % 3) % 3  # bases before the first whole in-frame codon
                cs, ce = 1 + ph, stop_off[2] + 1
            tt.write(f"{tid}\t{cs}\t{ce}\t+\t{L}\t{m}\t1\t{status}\t.\t.\n")
            break
print(outp, dict(counts), file=sys.stderr)
