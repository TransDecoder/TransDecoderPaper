#!/usr/bin/env python3
"""Label the transcripts of a Trinity assembly of simulated reads by where the annotated CDS of
their source transcript lies in them (evidence E30).

usage: label_assembly.py --bam Trinity.minimap2.bam --fasta Trinity.fasta
                         --gencode gencode_v50.tx.tsv.gz --kept kept.tsv --salmon quant.sf
                         --out <prefix>

Inputs
  bam      minimap2 -ax splice:hq -uf --secondary=no alignment of the assembled transcripts (sorted);
           the primary alignment gives each transcript's genomic structure. Transcripts are in sense
           orientation (stranded RF assembly; -uf), so the alignment strand is the transcript strand.
  gencode  per-transcript table from 01_gencode_tables.py (all GENCODE v50 transcripts; CDS segments
           include the stop codon)
  kept     the simulated transcripts (coding with a complete CDS, and lncRNAs)
  salmon   quant.sf of the real reads (TPM of the source transcripts)

Structure compatibility (assembled transcript T vs simulated transcript R, same chromosome and strand)
  shared span = [max(T start, R start), min(T end, R end)] of the aligned / annotated extents.
  T's exonic bases inside the span (aligned blocks; short alignment deletions 'D' are bridged, 'N' gaps
  are introns) must equal R's exonic bases inside the span, and they must overlap. So the introns are
  identical across the overlap, T has no intronic sequence of R inside the overlap (no retained intron,
  no alternative splice site), and T may extend beyond either end of R.
  Exception (END_TOL = 10 nt; used only when T has no exactly compatible transcript): T's first / last
  aligned base may extend up to 10 nt past R's exon
  boundary (a terminal overhang that minimap2 aligned into the intron instead of across the junction).
  Those bases are then treated like a soft clip (below).
Soft clips / overhangs: allowed in a UTR. Where the CDS runs off T's aligned part into <= 10 clipped or
  overhanging bases, the CDS is extended in frame to the last whole codon of T (these bases are T's own
  sequence, e.g. the start of the next exon); more than 10 such bases -> unresolved.

Labels (one per assembled transcript)
  complete / 5prime_partial / 3prime_partial / internal
           T is compatible with >= 1 simulated coding transcript whose CDS it overlaps, and every such
           transcript gives the same CDS projection (start and end in T coordinates). The CDS is
           projected through the alignment: cds_start = first whole codon of the reference frame in T
           (the ATG if present), cds_end = end of the stop codon if present, else of the last whole
           codon. complete = ATG and stop in T; 5prime_partial = stop only; 3prime_partial = ATG only;
           internal = neither.
  short_partial / short_complete
           as above but fewer than 50 sense codons (stop excluded) in T, or less than one whole codon of
           the CDS in T; left out of the truth (and of the false-positive count).
  ambiguous  compatible coding transcripts give different projections; left out of scoring.
  noncoding  compatible only with simulated lncRNAs, and T's exons overlap no CDS of any GENCODE v50
           transcript on the same strand. These are the negatives (absent from the truth file, so any
           prediction on them is a false positive in score_predictions.py).
  unresolved everything else, with a reason:
           unaligned, chimeric (supplementary alignment), large_softclip (> 30 nt clipped at an end),
           softclip_at_open_cds_end (> 10 nt clipped at an end where the CDS runs off T), no_overlap (no simulated
           transcript overlaps T on its strand), incompatible (overlaps, but no compatible structure),
           utr_only (compatible with coding transcripts but overlaps none of their CDSs),
           lnc_overlaps_cds (compatible only with lncRNAs, but overlaps a GENCODE CDS),
           frame_check_failed (projected CDS is not a whole number of codons, has an internal stop, or
           lacks the ATG / stop its type requires: indels or base errors in the assembly),
           cds_end_in_deletion (a CDS end falls in an alignment deletion).
Outputs
  <out>.labels.tsv   every assembled transcript: label, reason, structure and source details
  <out>.truth.tsv    scored positives, benchmark truth format (score_predictions.py) + extra columns
  <out>.negatives.txt / <out>.excluded.txt   noncoding IDs / IDs whose predictions are not scored
  <out>.counts.tsv   label and reason counts
"""
import argparse
import bisect
import collections
import gzip
import re
import subprocess
import sys

STOPS = {"TAA", "TAG", "TGA"}
MIN_CODONS = 50
MAX_CLIP = 30
END_TOL = 10  # nt of terminal overhang tolerated at T's aligned ends (see compatibility)
BIN = 100000
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def segs(s):
    return [] if s in (".", "") else [tuple(map(int, x.split("-"))) for x in s.split(",")]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def clip(iv, lo, hi):
    out = []
    for s, e in iv:
        s2, e2 = max(s, lo), min(e, hi)
        if s2 <= e2:
            out.append((s2, e2))
    return out


def overlaps(a, b):
    """any overlap between two sorted interval lists"""
    i = j = 0
    while i < len(a) and j < len(b):
        if a[i][1] < b[j][0]:
            i += 1
        elif b[j][1] < a[i][0]:
            j += 1
        else:
            return True
    return False


class Index:
    def __init__(self):
        self.bins = collections.defaultdict(list)

    def add(self, chrom, strand, s, e, item):
        for b in range(s // BIN, e // BIN + 1):
            self.bins[(chrom, strand, b)].append(item)

    def query(self, chrom, strand, s, e):
        seen = {}
        for b in range(s // BIN, e // BIN + 1):
            for it in self.bins.get((chrom, strand, b), ()):
                seen[id(it)] = it
        return list(seen.values())


def read_bam(path):
    """primary alignment per transcript: chrom, strand, M blocks [(gs, ge, qoff)], exons (D bridged),
    5'/3' soft clips in transcript orientation; and the set of transcripts with supplementary hits."""
    aln, supp = {}, set()
    p = subprocess.Popen(["samtools", "view", "-F", "0x104", path], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        x = line.split("\t", 6)
        flag = int(x[1])
        if flag & 2048:
            supp.add(x[0])
            continue
        g, q = int(x[3]), 0
        blocks, exons = [], []
        ops = [(int(n), o) for n, o in CIG.findall(x[5])]
        lead = ops[0][0] if ops[0][1] in "SH" else 0
        trail = ops[-1][0] if ops[-1][1] in "SH" else 0
        for n, o in ops:
            if o in "M=X":
                blocks.append((g, g + n - 1, q))
                if exons and exons[-1][1] >= g - 1:
                    exons[-1][1] = g + n - 1
                else:
                    exons.append([g, g + n - 1])
                g += n
                q += n
            elif o == "D":
                if exons:
                    exons[-1][1] = g + n - 1
                g += n
            elif o == "N":
                g += n
            elif o in "IS":
                q += n
            elif o == "H":
                q += n
        strand = "-" if flag & 16 else "+"
        c5, c3 = (lead, trail) if strand == "+" else (trail, lead)
        aln[x[0]] = dict(chrom=x[2], strand=strand, blocks=blocks, exons=[tuple(e) for e in exons],
                         clip5=c5, clip3=c3)
    p.wait()
    return aln, supp


def main():
    ap = argparse.ArgumentParser()
    for a in ("bam", "fasta", "gencode", "kept", "salmon", "out"):
        ap.add_argument("--" + a, required=True)
    args = ap.parse_args()

    tpm = {}
    with open(args.salmon) as fh:
        fh.readline()
        for line in fh:
            f = line.split("\t")
            tpm[f[0]] = float(f[3])
    kept = {}
    with open(args.kept) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            kept[f[0]] = f[1]
    ref_idx, cds_idx = Index(), Index()
    refs = {}
    with gzip.open(args.gencode, "rt") as fh:
        fh.readline()
        for line in fh:
            tid, gid, gname, gtype, ttype, tags, chrom, strand, ex, cd = line.rstrip("\n").split("\t")
            cds = segs(cd)
            if cds:
                for s, e in cds:
                    cds_idx.add(chrom, strand, s, e, (s, e))
            if tid in kept:
                exs = segs(ex)
                r = dict(tid=tid, gid=gid, gname=gname, cls=kept[tid], chrom=chrom, strand=strand,
                         exons=exs, cds=cds if kept[tid] == "coding" else [],
                         s=exs[0][0], e=exs[-1][1])
                if r["cds"]:
                    r["cdslen"] = sum(e - s + 1 for s, e in cds)
                refs[tid] = r
                ref_idx.add(chrom, strand, r["s"], r["e"], r)
    print(f"refs {len(refs)}", file=sys.stderr)

    seqs = read_fasta(args.fasta)
    aln, supp = read_bam(args.bam)
    print(f"assembled {len(seqs)}, primary aligned {len(aln)}, with supplementary {len(supp)}", file=sys.stderr)

    counts = collections.Counter()
    lab = open(args.out + ".labels.tsv", "w")
    tru = open(args.out + ".truth.tsv", "w")
    neg = open(args.out + ".negatives.txt", "w")
    exc = open(args.out + ".excluded.txt", "w")
    lab.write("transcript_id\tlabel\treason\ttranscript_length\tchrom\tstrand\tstart\tend\tn_exons\t"
              "n_compat_coding\tn_compat_lnc\tsource_transcripts\tgenes\ttpm\tcds_start\tcds_end\t"
              "n_codons\tcds_frac_present\n")
    tru.write("transcript_id\tcds_start\tcds_end\tstrand\ttranscript_length\tgene_id\ttransl_table\ttype\t"
              "utr5_source\tutr3_source\tsource_transcripts\tgene_name\ttpm\tcds_frac_present\tn_codons\n")

    for tid, seq in seqs.items():
        L = len(seq)
        row = dict(label="unresolved", reason="", chrom=".", strand=".", start=".", end=".", nex=".",
                   ncc=0, ncl=0, src=".", genes=".", tpm=".", cs=".", ce=".", ncod=".", frac=".")
        a = aln.get(tid)
        if a is None:
            row["reason"] = "unaligned"
        else:
            T_ex = a["exons"]
            ts, te = T_ex[0][0], T_ex[-1][1]
            row.update(chrom=a["chrom"], strand=a["strand"], start=ts, end=te, nex=len(T_ex))
            blocks = a["blocks"]
            gs = [b[0] for b in blocks]
            qlen = blocks[-1][2] + blocks[-1][1] - blocks[-1][0] + 1 + (a["clip3"] if a["strand"] == "+" else a["clip5"])
            if qlen != L:
                sys.exit(f"query length mismatch {tid}: {qlen} vs {L}")

            def toff(p):
                j = bisect.bisect_right(gs, p) - 1
                if j < 0 or p > blocks[j][1]:
                    return None
                g = blocks[j][2] + p - blocks[j][0]
                return g if a["strand"] == "+" else L - 1 - g

            cands = [r for r in ref_idx.query(a["chrom"], a["strand"], ts, te) if r["s"] <= te and r["e"] >= ts]
            compat, fuzzy = [], {}
            for r in cands:
                lo, hi = max(ts, r["s"]), min(te, r["e"])
                ct, cr = clip(T_ex, lo, hi), clip(r["exons"], lo, hi)
                if ct and ct == cr:
                    compat.append(r)
                if r not in compat and END_TOL and ct and len(ct) == len(cr):
                    # tolerate <= END_TOL nt at T's own aligned ends (terminal overhangs that minimap2
                    # aligned into the intron instead of across the splice junction)
                    dl = cr[0][0] - ct[0][0]
                    dr = ct[-1][1] - cr[-1][1]
                    inner_ok = all(ct[i][1] == cr[i][1] for i in range(len(ct) - 1)) and \
                        all(ct[i][0] == cr[i][0] for i in range(1, len(ct)))
                    left_ok = dl == 0 or (0 < dl <= END_TOL and ct[0][0] == ts)
                    right_ok = dr == 0 or (0 < dr <= END_TOL and ct[-1][1] == te)
                    if inner_ok and left_ok and right_ok and (dl or dr):
                        compat.append(r)
                        fuzzy[r["tid"]] = (dl, dr)
            if any(r["tid"] not in fuzzy for r in compat):  # end-tolerant matches only when no exact one
                compat = [r for r in compat if r["tid"] not in fuzzy]
                fuzzy = {}
            cod = [r for r in compat if r["cls"] == "coding"]
            lnc = [r for r in compat if r["cls"] == "lncRNA"]
            row.update(ncc=len(cod), ncl=len(lnc))
            cod_cds = [r for r in cod if overlaps(T_ex, r["cds"])]
            if tid in supp:
                row["reason"] = "chimeric"
            elif max(a["clip5"], a["clip3"]) > MAX_CLIP:
                row["reason"] = "large_softclip"
            elif not cands:
                row["reason"] = "no_overlap"
            elif not compat:
                row["reason"] = "incompatible"
            elif cod_cds:
                projs = collections.defaultdict(list)
                fails = collections.Counter()
                fracs = {}
                for r in cod_cds:
                    if r["tid"] in fuzzy:  # trim the overhang; treat it like a soft clip at that end
                        dl, dr = fuzzy[r["tid"]]
                        Tt = [list(x) for x in T_ex]
                        Tt[0][0] += dl
                        Tt[-1][1] -= dr
                        Tt = [tuple(x) for x in Tt]
                        a2 = dict(a)
                        e5, e3 = (dl, dr) if a["strand"] == "+" else (dr, dl)
                        a2["clip5"], a2["clip3"] = a["clip5"] + e5, a["clip3"] + e3
                        pr = project(r, Tt, toff, a2, L)
                    else:
                        pr = project(r, T_ex, toff, a, L)
                    if isinstance(pr, str):
                        fails[pr] += 1
                        projs[("fail", pr)].append(r)
                    else:  # group by the projection in T only (cds_start, cds_end, type)
                        projs[pr[:3]].append(r)
                        fracs.setdefault(pr[:3], pr[3])
                if len(projs) > 1:
                    row["label"], row["reason"] = "ambiguous", f"{len(projs)}_projections"
                    srcs = cod_cds
                else:
                    (key, srcs), = projs.items()
                    if key[0] == "fail":
                        row["reason"] = key[1]
                        if key[1] == "cds_overlap_lt_1_codon":
                            row["label"] = "short_partial"
                    else:
                        cs, ce, typ = key
                        c = seq[cs - 1:ce]
                        has_start = typ in ("complete", "3prime_partial")
                        has_stop = typ in ("complete", "5prime_partial")
                        ok = (len(c) % 3 == 0 and (not has_start or c[:3] == "ATG") and
                              (not has_stop or c[-3:] in STOPS) and
                              not any(c[i:i + 3] in STOPS for i in range(0, len(c) - 3 * has_stop, 3)))
                        ncod = len(c) // 3 - (1 if has_stop else 0)
                        row.update(cs=cs, ce=ce, ncod=ncod, frac=f"{fracs[key]:.4f}")
                        if not ok:
                            row["reason"] = "frame_check_failed"
                        elif ncod < MIN_CODONS:
                            row["label"] = "short_complete" if typ == "complete" else "short_partial"
                            row["reason"] = typ
                        else:
                            row["label"], row["reason"] = typ, "resolved"
                row.update(src=",".join(sorted(r["tid"] for r in srcs)),
                           genes=",".join(sorted({r["gname"] for r in srcs})),
                           tpm=f"{sum(tpm.get(r['tid'], 0.0) for r in srcs):.4f}")
            elif cod:
                row["reason"] = "utr_only"
                row.update(src=",".join(sorted(r["tid"] for r in cod)),
                           genes=",".join(sorted({r["gname"] for r in cod})),
                           tpm=f"{sum(tpm.get(r['tid'], 0.0) for r in cod):.4f}")
            else:  # lncRNA only
                row.update(src=",".join(sorted(r["tid"] for r in lnc)),
                           genes=",".join(sorted({r["gname"] for r in lnc})),
                           tpm=f"{sum(tpm.get(r['tid'], 0.0) for r in lnc):.4f}")
                cdsov = [iv for iv in cds_idx.query(a["chrom"], a["strand"], ts, te)]
                if overlaps(T_ex, sorted(set(cdsov))):
                    row["reason"] = "lnc_overlaps_cds"
                else:
                    row["label"], row["reason"] = "noncoding", "lncRNA_only"
        counts[(row["label"], row["reason"] if row["label"] in ("unresolved", "ambiguous") else "")] += 1
        lab.write("\t".join(map(str, [tid, row["label"], row["reason"], L, row["chrom"], row["strand"],
                                      row["start"], row["end"], row["nex"], row["ncc"], row["ncl"],
                                      row["src"], row["genes"], row["tpm"], row["cs"], row["ce"],
                                      row["ncod"], row["frac"]])) + "\n")
        if row["label"] in ("complete", "5prime_partial", "3prime_partial", "internal"):
            typ = row["label"]
            u5 = "assembled" if typ in ("complete", "3prime_partial") else "none"
            u3 = "assembled" if typ in ("complete", "5prime_partial") else "none"
            gid = ",".join(sorted({refs[s]["gid"] for s in row["src"].split(",")}))
            tru.write("\t".join(map(str, [tid, row["cs"], row["ce"], "+", L, gid, 1, typ, u5, u3,
                                          row["src"], row["genes"], row["tpm"], row["frac"], row["ncod"]])) + "\n")
        elif row["label"] == "noncoding":
            neg.write(tid + "\n")
        else:
            exc.write(f"{tid}\t{row['label']}\t{row['reason']}\n")
    for f in (lab, tru, neg, exc):
        f.close()
    with open(args.out + ".counts.tsv", "w") as fh:
        fh.write("label\treason\tn\n")
        for (l, r), n in sorted(counts.items()):
            fh.write(f"{l}\t{r}\t{n}\n")
            print(l, r, n, sep="\t", file=sys.stderr)


def project(r, T_ex, toff, a, L):
    """Project r's CDS onto T. Returns (cds_start, cds_end, type, frac) in T coordinates (1-based,
    sense; frac = fraction of r's CDS bases inside T), or a failure reason string."""
    cds = r["cds"]  # ascending genomic segments, incl. stop
    n = r["cdslen"]
    plus = r["strand"] == "+"
    # covered CDS genomic intervals (T exons, D bridged, within the CDS)
    cov = []
    for s, e in cds:
        for x, y in T_ex:
            s2, e2 = max(s, x), min(e, y)
            if s2 <= e2:
                cov.append((s2, e2))
    cov.sort()
    glo, ghi = cov[0][0], cov[-1][1]

    def cds_off(g):  # 0-based offset of genomic position g in the CDS, 5'->3'
        o = 0
        for s, e in cds:
            if s <= g <= e:
                o += g - s
                break
            o += e - s + 1
        return o if plus else n - 1 - o

    def cds_gen(k):  # genomic position of CDS offset k (5'->3')
        k2 = k if plus else n - 1 - k
        for s, e in cds:
            if k2 <= e - s:
                return s + k2
            k2 -= e - s + 1
        raise ValueError

    g5, g3 = (glo, ghi) if plus else (ghi, glo)
    k5, k3 = cds_off(g5), cds_off(g3)
    if k3 < k5:
        return "projection_error"
    has_start, has_stop = k5 == 0, k3 == n - 1
    # soft clips where the CDS runs off T's aligned part
    t5_is_cds = (a["strand"] == "+" and glo == T_ex[0][0]) or (a["strand"] == "-" and ghi == T_ex[-1][1])
    t3_is_cds = (a["strand"] == "+" and ghi == T_ex[-1][1]) or (a["strand"] == "-" and glo == T_ex[0][0])
    open5 = not has_start and t5_is_cds and a["clip5"] > 0
    open3 = not has_stop and t3_is_cds and a["clip3"] > 0
    if (open5 and a["clip5"] > END_TOL) or (open3 and a["clip3"] > END_TOL):
        return "softclip_at_open_cds_end"
    k5w = k5 + (3 - k5 % 3) % 3
    k3w = k3 if has_stop else k3 - (k3 + 1) % 3
    if k3w - k5w + 1 < 3:
        return "cds_overlap_lt_1_codon"
    o5, o3 = toff(cds_gen(k5w)), toff(cds_gen(k3w))
    if o5 is None or o3 is None:
        return "cds_end_in_deletion"
    # CDS open at an end with <= END_TOL clipped / overhanging bases: those bases are the transcript's own
    # continuation, so extend the CDS in frame to the last whole codon in T (checked later by translation)
    if open5:
        o5 -= 3 * (o5 // 3)
    if open3:
        o3 += 3 * ((L - 1 - o3) // 3)
    typ = {(True, True): "complete", (False, True): "5prime_partial",
           (True, False): "3prime_partial", (False, False): "internal"}[(has_start, has_stop)]
    return (o5 + 1, o3 + 1, typ, (k3 - k5 + 1) / n)


if __name__ == "__main__":
    main()
