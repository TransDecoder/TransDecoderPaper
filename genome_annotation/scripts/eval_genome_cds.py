#!/usr/bin/env python3
"""Evaluate genome CDS models (TransDecoder genome-map GFF3) against MANE and GENCODE.

usage: eval_genome_cds.py --pred <genome.gff3> --structures <assembly .gtf/.bam>
                          --mane <MANE gtf[.gz]> --gencode <GENCODE gtf[.gz]> --label <name> -o <out_prefix>
       (PYTHONPATH must include a PyTransDecoder tree, for pytransdecoder.genome)

A CDS is its genomic segments including the stop codon (GTF CDS + stop_codon merged).

Match levels between a predicted CDS P and a reference CDS R (same chromosome and strand):
  exact       identical segments
  compatible  P is an in-frame, contiguous piece of R: P's segment boundaries inside R's span are R's
              intron boundaries, and P's codons are R's codons (a correct partial CDS)

Reference-centric (MANE Select coding transcripts on the primary assembly; GENCODE protein-coding transcripts with complete CDS):
  assembled   an assembled transcript contains R's whole CDS with all of R's CDS introns
              (upper bound for any CDS caller on this assembly)
  exact       some prediction matches R exactly
  compatible  some prediction (complete or partial) is compatible with R
Prediction-centric, over distinct predicted CDSs:
  exact match to a MANE CDS (among predictions overlapping a MANE CDS on the same strand);
  exact / compatible with any GENCODE CDS; no same-strand overlap with any GENCODE CDS.

Outputs: <out_prefix>.summary.tsv (metric, value), <out_prefix>.ref_status.tsv
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

from pytransdecoder.genome import read_structures

BIN = 100_000
TAG_RE = re.compile(r'tag "([^"]+)"')


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def merge(segs):
    out = []
    for a, b in sorted(segs):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return tuple((a, b) for a, b in out)


def attr(text, key):
    m = re.search(key + r' "([^"]+)"', text)
    return m.group(1) if m else ""


def read_ref_gtf(path, seqids):
    """transcript_id -> dict(chrom, strand, segs, phase5, gene, gene_name, ttype, tags)"""
    tx = {}
    for line in opener(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] not in ("CDS", "stop_codon"):
            continue
        tid = attr(f[8], "transcript_id")
        t = tx.get(tid)
        if t is None:
            t = tx[tid] = dict(chrom=f[0], strand=f[6], segs=[], phases={}, gene=attr(f[8], "gene_id"),
                               gene_name=attr(f[8], "gene_name"), ttype=attr(f[8], "transcript_type"),
                               tags=set(TAG_RE.findall(f[8])))
        a, b = int(f[3]), int(f[4])
        t["segs"].append((a, b))
        if f[2] == "CDS":
            t["phases"][(a, b)] = int(f[7]) if f[7] in "012" else 0
    refs, skipped = {}, 0
    for tid, t in tx.items():
        if seqids is not None and t["chrom"] not in seqids:
            skipped += 1
            continue
        cds = [s for s in t["phases"]]
        five = min(cds) if t["strand"] == "+" else max(cds)
        t["phase5"] = t["phases"][five]
        t["segs"] = merge(t["segs"])
        del t["phases"]
        refs[tid] = t
    return refs, skipped


def read_pred_gff3(path):
    """model_id -> dict(chrom, strand, segs, phase5, orf_type)"""
    models = {}
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] == "mRNA":
            mid = re.search(r"ID=([^;]+)", f[8]).group(1)
            m = re.search(r"ORF type:(\w+)", f[8])
            models[mid] = dict(chrom=f[0], strand=f[6], segs=[], phases={}, orf_type=m.group(1) if m else "NA")
        elif f[2] == "CDS":
            mid = re.search(r"Parent=([^;]+)", f[8]).group(1)
            a, b = int(f[3]), int(f[4])
            models[mid]["segs"].append((a, b))
            models[mid]["phases"][(a, b)] = int(f[7])
    for m in models.values():
        five = min(m["segs"]) if m["strand"] == "+" else max(m["segs"])
        m["phase5"] = m["phases"][five]
        m["segs"] = merge(m["segs"])
        del m["phases"]
    return {k: v for k, v in models.items() if v["segs"]}


def contained(sub, sup):
    """Index j of sup's segment holding sub's first segment if sub lies inside sup with sub's
    internal boundaries equal to sup's (consecutive segments), else None."""
    k = len(sub)
    for j, (a, b) in enumerate(sup):
        if a <= sub[0][0] <= b:
            break
    else:
        return None
    if j + k > len(sup):
        return None
    if k == 1:
        return j if sub[0][1] <= sup[j][1] else None
    if sub[0][1] != sup[j][1]:
        return None
    for i in range(1, k - 1):
        if sub[i] != sup[j + i]:
            return None
    last, s = sub[-1], sup[j + k - 1]
    return j if last[0] == s[0] and last[1] <= s[1] else None


def compatible(p, r):
    j = contained(p["segs"], r["segs"])
    if j is None:
        return False
    rs, ps = r["segs"], p["segs"]
    if r["strand"] == "+":
        offset = sum(b - a + 1 for a, b in rs[:j]) + ps[0][0] - rs[j][0]
    else:
        idx = j + len(ps) - 1
        offset = sum(b - a + 1 for a, b in rs[idx + 1:]) + rs[idx][1] - ps[-1][1]
    return (offset + p["phase5"]) % 3 == r["phase5"] % 3


class Index:
    def __init__(self, items):
        self.bins = defaultdict(list)
        for key, it in items.items():
            lo, hi = it["segs"][0][0], it["segs"][-1][1]
            for b in range(lo // BIN, hi // BIN + 1):
                self.bins[(it["chrom"], it["strand"], b)].append((lo, hi, key))

    def overlapping(self, it):
        lo, hi = it["segs"][0][0], it["segs"][-1][1]
        seen = set()
        for b in range(lo // BIN, hi // BIN + 1):
            for a, z, key in self.bins.get((it["chrom"], it["strand"], b), ()):
                if a <= hi and lo <= z and key not in seen:
                    seen.add(key)
                    yield key


def cds_overlap(p, r):
    i = j = 0
    while i < len(p) and j < len(r):
        if p[i][1] < r[j][0]:
            i += 1
        elif r[j][1] < p[i][0]:
            j += 1
        else:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pred", required=True)
    ap.add_argument("--structures", required=True)
    ap.add_argument("--mane", required=True)
    ap.add_argument("--gencode", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("-o", "--out-prefix", required=True)
    args = ap.parse_args()

    structs = read_structures(Path(args.structures))
    assembly = {}
    for tid, s in structs.items():
        if s.strand in "+-" and s.exons:
            assembly[tid] = dict(chrom=s.seqid, strand=s.strand, segs=tuple(sorted(s.exons)))
    seqids = {s["chrom"] for s in assembly.values()}
    print(f"[{args.label}] {len(assembly)} assembled transcripts", file=sys.stderr)

    gencode, _ = read_ref_gtf(args.gencode, None)
    mane, mane_skipped = read_ref_gtf(args.mane, {r["chrom"] for r in gencode.values()})
    mane = {t: r for t, r in mane.items() if "MANE_Select" in r["tags"]}   # not MANE Plus Clinical
    preds = read_pred_gff3(args.pred)
    print(f"[{args.label}] {len(preds)} predicted CDS models; MANE {len(mane)} "
          f"(+{mane_skipped} off primary assembly); GENCODE CDS transcripts {len(gencode)}", file=sys.stderr)

    # distinct predicted CDSs (several transcripts can give the same CDS); complete wins the type
    distinct = {}
    for mid, p in preds.items():
        key = (p["chrom"], p["strand"], p["segs"])
        if key not in distinct or p["orf_type"] == "complete":
            distinct[key] = dict(p, model_id=mid)
    pred_by_key = {f"p{i}": p for i, p in enumerate(distinct.values())}
    exact_keys = {(p["chrom"], p["strand"], p["segs"]) for p in pred_by_key.values()}

    asm_index, pred_index = Index(assembly), Index(pred_by_key)
    gencode_full = {t: r for t, r in gencode.items() if r["ttype"] == "protein_coding"
                    and not ({"cds_start_NF", "cds_end_NF"} & r["tags"])}

    summary = [("label", args.label), ("assembled_transcripts", len(assembly)),
               ("predicted_cds_models", len(preds)), ("distinct_predicted_cds", len(pred_by_key))]
    for t in ("complete", "5prime_partial", "3prime_partial", "internal"):
        summary.append((f"distinct_predicted_{t}", sum(p["orf_type"] == t for p in pred_by_key.values())))

    status_rows = []
    for setname, refs in (("MANE", mane), ("GENCODE_complete", gencode_full)):
        n = asm = ex = comp = ex_given_asm = 0
        for tid, r in refs.items():
            n += 1
            in_asm = any(assembly[k]["segs"] and contained(r["segs"], assembly[k]["segs"]) is not None
                         for k in asm_index.overlapping(r))
            is_exact = (r["chrom"], r["strand"], r["segs"]) in exact_keys
            is_comp = is_exact or any(compatible(pred_by_key[k], r) for k in pred_index.overlapping(r))
            asm += in_asm
            ex += is_exact
            comp += is_comp
            ex_given_asm += is_exact and in_asm
            status_rows.append((setname, tid, r["gene_name"], int(in_asm), int(is_exact), int(is_comp)))
        summary += [(f"{setname}_n", n), (f"{setname}_assembled", asm),
                    (f"{setname}_exact", ex), (f"{setname}_compatible", comp),
                    (f"{setname}_sens_exact", f"{ex / n:.4f}"),
                    (f"{setname}_sens_compatible", f"{comp / n:.4f}"),
                    (f"{setname}_assembled_frac", f"{asm / n:.4f}"),
                    (f"{setname}_sens_exact_given_assembled", f"{ex_given_asm / max(asm, 1):.4f}")]

    # prediction-centric
    mane_index, gc_index = Index(mane), Index(gencode)
    mane_keys = {(r["chrom"], r["strand"], r["segs"]) for r in mane.values()}
    gc_keys = {(r["chrom"], r["strand"], r["segs"]) for r in gencode.values()}
    for subset in ("complete", "partial", "all"):
        ps = [p for p in pred_by_key.values()
              if subset == "all" or (p["orf_type"] == "complete") == (subset == "complete")]
        at_mane = mane_hit = gc_exact = gc_comp = novel = 0
        for p in ps:
            key = (p["chrom"], p["strand"], p["segs"])
            if any(cds_overlap(p["segs"], mane[k]["segs"]) for k in mane_index.overlapping(p)):
                at_mane += 1
                mane_hit += key in mane_keys
            ov = [k for k in gc_index.overlapping(p) if cds_overlap(p["segs"], gencode[k]["segs"])]
            if not ov:
                novel += 1
            if key in gc_keys:
                gc_exact += 1
                gc_comp += 1
            elif any(compatible(p, gencode[k]) for k in ov):
                gc_comp += 1
        m = max(len(ps), 1)
        summary += [(f"pred_{subset}_n", len(ps)),
                    (f"pred_{subset}_overlapping_MANE", at_mane),
                    (f"pred_{subset}_MANE_precision_exact", f"{mane_hit / max(at_mane, 1):.4f}"),
                    (f"pred_{subset}_GENCODE_exact_frac", f"{gc_exact / m:.4f}"),
                    (f"pred_{subset}_GENCODE_compatible_frac", f"{gc_comp / m:.4f}"),
                    (f"pred_{subset}_no_GENCODE_CDS_overlap_frac", f"{novel / m:.4f}")]

    with open(f"{args.out_prefix}.summary.tsv", "w") as fh:
        for k, v in summary:
            fh.write(f"{k}\t{v}\n")
    with open(f"{args.out_prefix}.ref_status.tsv", "w") as fh:
        fh.write("set\ttranscript_id\tgene_name\tassembled\texact\tcompatible\n")
        for row in status_rows:
            fh.write("\t".join(map(str, row)) + "\n")
    for k, v in summary:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
