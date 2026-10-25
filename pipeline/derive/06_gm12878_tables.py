#!/usr/bin/env python3
"""Compact GM12878 tables for notebooks/06_gm12878_application.Rmd (evidence E03, E25).

usage: 06_gm12878_tables.py <work_dir=$P/work/gm12878_td> <MANE gtf.gz> <out_prefix>

For the earlier runs (final_abinitio, PyTD 48c870c) and the release runs (final_v7_abinitio, PyTD 3a793ff)
of LRAA, StringTie and Trinity GG, writes <out_prefix><name>:
  gm12878_class_by_length.tsv  run, assembly, length_bin, n, exact, compatible, other_overlap, novel
      Distinct predicted CDSs by predicted protein length (pred_status.tsv `aa`) and GENCODE class:
      exact > compatible > other_overlap (overlaps a GENCODE CDS) > novel (no same-strand GENCODE CDS
      overlap). Same as figures/fig6_gm12878/derive_class_by_length.py. <100 includes the few <50 aa.
  gm12878_genome_map.tsv       run, assembly, transcript_orfs, projected, skipped_unaligned_cds,
      skipped_no_structure  (from each run's genome_map.log)
  gm12878_mane_not_exact.tsv   run, assembly, assembled, exact_and_assembled, assembled_not_exact,
      same_stop, same_stop_start_upstream, same_stop_start_downstream
      For MANE Select CDSs that were assembled but not predicted exactly (ref_status.tsv), whether some
      distinct predicted CDS that ends in a stop codon (complete or 5' partial) has the same stop codon
      (same chromosome, strand and 3' end) and so differs only at the start. Upstream/downstream: the
      predicted start lies 5' / 3' of the MANE start.
"""
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "genome_annotation/scripts"))
from eval_genome_cds import read_pred_gff3, read_ref_gtf  # noqa: E402  (needs PyTD on PYTHONPATH)

WORK, MANE_GTF, OUT = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
RUNS = {"earlier": ("final_abinitio", "eval/status"), "release": ("final_v7_abinitio", "eval")}
ASM = ["lraa", "stringtie", "trinity"]
BINS = [("lt100", 0, 100), ("100-149", 100, 150), ("150-199", 150, 200), ("200-299", 200, 300),
        ("ge300", 300, 10 ** 9)]
CLASSES = ["exact", "compatible", "other_overlap", "novel"]


def cls(r):
    if r["gencode_exact"] == "1":
        return "exact"
    if r["gencode_compatible"] == "1":
        return "compatible"
    return "novel" if r["no_gencode_cds_overlap"] == "1" else "other_overlap"


def ends(x):
    """(start, stop-end) genomic coordinates of a CDS in transcript orientation."""
    s = x["segs"]
    return (s[0][0], s[-1][1]) if x["strand"] == "+" else (s[-1][1], s[0][0])


mane, _ = read_ref_gtf(MANE_GTF, None)
mane = {t: r for t, r in mane.items() if "MANE_Select" in r["tags"]}

with open(f"{OUT}gm12878_class_by_length.tsv", "w") as fl, \
        open(f"{OUT}gm12878_genome_map.tsv", "w") as fg, \
        open(f"{OUT}gm12878_mane_not_exact.tsv", "w") as fm:
    fl.write("run\tassembly\tlength_bin\tn\t" + "\t".join(CLASSES) + "\n")
    fg.write("run\tassembly\ttranscript_orfs\tprojected\tskipped_unaligned_cds\tskipped_no_structure\n")
    fm.write("run\tassembly\tassembled\texact_and_assembled\tassembled_not_exact\tsame_stop\t"
             "same_stop_start_upstream\tsame_stop_start_downstream\n")
    for run, (rdir, sdir) in RUNS.items():
        for asm in ASM:
            # class by length
            counts = {b[0]: Counter() for b in BINS}
            with open(WORK / sdir / f"{asm}_{rdir}.pred_status.tsv") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    aa = int(r["aa"])
                    counts[next(n for n, lo, hi in BINS if lo <= aa < hi)][cls(r)] += 1
            for name, _, _ in BINS:
                c = counts[name]
                fl.write(f"{run}\t{asm}\t{name}\t{sum(c.values())}\t" + "\t".join(str(c[k]) for k in CLASSES) + "\n")
            # genome-map log
            log = (WORK / asm / rdir / "genome_map.log").read_text()
            m = re.search(r"(\d+) / (\d+) transcript ORFs propagated", log)
            fg.write(f"{run}\t{asm}\t{m.group(2)}\t{m.group(1)}\t"
                     f"{len(re.findall(r'WARNING.*not entirely in its aligned part', log))}\t"
                     f"{len(re.findall(r'WARNING.*has no genome structure', log))}\n")
            # assembled MANE CDSs not predicted exactly: same stop codon?
            stops = defaultdict(list)
            for p in read_pred_gff3(WORK / asm / rdir / "transcripts.fa.transdecoder.genome.gff3").values():
                if p["orf_type"] in ("complete", "5prime_partial"):
                    st, en = ends(p)
                    stops[(p["chrom"], p["strand"], en)].append(st)
            n_asm = n_ex = n_ne = same = up = down = 0
            with open(WORK / sdir / f"{asm}_{rdir}.ref_status.tsv") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r["set"] != "MANE" or r["assembled"] != "1":
                        continue
                    n_asm += 1
                    if r["exact"] == "1":
                        n_ex += 1
                        continue
                    n_ne += 1
                    ref = mane[r["transcript_id"]]
                    rst, ren = ends(ref)
                    starts = stops.get((ref["chrom"], ref["strand"], ren), [])
                    if starts:
                        same += 1
                        # nearest predicted start to the MANE start decides up/downstream
                        s = min(starts, key=lambda x: abs(x - rst))
                        upstream = (s < rst) if ref["strand"] == "+" else (s > rst)
                        up += upstream
                        down += not upstream
            fm.write(f"{run}\t{asm}\t{n_asm}\t{n_ex}\t{n_ne}\t{same}\t{up}\t{down}\n")
            print(run, asm, "done", file=sys.stderr)
