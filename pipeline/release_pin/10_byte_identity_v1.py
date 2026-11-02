#!/usr/bin/env python3
"""Byte comparison of pinned-commit runs (PyTD d3892b5, work/release_pin) with the runs of earlier
commits that the paper's tables were made from (E40).
  bench rel       vs work/bench_v3/runs/<org>/<input>/pytd_v7            (3a793ff; 108 inputs)
  bench precise   vs pytd_v7_precise, sensitive vs pytd_v7_sensitive     (3a793ff; 20 inputs each)
  bench retonly   vs pytd_v7_homology                                    (3a793ff; 84 inputs)
  gm rel/precise/sensitive vs work/gm12878_td/<asm>/final_v7_abinitio[_..] (3a793ff)
  gm trinity retonly vs final_v7_homology                                (3a793ff)
Files: final GFF3, BED, CDS, pep; for the benchmark also preds.tsv and LongOrfs GFF3.
usage: 10_byte_identity_v1.py > data/release_pin/byte_identity.tsv
"""
import filecmp
import os

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
NEW, OLD = f"{P}/work/release_pin/bench", f"{P}/work/bench_v3/runs"
FILES = ["input.fa.transdecoder.gff3", "input.fa.transdecoder.bed", "input.fa.transdecoder.cds",
         "input.fa.transdecoder.pep", "preds.tsv", "input.fa.transdecoder_dir/longest_orfs.gff3"]
GMF = ["transcripts.fa.transdecoder.gff3", "transcripts.fa.transdecoder.bed", "transcripts.fa.transdecoder.cds",
       "transcripts.fa.transdecoder.pep", "transcripts.fa.transdecoder.genome.gff3"]
print("set\tasm_or_org\tinput\tlabel\treference\tn_files\tn_identical\tdiffering")


def cmp(a, b, files):
    diff = []
    for f in files:
        if not (os.path.exists(f"{a}/{f}") and os.path.exists(f"{b}/{f}")) or \
                not filecmp.cmp(f"{a}/{f}", f"{b}/{f}", shallow=False):
            diff.append(os.path.basename(f))
    return diff


for org in sorted(os.listdir(NEW)):
    for inp in sorted(os.listdir(f"{NEW}/{org}")):
        for lab, ref in [("rel", "pytd_v7"), ("precise", "pytd_v7_precise"), ("sensitive", "pytd_v7_sensitive"),
                         ("retonly", "pytd_v7_homology")]:
            a = f"{NEW}/{org}/{inp}/{lab}"
            if os.path.exists(f"{a}/DONE"):
                d = cmp(a, f"{OLD}/{org}/{inp}/{ref}", FILES)
                print(f"bench\t{org}\t{inp}\t{lab}\t{ref}\t{len(FILES)}\t{len(FILES) - len(d)}\t{','.join(d) or '.'}")
GM_NEW, GM_OLD = f"{P}/work/release_pin/gm", f"{P}/work/gm12878_td"
for asm in ("lraa", "stringtie", "trinity"):
    for lab, ref in [("rel", "final_v7_abinitio"), ("precise", "final_v7_abinitio_precise"),
                     ("sensitive", "final_v7_abinitio_sensitive"), ("retonly", "final_v7_homology"),
                     ("a0.4_w17", "final_v7_abinitio")]:
        a, b = f"{GM_NEW}/{asm}/{lab}", f"{GM_OLD}/{asm}/{ref}"
        if os.path.exists(f"{a}/DONE") and os.path.isdir(b):
            d = cmp(a, b, GMF)
            print(f"gm\t{asm}\t.\t{lab}\t{ref}\t{len(GMF)}\t{len(GMF) - len(d)}\t{','.join(d) or '.'}")
