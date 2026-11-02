#!/usr/bin/env python3
"""Ribo-seq and MS support of the GM12878 ORFs that the q-value test removes, at the pinned release code (E40).
Baseline: release code with the q-value test off (m50_noq); test: release default (rel). ORFs are identified by
transcript, strand and stop. Rows per assembly: all (default predictions), kept, removed (baseline ORFs with no default
prediction with that stop), added. Otherwise as pipeline/derive/08_rule_support_v7.py (kept).
usage: 14_rule_support_v1.py > data/gm12878/release_pin/rule_support.tsv
"""
import csv
import os
import re

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
G, W = f"{P}/work/release_pin/gm", f"{P}/work/release_pin/gm_scoring"
MIN_PSITES = 10


def stops(gff3):
    out = {}
    for line in open(gff3):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        out[re.search(r"Parent=([^;]+)", f[8]).group(1)] = (f[0], f[6], int(f[4]) if f[6] == "+" else int(f[3]))
    return out


def support(stop_of, tag):
    sup = {}
    for r in csv.DictReader(open(f"{W}/ribo/gm.{tag}.orf_support.tsv"), delimiter="\t"):
        if r["model_id"] in stop_of:
            sup.setdefault(stop_of[r["model_id"]], {})["ribo"] = (int(r["psites"]), r["supported"] == "1")
    for r in csv.DictReader(open(f"{W}/proteomics/ms_peptide_support.{tag}.model_hits.tsv"), delimiter="\t"):
        if r["model_id"] in stop_of:
            sup.setdefault(stop_of[r["model_id"]], {})["ms"] = int(r["n_peptides"]) > 0
    return sup


def row(asm, name, keys, sup):
    rb = [sup[k]["ribo"][1] for k in keys if "ribo" in sup.get(k, {}) and sup[k]["ribo"][0] >= MIN_PSITES]
    ms = sum(1 for k in keys if sup.get(k, {}).get("ms"))
    n = len(keys)
    print(f"{asm}\t{name}\t{n}\t{len(rb)}\t{sum(rb)}\t{sum(rb) / len(rb) if rb else 0:.4f}\t{ms}\t{ms / n if n else 0:.4f}")


print("asm\tset\tn\tribo_n\tribo_sup\tribo_rate\tms_sup\tms_rate")
for asm in ("lraa", "stringtie", "trinity"):
    bs = stops(f"{G}/{asm}/m50_noq/transcripts.fa.transdecoder.gff3")
    ds = stops(f"{G}/{asm}/rel/transcripts.fa.transdecoder.gff3")
    bu, du = support(bs, f"{asm}_m50_noq"), support(ds, f"{asm}_rel")
    b, d = set(bs.values()), set(ds.values())
    row(asm, "baseline_no_q_test", b, bu)
    row(asm, "all", d, du)
    row(asm, "kept", d & b, bu)
    row(asm, "removed", b - d, bu)
    row(asm, "added", d - b, du)
