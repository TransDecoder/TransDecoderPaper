#!/usr/bin/env python3
"""Notebook 08: Ribo-seq and MS support of GM12878 ORFs that the release q-value test removes.

Release analogue of data/benchmark/rules/gm12878_rule_support.tsv (E20/E22, gm_rule_support.py), using
the release runs (E25) instead of the development decoy runs.

Baseline: the earlier final_abinitio runs (fixed short-ORF test, E03), with their Ribo-seq support
(ribo_gm12878/gm12878.<asm>.orf_support.tsv, E01) and MS peptide hits
(proteomics/ms_peptide_support.<asm>.model_hits.tsv, E02).
Release runs: <asm>/final_v7_abinitio[_precise|_sensitive], support from final_v7_scoring/ (E25).
ORFs are identified by transcript, strand and stop (3' end of the CDS in transcript coordinates).

Rows per assembly and setting:
  all      release predictions, with their own release support tables
  kept     baseline ORFs with the same stop in the release run (baseline support)
  dropped  baseline ORFs without a release prediction with that stop (baseline support)
  added    release ORFs whose stop is not in the baseline (release support)
plus one 'baseline all' row per assembly.
Columns: n ORFs (gene models); ribo_n = ORFs with >= 10 P-sites; ribo_sup = of those, passing the
frame test; ms_sup = ORFs with >= 1 MS peptide (rate over n).

usage: 08_rule_support_v7.py > data/gm12878/final_v7/08_rule_support_v7.tsv
"""
import csv
import os
import re
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
GM = os.path.join(P, "work/gm12878_td")
V7 = os.path.join(GM, "final_v7_scoring")
MIN_PSITES = 10
RUNS = {"lraa": ["default", "precise", "sensitive"], "stringtie": ["default", "precise", "sensitive"],
        "trinity": ["default"]}


def orf_stops(gff3):
    """ORF id -> (transcript, strand, 3' end) from CDS lines of a transcript-coordinate GFF3."""
    out = {}
    with open(gff3) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            orf = re.search(r"Parent=([^;]+)", f[8]).group(1)
            out[orf] = (f[0], f[6], int(f[4]) if f[6] == "+" else int(f[3]))
    return out


def support(stop_of, ribo_tsv, ms_tsv):
    """stop key -> {'ribo': (psites, supported), 'ms': bool}"""
    sup = {}
    with open(ribo_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["model_id"] in stop_of:
                sup.setdefault(stop_of[r["model_id"]], {})["ribo"] = (int(r["psites"]), r["supported"] == "1")
    with open(ms_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["model_id"] in stop_of:
                sup.setdefault(stop_of[r["model_id"]], {})["ms"] = int(r["n_peptides"]) > 0
    return sup


def row(asm, setting, name, keys, sup):
    rb = [sup[k]["ribo"][1] for k in keys if "ribo" in sup.get(k, {}) and sup[k]["ribo"][0] >= MIN_PSITES]
    ms = sum(1 for k in keys if sup.get(k, {}).get("ms"))
    n = len(keys)
    print(f"{asm}\t{setting}\t{name}\t{n}\t{len(rb)}\t{sum(rb)}\t{sum(rb) / len(rb) if rb else 0:.4f}\t"
          f"{ms}\t{ms / n if n else 0:.4f}")


def main():
    print("asm\tsetting\tset\tn\tribo_n\tribo_sup\tribo_rate\tms_sup\tms_rate")
    for asm, settings in RUNS.items():
        b_stop = orf_stops(os.path.join(GM, asm, "final_abinitio", "transcripts.fa.transdecoder.gff3"))
        b_sup = support(b_stop, os.path.join(GM, "ribo_gm12878", f"gm12878.{asm}.orf_support.tsv"),
                        os.path.join(GM, "proteomics", f"ms_peptide_support.{asm}.model_hits.tsv"))
        base = set(b_stop.values())
        row(asm, "baseline_fixed_short_orf_test", "all", base, b_sup)
        for s in settings:
            run = "final_v7_abinitio" + ("" if s == "default" else "_" + s)
            tag = asm + ("" if s == "default" else "_" + s)
            v_stop = orf_stops(os.path.join(GM, asm, run, "transcripts.fa.transdecoder.gff3"))
            v_sup = support(v_stop, os.path.join(V7, "ribo", f"gm12878_v7.{tag}.orf_support.tsv"),
                            os.path.join(V7, "proteomics", f"ms_peptide_support.{tag}.model_hits.tsv"))
            v = set(v_stop.values())
            row(asm, s, "all", v, v_sup)
            row(asm, s, "kept", v & base, b_sup)
            row(asm, s, "dropped", base - v, b_sup)
            row(asm, s, "added", v - base, v_sup)


if __name__ == "__main__":
    sys.exit(main())
