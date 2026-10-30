#!/usr/bin/env python3
"""LaTeX body of the supplementary table of the GM12878 homology series (E06 release rerun).

Reads data/gm12878/gm12878_homology_series_v7h.tsv and prints the table rows (three assembly blocks,
twelve configurations each) used in manuscript/supplement.tex (tab:hseries).

usage (from the paper repo root):
  python3 pipeline/derive/08_homology_series_v7h_table.py > data/gm12878/gm12878_homology_series_v7h.tex
"""
import csv

ASM = [("lraa", "LRAA (long reads)"), ("stringtie", "StringTie (short reads)"),
       ("trinity", "Trinity genome-guided (short reads)")]
CLADE = {"9606": "human", "9443": "primates", "40674": "mammals", "7742": "vertebrates", "33208": "animals"}
ORDER = ["abinitio", "pfam"] + [f"blast_{c}_nopfam" for c in CLADE] + [f"blast_{c}_pfam" for c in CLADE]


def label(cfg):
    if cfg == "abinitio":
        return "Ab initio"
    if cfg == "pfam":
        return "Pfam only"
    _, taxid, pf = cfg.split("_")
    return f"Swiss-Prot $-$ {CLADE[taxid]}" + (" + Pfam" if pf == "pfam" else "")


rows = list(csv.DictReader(open("data/gm12878/gm12878_homology_series_v7h.tsv"), delimiter="\t"))
t = {(r["assembly"], r["config"]): r for r in rows}
out = []
for a, name in ASM:
    out.append(r"\multicolumn{7}{l}{\textit{" + name + r"}} \\")
    base = float(t[(a, "abinitio")]["MANE_sens_exact_given_assembled"])
    for cfg in ORDER:
        r = t[(a, cfg)]
        ega = float(r["MANE_sens_exact_given_assembled"])
        d = "" if cfg == "abinitio" else f"{100 * (ega - base):+.1f}"
        out.append(" & ".join([
            r"\quad " + label(cfg),
            f"{int(r['distinct_predicted_cds']):,}".replace(",", "{,}"),
            f"{float(r['MANE_sens_exact']):.3f}",
            f"{ega:.3f}",
            d.replace("-", "$-$"),
            f"{float(r['MANE_sens_compatible']):.3f}",
            f"{100 * float(r['pred_complete_no_GENCODE_CDS_overlap_frac']):.1f}",
        ]) + r" \\")
    if a != "trinity":
        out.append(r"\addlinespace")
print("\n".join(out))
