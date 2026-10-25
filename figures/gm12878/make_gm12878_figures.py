#!/usr/bin/env python3
"""GM12878 application figures, drawn from committed result tables only.

usage: make_gm12878_figures.py [results_dir] [out_dir]
  results_dir default: data/gm12878 (relative to the repo root)
  out_dir     default: figures/gm12878

Figures (PDF + PNG):
  gm12878_min_length   minimum ORF length (post-hoc --min-aa) vs sensitivity and precision measures
  gm12878_evidence     Ribo-seq footprint coverage and MS peptide support by prediction class;
                       footprint coverage by RNA expression (StringTie)
  gm12878_homology     homology restriction series (clade-excluded Swiss-Prot, Pfam)
Evidence records: E01, E02, E04, E06.
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
RES = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "data/gm12878"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else REPO / "figures/gm12878"
ASM = ["lraa", "stringtie", "trinity"]
ASM_LABEL = {"lraa": "LRAA (long reads)", "stringtie": "StringTie", "trinity": "Trinity (genome-guided)"}
COLOR = {"lraa": "#1b9e77", "stringtie": "#d95f02", "trinity": "#7570b3"}
CLASSES = ["exact", "compatible", "other_overlap", "novel"]
CLASS_LABEL = {"exact": "GENCODE\nexact", "compatible": "GENCODE\ncompatible",
               "other_overlap": "other CDS\noverlap", "novel": "unannotated"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def kv(path):
    return {r[0]: r[1] for r in csv.reader(open(path), delimiter="\t") if len(r) == 2}


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def min_length():
    mins = [50, 60, 70, 80, 90, 100, 150, 200, 300]
    panels = [("MANE_sens_exact_given_assembled", "MANE exact | assembled"),
              ("MANE_sens_compatible", "MANE compatible (sensitivity)"),
              ("pred_complete_GENCODE_exact_frac", "complete predictions:\nGENCODE exact fraction"),
              ("pred_complete_no_GENCODE_CDS_overlap_frac", "complete predictions:\nno GENCODE CDS overlap")]
    fig, axs = plt.subplots(1, 4, figsize=(9, 2.3))
    for ax, (key, title) in zip(axs, panels):
        for a in ASM:
            y = [float(kv(RES / f"eval_summaries/minaa/{a}_min{m}.summary.tsv")[key]) for m in mins]
            ax.plot(mins, y, "o-", ms=3, color=COLOR[a], label=ASM_LABEL[a])
        ax.axvline(100, color="grey", lw=0.6, ls=":")
        ax.set_title(title)
        ax.set_xlabel("minimum ORF length (aa)")
    axs[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "gm12878_min_length")


def evidence():
    ribo = {(r["assembly"], r["class"]): r for r in rows(RES / "gm12878_riboseq_support.tsv")}
    ms = {(r["assembly"], r["class"]): r for r in rows(RES / "gm12878_ms_peptide_support.tsv") if r["size"] == "all"}
    expr = rows(RES / "gm12878_riboseq_stringtie_by_expression.tsv")
    fig, axs = plt.subplots(1, 3, figsize=(10, 2.6), gridspec_kw={"width_ratios": [1.3, 1.3, 1]})
    w = 0.26
    for ax, title, val in (
            (axs[0], "Ribo-seq: ≥10 P-sites", lambda a, c: int(ribo[(a, c)]["with_min_psites"]) / int(ribo[(a, c)]["n"])),
            (axs[1], "MS: ≥1 peptide", lambda a, c: float(ms[(a, c)]["frac_ge1"]))):
        for i, a in enumerate(ASM):
            ax.bar([j + (i - 1) * w for j in range(len(CLASSES))], [val(a, c) for c in CLASSES], w,
                   color=COLOR[a], label=ASM_LABEL[a])
        ax.set_xticks(range(len(CLASSES)), [CLASS_LABEL[c] for c in CLASSES])
        ax.set_ylim(0, 1)
        ax.set_ylabel("fraction of predicted CDSs")
        ax.set_title(title)
    axs[0].legend(frameon=False, fontsize=7, loc="upper right")
    bins = ["a<1", "b1-5", "c5-20", "d20-100", "e>=100"]
    ax = axs[2]
    for c, col in (("exact", "#333333"), ("novel", "#e7298a")):
        d = {r["tpm_bin"]: float(r["frac_ge10_psites"]) for r in expr if r["class"] == c}
        ax.plot(range(len(bins)), [d[b] for b in bins], "o-", color=col, ms=3, label=CLASS_LABEL[c].replace("\n", " "))
    ax.set_xticks(range(len(bins)), ["<1", "1–5", "5–20", "20–100", "≥100"])
    ax.set_xlabel("source transcript TPM (StringTie)")
    ax.set_ylabel("fraction with ≥10 P-sites")
    ax.set_ylim(0, 1)
    ax.set_title("Ribo-seq coverage vs expression")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save(fig, "gm12878_evidence")


def homology():
    t = {(r["assembly"], r["config"]): r for r in rows(RES / "gm12878_homology_series.tsv")}
    clades = [("9606", "human"), ("9443", "primates"), ("40674", "mammals"), ("7742", "vertebrates"),
              ("33208", "Metazoa")]
    panels = [("MANE_sens_compatible", "MANE compatible"),
              ("MANE_sens_exact_given_assembled", "MANE exact | assembled"),
              ("distinct_predicted_cds", "distinct predicted CDSs (×1000)")]
    fig, axs = plt.subplots(1, 3, figsize=(9, 2.5))
    x = range(len(clades))
    for ax, (key, title) in zip(axs, panels):
        f = (lambda v: float(v) / 1000) if key == "distinct_predicted_cds" else float
        for a in ASM:
            ax.plot(x, [f(t[(a, f"blast_{c}_nopfam")][key]) for c, _ in clades], "o-", ms=3, color=COLOR[a],
                    label=ASM_LABEL[a])
            ax.plot(x, [f(t[(a, f"blast_{c}_pfam")][key]) for c, _ in clades], "s--", ms=2.5, lw=0.8,
                    color=COLOR[a])
            ax.axhline(f(t[(a, "abinitio")][key]), color=COLOR[a], lw=0.6, ls=":")
        ax.set_xticks(x, [n for _, n in clades], rotation=30)
        ax.set_xlabel("Swiss-Prot with clade removed")
        ax.set_title(title)
    axs[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.text(0.5, -0.02, "solid: Swiss-Prot only;  dashed: Swiss-Prot + Pfam;  dotted: ab initio", ha="center",
             fontsize=7)
    save(fig, "gm12878_homology")


if __name__ == "__main__":
    min_length()
    evidence()
    homology()
    print(f"figures written to {OUT}")
