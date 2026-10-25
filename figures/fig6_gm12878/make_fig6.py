#!/usr/bin/env python3
"""Figure 6. PyTransDecoder on three GM12878 assemblies, tested against Ribo-seq and MS.

usage: make_fig6.py [results_dir] [out_dir]
  results_dir default: data/gm12878 (relative to the repo root)
  out_dir     default: figures/fig6_gm12878
Panels:
  a  MANE Select CDSs assembled and predicted exactly, all genes and expressed genes
     (gm12878_expressed_sensitivity.tsv); exact given assembled printed above each bar
  b  class composition of distinct predicted CDSs by predicted length (data/pred_class_by_length.tsv,
     from derive_class_by_length.py)
  c  Ribo-seq: >=10 P-sites and frame-0 test by class (gm12878_riboseq_support.tsv)
  d  MS: >=1 peptide by class, >=100 aa (bars) and <100 aa (markers) (gm12878_ms_peptide_support.tsv)
  e  Ribo-seq coverage by source-transcript TPM, StringTie (gm12878_riboseq_stringtie_by_expression.tsv)
Evidence records: E01, E02, E03, E04, E05.
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPO = HERE.parents[1]
RES = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "data/gm12878"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
ASM = ["lraa", "stringtie", "trinity"]
SHORT = {"lraa": "LRAA", "stringtie": "StringTie", "trinity": "Trinity"}
CLASSES = ["exact", "compatible", "other_overlap", "novel"]
CLASS_LABEL = {"exact": "GENCODE\nexact", "compatible": "GENCODE\ncompatible",
               "other_overlap": "other CDS\noverlap", "novel": "unannotated"}
SHORT_CLASS = {"exact": "exact", "compatible": "compatible", "other_overlap": "other\noverlap",
               "novel": "unanno-\ntated"}
SHORT_CLASS["compatible"] = "compat-\nible"
CLASS_COLOR = {"exact": "#08519c", "compatible": "#6baed6", "other_overlap": "#bdbdbd",
               "novel": "#e7298a"}


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def panel_a(ax):
    d = {(r["label"], r["tpm_min"]): r for r in rows(RES / "gm12878_expressed_sensitivity.tsv")
         if r["set"] == "MANE"}
    tpms = [("0", "all genes"), ("1", "TPM ≥ 1"), ("10", "TPM ≥ 10")]
    w = 0.26
    for i, a in enumerate(ASM):
        for j, (t, _) in enumerate(tpms):
            r = d[(a, t)]
            x = i + (j - 1) * (w + 0.02)
            asm_f, ex = float(r["assembled_frac"]), float(r["sens_exact"])
            ax.bar(x, asm_f, w, color="white", edgecolor=style.ASM_COLOR[a], lw=0.8)
            ax.bar(x, ex, w, color=style.ASM_COLOR[a], alpha=[1, 0.7, 0.45][j], lw=0)
            ax.text(x, asm_f + 0.015, f"{float(r['sens_exact_given_assembled']):.2f}",
                    ha="center", va="bottom", fontsize=5.5, rotation=90)
            ax.text(x, -0.03, ["all", "≥1", "≥10"][j], ha="center", va="top", fontsize=4.8)
    ax.set_xticks(range(3))
    ax.set_xticklabels([SHORT[a] for a in ASM])
    ax.tick_params(axis="x", length=0, pad=9)
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.spines["left"].set_bounds(0, 1)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylabel("fraction of MANE Select CDSs")
    ax.set_title("MANE Select recovery")
    ax.legend(handles=[Patch(facecolor="white", edgecolor="k", lw=0.8, label="assembled"),
                       Patch(facecolor="0.3", label="exact")],
              loc="upper center", ncol=2, handlelength=1, columnspacing=0.8)


def panel_b(ax):
    d = rows(HERE / "data/pred_class_by_length.tsv")
    bins = ["lt100", "100-149", "150-199", "200-299", "ge300"]
    blabel = ["<100", "100–\n149", "150–\n199", "200–\n299", "≥300"]
    gap = 0.8
    xt, xl = [], []
    for i, a in enumerate(ASM):
        for j, b in enumerate(bins):
            r = next(r for r in d if r["assembly"] == a and r["length_bin"] == b)
            n = int(r["n"])
            x = i * (len(bins) + gap) + j
            bottom = 0
            for c in CLASSES:
                f = int(r[c]) / n
                ax.bar(x, f, 0.85, bottom=bottom, color=CLASS_COLOR[c], lw=0)
                bottom += f
            xt.append(x)
            xl.append(blabel[j])
        ax.text(i * (len(bins) + gap) + 2, -0.17, SHORT[a], ha="center", va="top",
                color=style.ASM_COLOR[a], fontweight="bold", transform=ax.get_xaxis_transform())
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=5)
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.7, xt[-1] + 0.7)
    ax.set_ylabel("fraction of predicted CDSs")
    ax.set_xlabel("predicted protein length (aa)", labelpad=14)
    lab = {"exact": "GENCODE exact", "compatible": "compatible", "other_overlap": "other overlap",
           "novel": "unannotated"}
    ax.legend(handles=[Patch(color=CLASS_COLOR[c], label=lab[c]) for c in CLASSES], loc="lower right", bbox_to_anchor=(1.0, 1.0),
              ncol=4, handlelength=1, columnspacing=1.0)


def class_bars(ax, get, ylabel, title):
    w = 0.27
    for i, a in enumerate(ASM):
        for j, c in enumerate(CLASSES):
            get(ax, a, c, j + (i - 1) * w, w)
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels([SHORT_CLASS[c] for c in CLASSES], fontsize=5.5)
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.spines["left"].set_bounds(0, 1)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def panel_c(ax):
    d = {(r["assembly"], r["class"]): r for r in rows(RES / "gm12878_riboseq_support.tsv")}

    def get(ax, a, c, x, w):
        r = d[(a, c)]
        n = int(r["n"])
        ax.bar(x, int(r["with_min_psites"]) / n, w, color=style.ASM_COLOR[a], alpha=0.45, lw=0)
        ax.bar(x, int(r["supported"]) / n, w, color=style.ASM_COLOR[a], lw=0)
        if c == "novel":
            ax.text(x, int(r["with_min_psites"]) / n + 0.02, f"n={n:,}", rotation=90,
                    ha="center", va="bottom", fontsize=5)

    class_bars(ax, get, "fraction of predicted CDSs", "Ribo-seq footprints")
    ax.legend(handles=[Patch(color="0.6", label="≥10 P-sites"),
                       Patch(color="0.2", label="frame 0")],
              loc="upper center", ncol=2, handlelength=1, columnspacing=1)


def panel_d(ax):
    d = {(r["assembly"], r["class"], r["size"]): r for r in rows(RES / "gm12878_ms_peptide_support.tsv")}

    def get(ax, a, c, x, w):
        ax.bar(x, float(d[(a, c, "long")]["frac_ge1"]), w, color=style.ASM_COLOR[a], lw=0)
        ax.plot(x, float(d[(a, c, "short")]["frac_ge1"]), marker="D", ms=2.6, color="k",
                mfc="white", mew=0.6, ls="none")

    class_bars(ax, get, "fraction with ≥ 1 peptide", "MS peptides")
    ax.legend(handles=[Patch(color="0.4", label="≥ 100 aa"),
                       Line2D([], [], marker="D", ms=2.6, color="k", mfc="white", mew=0.6, ls="none",
                              label="< 100 aa")],
              loc="upper center", ncol=2, handlelength=1, columnspacing=1)


def panel_e(ax):
    d = rows(RES / "gm12878_riboseq_stringtie_by_expression.tsv")
    tb = ["a<1", "b1-5", "c5-20", "d20-100", "e>=100"]
    tl = ["<1", "1–5", "5–20", "20–100", "≥100"]
    for c in CLASSES:
        y = [float(next(r for r in d if r["class"] == c and r["tpm_bin"] == t)["frac_ge10_psites"])
             for t in tb]
        ax.plot(range(len(tb)), y, "o-", color=CLASS_COLOR[c] if c != "other_overlap" else "0.55",
                label=SHORT_CLASS[c].replace("-\n", "").replace("\n", " "), ms=2.5, lw=1)
    ax.set_xticks(range(len(tb)))
    ax.set_xticklabels(tl)
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.spines["left"].set_bounds(0, 1)
    ax.set_xlabel("source transcript TPM")
    ax.set_ylabel("fraction with ≥ 10 P-sites")
    ax.set_title("Ribo-seq by expression (StringTie)")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.12), handlelength=1.2, labelspacing=0.3)


def main():
    fig = plt.figure(figsize=(style.FULL_W, 5.0))
    gs = GridSpec(2, 6, figure=fig, height_ratios=[1, 1], hspace=0.7, wspace=1.7)
    axa = fig.add_subplot(gs[0, 0:2])
    axb = fig.add_subplot(gs[0, 2:6])
    axc = fig.add_subplot(gs[1, 0:2])
    axd = fig.add_subplot(gs[1, 2:4])
    axe = fig.add_subplot(gs[1, 4:6])
    panel_a(axa)
    panel_b(axb)
    panel_c(axc)
    panel_d(axd)
    panel_e(axe)
    for ax, l in zip([axa, axb, axc, axd, axe], "abcde"):
        style.panel_label(ax, l, x=-0.30 if ax is not axb else -0.13, y=1.04)
    pos = axc.get_position()
    fig.legend(handles=[Patch(color=style.ASM_COLOR[a], label=style.ASM_LABEL[a]) for a in ASM],
               loc="lower left", bbox_to_anchor=(pos.x0, pos.y1 + 0.045), ncol=3, handlelength=1)
    style.save(fig, OUT, "fig6_gm12878")


if __name__ == "__main__":
    main()
