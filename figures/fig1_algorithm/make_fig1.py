#!/usr/bin/env python3
"""Figure 1: the TransDecoder algorithm.

usage: make_fig1.py [out_dir]     (default: figures/fig1_algorithm)

(a) workflow schematic (numbers and rules as in manuscript/sections/methods_algorithm.tex)
(b) GC-calibrated long-ORF length L*(GC), from analyses/long_orf_threshold/final_thresholds.tsv
    (q0.999), with the benchmark organisms at their transcriptome G+C (E07 Table 1)
(c) start-codon refinement: the worked 90-nt example of methods_algorithm.tex
Evidence records: E08, E09, E10, E04, E13.
"""
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

REPO = HERE.parents[1]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE
THRESH = REPO / "analyses/long_orf_threshold/final_thresholds.tsv"

INK, MUTED = "#222222", "#666666"
BOX_FC = {"stage": "#e8eef6", "model": "#eef6ea", "select": "#fbf1e3", "out": "#f2f2f2",
          "opt": "#ffffff"}
BOX_EC = {"stage": "#5b7fae", "model": "#5f9c4e", "select": "#c8892d", "out": "#888888",
          "opt": "#888888"}

# ------------------------------------------------------------------------------------------- (a)


def box(ax, x, y, w, h, title, body, kind, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
                                fc=BOX_FC[kind], ec=BOX_EC[kind], lw=0.8,
                                ls=(0, (3, 2)) if dashed else "-", mutation_aspect=1))
    ax.text(x + 0.008, y + h - 0.025, title, fontsize=6.5, fontweight="bold", va="top", ha="left",
            color=INK)
    ax.text(x + 0.008, y + h - 0.095, body, fontsize=5.8, va="top", ha="left", color=INK,
            linespacing=1.35)


def arrow(ax, x0, y0, x1, y1, dashed=False, conn="arc3"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#444444", shrinkA=0, shrinkB=0,
                                mutation_scale=7, ls=(0, (3, 2)) if dashed else "-",
                                connectionstyle=conn))


def panel_a(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    w, h, gap = 0.226, 0.40, 0.032
    x = [0.0 + i * (w + gap) for i in range(4)]
    y1, y2 = 0.56, 0.04

    # row 1: LongOrfs and training
    box(ax, x[0], y1, w, h, "Transcripts",
        "assembled transcripts (FASTA),\n"
        "or genome + transcript GTF\n"
        "-S: strand-specific data\n"
        "(forward strand only)", "out")
    box(ax, x[1], y1, w, h, "LongOrfs",
        "each stop paired with the most\n"
        "upstream in-frame ATG, both strands\n"
        "ORFs may run off either end:\n"
        "complete, 5' or 3' partial, internal\n"
        "keep ORFs of $\\geq$ 50 codons ($m$)", "stage")
    box(ax, x[2], y1, w, h, "Training set",
        "candidates: 10T longest ORFs\n"
        "(T = 500), CDS $\\leq$ 5 kb\n"
        "drop low-complexity and\n"
        "near-identical proteins\n"
        "keep the T longest", "model")
    box(ax, x[3], y1, w, h, "Self-trained models",
        "coding: 3-periodic, 5th-order\n"
        "Markov chain vs. background\n"
        "$\\mu$: median frame contrast per\n"
        "base of the training ORFs\n"
        "Kozak: start vs. internal ATGs", "model")

    # row 2: Predict
    box(ax, x[0], y2, w, h, "Score each ORF",
        "$S_1$: own frame\n"
        "$S_2, S_3$: shifted +1, +2 nt\n"
        "$S_4$-$S_6$: reverse complement\n"
        "frame contrast\n"
        "$C = S_1 - \\max(S_2, S_3)$", "select")
    box(ax, x[1], y2, w, h, "Keep if any rule holds",
        "1  homology hit (optional)\n"
        "2  CDS $\\geq L^\\ast$(GC)  (b)\n"
        "3  $S_1 > 0$ and $C > 0$;\n"
        "    if < 100 codons, also\n"
        "    $C \\geq 200\\,\\mu$ (short-ORF rule)", "select")
    box(ax, x[2], y2, w, h, "Rank",
        "by homology hits, then $S_1$,\n"
        "then length\n"
        "report the top ORF per\n"
        "transcript (or several non-\n"
        "overlapping ORFs, optional)", "select")
    box(ax, x[3], y2, w, h, "Refine starts (c), output",
        "open 5' end: trim to the first\n"
        "in-frame ATG if UTR > coding\n"
        "GFF3, BED, CDS and protein\n"
        "optional: project CDS onto\n"
        "the genome through the exons", "stage")

    for i in range(3):
        arrow(ax, x[i] + w, y1 + h / 2, x[i + 1], y1 + h / 2)
        arrow(ax, x[i] + w, y2 + h / 2, x[i + 1], y2 + h / 2)
    # row 1 -> row 2 (models feed scoring)
    xr, xl, ym = x[3] + w / 2, x[0] + w / 2, y1 - 0.035
    ax.plot([xr, xr, xl], [y1, ym, ym], color="#444444", lw=0.8)
    arrow(ax, xl, ym, xl, y2 + h)
    # optional homology search
    ax.text(x[1] + w / 2, y2 + h + 0.012, "optional: DIAMOND/BLASTP vs. proteins, hmmsearch vs. Pfam",
            fontsize=5.8, color=MUTED, ha="center", va="bottom", style="italic")
    ax.text(0.5, -0.04, "Perl TransDecoder v7.0.0 and PyTransDecoder v0.2.0 implement the same "
            "algorithm and give identical predictions", fontsize=6.3, color=MUTED, ha="center",
            va="top")


# ------------------------------------------------------------------------------------------- (b)


def panel_b(ax):
    rows = list(csv.DictReader(open(THRESH), delimiter="\t"))
    gc = [int(r["gc"]) for r in rows]
    lstar = [float(r["q0.999"]) for r in rows]
    # lookup: G+C rounded up to the next 5% step; below 25% the 25% value; above 80% no rule
    xs = [15] + gc
    ys = [lstar[0]] + lstar
    ax.step(xs, ys, where="pre", color=style.TD_C, lw=1.2, zorder=2)
    ax.axvspan(80, 90, color="#eeeeee", zorder=0)
    ax.text(85, 1800, "> 80% G+C:\nno length\nrule", fontsize=5.8, color=MUTED, ha="center",
            va="center")

    def lookup(g):
        step = max(25, 5 * math.ceil(g / 5))
        return dict(zip(gc, lstar))[step]

    label_off = {"pfalciparum": (3, -260, "center"), "spombe": (-2, 330, "right"),
                 "athal": (1, -330, "left"), "human": (-3, 420, "right"),
                 "creinhardtii": (-2, 330, "right")}
    for org in style.ORGS:
        g = style.ORG_GC[org]
        if org == "tthermophila":
            continue
        y = lookup(g)
        ax.plot(g, y, "o", ms=4, color=style.ORG_COLOR[org], mec="white", mew=0.5, zorder=3)
    labels = {"pfalciparum": "P. falciparum", "spombe": "S. pombe", "athal": "A. thaliana",
              "human": "human, mouse,\n$\\it{D.~melanogaster}$", "creinhardtii": "C. reinhardtii"}
    anchor = {"pfalciparum": "pfalciparum", "spombe": "spombe", "athal": "athal",
              "human": "dmel", "creinhardtii": "creinhardtii"}
    for key, text in labels.items():
        g = style.ORG_GC[anchor[key]]
        dx, dy, ha = label_off[key]
        ax.annotate(text, xy=(g, lookup(g)), xytext=(g + dx, lookup(g) + dy), fontsize=5.8,
                    ha=ha, va="bottom" if dy > 0 else "top", color=INK, style="italic" if "\n" not in text else "normal",
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED, shrinkA=0, shrinkB=2))
    ax.text(16, 2350, "$\\it{T.~thermophila}$ (26% G+C): rule off\n(ciliate code: TGA is the only stop)",
            fontsize=5.8, color=MUTED, ha="left", va="top")
    ax.set_xlim(15, 90)
    ax.set_ylim(0, 2600)
    ax.set_xticks(range(20, 91, 10))
    ax.set_xlabel("transcriptome G+C (%)")
    ax.set_ylabel("$L^\\ast$, CDS length kept\nwithout a coding test (nt)")


# ------------------------------------------------------------------------------------------- (c)


def panel_c(ax_t, ax_b):
    # transcript sketch
    ax = ax_t
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.axis("off")
    y = 4.2
    ax.plot([2, 6], [y, y], color="#888888", lw=1.2, ls=(0, (1, 1)))
    ax.add_patch(Rectangle((6, y - 1.6), 26, 3.2, fc="#fbe3c6", ec="#c8892d", lw=0.7))
    ax.add_patch(Rectangle((32, y - 1.6), 60, 3.2, fc="#c6d8ee", ec=style.TD_C, lw=0.7))
    ax.plot([92, 98], [y, y], color="#888888", lw=1.2)
    ax.text(19, y, "u  (90 nt)", ha="center", va="center", fontsize=6.3)
    ax.text(62, y, "ORF continues to the stop codon", ha="center", va="center", fontsize=6)
    ax.text(32, y + 2.2, "first in-frame ATG", ha="center", va="bottom", fontsize=6)
    ax.plot([32, 32], [y - 1.6, y + 2.0], color=INK, lw=0.8)
    ax.text(4, y + 2.2, "open 5' end", ha="center", va="bottom", fontsize=6, color=MUTED)
    ax.text(50, 0.2, "trim to the ATG if  UTR(u) = $\\alpha$|u| + $w$K  >  coding(u) = C(u)/$\\mu$",
            ha="center", va="bottom", fontsize=6.3)

    # scores for the worked example: |u| = 90, coding(u) = 45, alpha = 0.4, w = 17
    ax = ax_b
    alpha, w, u, coding = 0.4, 17, 90, 45
    base = alpha * u
    cases = [("K = +1", 1), ("K = 0", 0), ("K = -1", -1)]
    xs = range(len(cases))
    for i, (lab, k) in enumerate(cases):
        ax.bar(i, base, width=0.55, color="#fbe3c6", ec="#c8892d", lw=0.6)
        bonus = w * k
        if bonus > 0:
            ax.bar(i, bonus, bottom=base, width=0.55, color="#c8892d", lw=0)
        elif bonus < 0:
            ax.bar(i, -bonus, bottom=base + bonus, width=0.55, color="white", ec="#c8892d",
                   hatch="////", lw=0.6)
        total = base + bonus
        ax.text(i, max(total, base) + 2, f"{total:.0f}", ha="center", va="bottom", fontsize=6.3)
        ax.text(i, -9, "trim" if total > coding else "keep\npartial", ha="center", va="top",
                fontsize=6, color=style.TD_C if total > coding else MUTED, fontweight="bold")
    ax.axhline(coding, color=style.TD_C, lw=1.0, ls="--")
    ax.text(2.9, coding + 1.5, "coding(u) = 45", color=style.TD_C, fontsize=6, ha="right",
            va="bottom")
    ax.text(2.35, base / 2, "$\\alpha$|u| = 36", fontsize=5.8, ha="left", va="center",
            color=MUTED)
    ax.text(0, base + w / 2 - 2, "+$w$K", fontsize=5.8, ha="center", va="center", color="white")
    ax.text(2.33, base - w / 2, "$w$K < 0", fontsize=5.8, ha="left", va="center", color=MUTED)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([c[0] for c in cases])
    ax.tick_params(axis="x", length=0, pad=1)
    ax.set_xlim(-0.5, 2.9)
    ax.set_ylim(0, 62)
    ax.set_ylabel("UTR(u)\n(bases of coding signal)")
    ax.set_title("example: |u| = 90, coding(u) = 45;  $\\alpha$ = 0.4, $w$ = 17 bases per bit",
                 fontsize=5.8, color=MUTED, pad=3)


def main():
    fig = plt.figure(figsize=(style.FULL_W, 125 * style.MM))
    ax_a = fig.add_axes([0.02, 0.55, 0.96, 0.42])
    ax_b = fig.add_axes([0.085, 0.07, 0.36, 0.36])
    ax_ct = fig.add_axes([0.53, 0.33, 0.45, 0.13])
    ax_cb = fig.add_axes([0.60, 0.10, 0.33, 0.21])
    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_ct, ax_cb)
    fig.text(0.005, 0.975, "a", fontsize=9, fontweight="bold", va="top")
    fig.text(0.005, 0.47, "b", fontsize=9, fontweight="bold", va="top")
    fig.text(0.50, 0.47, "c", fontsize=9, fontweight="bold", va="top")
    style.save(fig, OUT, "fig1_algorithm")


if __name__ == "__main__":
    main()
