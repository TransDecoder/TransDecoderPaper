#!/usr/bin/env python3
"""Figure 7. Parameter dependence on real assemblies (GM12878).

usage: make_fig7.py [repo_root] [out_dir]
Panels:
  a  start refinement: alpha series at w = 17 per bit (= 25 per nat in the grid) and w = 0;
     x = MANE exact given assembled, y = exact-start share of GENCODE-compatible complete calls
     (analyses/start_param_sensitivity/results/gm12878.tsv; LRAA and StringTie only)
  b  post-hoc minimum length (--min-aa) on the final-default predictions: MANE exact given assembled
     and fraction of predictions with no GENCODE CDS overlap (genome_annotation/results/
     gm12878_minaa_sweep.tsv, headerless, columns in evidence/E04); pre-defined LongOrfs minimum
     runs overlaid (eval_summaries/{lraa,stringtie}_{kozak_add,m50,m50_sc200}.summary.tsv)
  c  homology: change in MANE exact given assembled over ab initio vs clade removed from Swiss-Prot,
     with and without Pfam; Pfam alone at right (gm12878_homology_series.tsv)
  d  homology: change in distinct predicted CDSs, same x; Trinity UTY/Alu artifact annotated (E06)
Evidence records: E04, E06, E13.
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

REPO = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[1]
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
RES = REPO / "genome_annotation/results"
ASM = ["lraa", "stringtie", "trinity"]
SHORT = {"lraa": "LRAA", "stringtie": "StringTie", "trinity": "Trinity"}
CLADES = [("9606", "human"), ("9443", "primates"), ("40674", "mammals"), ("7742", "vertebrates"),
          ("33208", "animals")]


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kv(path):
    with open(path) as fh:
        return {r[0]: r[1] for r in csv.reader(fh, delimiter="\t") if len(r) == 2}


def panel_a(ax):
    d = {(r["assembly"], r["config"]): r for r in
         rows(REPO / "analyses/start_param_sensitivity/results/gm12878.tsv")}
    alphas = ["0.2", "0.4", "0.6", "1.0"]
    for a in ["lraa", "stringtie"]:
        for w, ls, mfc in [("25", "-", None), ("0", "--", "white")]:
            xs = [float(d[(a, f"a{al}_w{w}")]["MANE_sens_exact_given_assembled"]) for al in alphas]
            ys = [float(d[(a, f"a{al}_w{w}")]["pred_complete_start_precision"]) for al in alphas]
            ax.plot(xs, ys, ls=ls, marker="o", ms=3, color=style.ASM_COLOR[a],
                    mfc=mfc or style.ASM_COLOR[a], lw=1)
            if w == "25":
                ax.plot(xs[1], ys[1], marker="o", ms=6.5, mfc="none", mec="k", mew=0.8)
                for k, (al, x, y) in enumerate(zip(alphas, xs, ys)):
                    if k in (0, 3):
                        ax.annotate(f"α = {al.rstrip('0').rstrip('.')}", (x, y),
                                    xytext=(-3, 3) if k == 0 else (3, -1),
                                    textcoords="offset points", ha="right" if k == 0 else "left",
                                    va="bottom" if k == 0 else "top", fontsize=5.5,
                                    color=style.ASM_COLOR[a])
    ax.set_xlim(0.83, 0.975)
    ax.set_xlabel("MANE exact given assembled")
    ax.set_ylabel("exact-start share of\ncompatible complete calls")
    ax.set_title("Start refinement (α)")
    ax.legend(handles=[Line2D([], [], color="0.3", marker="o", ms=3, label="w = 17/bit"),
                       Line2D([], [], color="0.3", ls="--", marker="o", ms=3, mfc="white",
                              label="w = 0"),
                       Line2D([], [], ls="none", marker="o", ms=6.5, mfc="none", mec="k", mew=0.8,
                              label="default")],
              loc="center left", bbox_to_anchor=(0.0, 0.5), handlelength=1.8, labelspacing=0.3)


def minaa_table():
    out = {}
    with open(RES / "gm12878_minaa_sweep.tsv") as fh:
        for r in csv.reader(fh, delimiter="\t"):
            # assembly, min aa, distinct CDSs, MANE sens exact, MANE exact given assembled,
            # GENCODE complete sens compatible, all-pred GENCODE exact frac, compatible frac,
            # no-GENCODE-overlap frac, complete-pred GENCODE exact frac (evidence/E04)
            out.setdefault(r[0], []).append((int(r[1]), float(r[4]), float(r[8])))
    return out


PREDEF = [("m50", "s", "−m 50, no rule"), ("m50_sc200", "*", "−m 50 + rule (default)"),
          ("kozak_add", "^", "−m 100")]


def panel_b(ax, col, ylabel, title):
    t = minaa_table()
    key = {1: "MANE_sens_exact_given_assembled", 2: "pred_all_no_GENCODE_CDS_overlap_frac"}[col]
    for a in ASM:
        xs = [r[0] for r in t[a]]
        ys = [r[col] for r in t[a]]
        ax.plot(xs, ys, "o-", ms=2.5, lw=1, color=style.ASM_COLOR[a])
        if a == "trinity":
            continue
        for cfg, m, _ in PREDEF:
            v = float(kv(RES / f"eval_summaries/{a}_{cfg}.summary.tsv")[key])
            x = 100 if cfg == "kozak_add" else 50
            ax.plot(x, v, marker=m, ms=5 if m == "*" else 3.8, mfc="white", mec=style.ASM_COLOR[a],
                    mew=0.9, ls="none", zorder=5)
    ax.axvline(50, color="0.5", lw=0.6, ls=":")
    ax.set_xticks([50, 100, 150, 200, 300])
    ax.set_xlabel("minimum protein length (aa)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def homology_rows():
    return {(r["assembly"], r["config"]): r for r in rows(RES / "gm12878_homology_series.tsv")}


def panel_cd(ax, field, scale, ylabel, title, annotate=False):
    d = homology_rows()
    xs = list(range(len(CLADES)))
    for a in ASM:
        base = float(d[(a, "abinitio")][field])
        for pf, ls, mfc in [("nopfam", "-", None), ("pfam", "--", "white")]:
            ys = [(float(d[(a, f"blast_{t}_{pf}")][field]) - base) * scale for t, _ in CLADES]
            ax.plot(xs, ys, ls=ls, marker="o", ms=3, lw=1, color=style.ASM_COLOR[a],
                    mfc=mfc or style.ASM_COLOR[a])
        ax.plot(len(CLADES) + 0.4, (float(d[(a, "pfam")][field]) - base) * scale, marker="D", ms=3.5,
                color=style.ASM_COLOR[a], ls="none")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xticks(xs + [len(CLADES) + 0.4])
    ax.set_xticklabels([f"−{c}" for _, c in CLADES] + ["Pfam\nonly"], rotation=30, ha="right",
                       rotation_mode="anchor")
    ax.set_xlim(-0.4, len(CLADES) + 0.8)
    ax.set_xlabel("Swiss-Prot with clade removed")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if annotate:
        # 4,465 Trinity CDSs found with -human but not -primates have chimp UTY as best hit (E06)
        h = (float(d[("trinity", "blast_9606_nopfam")][field]) - float(d[("trinity", "abinitio")][field]))
        ax.annotate("4,465 short unannotated ORFs kept only\nby homology to the Alu-encoded\n"
                    "C-terminus of chimpanzee UTY",
                    xy=(0.08, h * scale), xytext=(1.0, h * scale * 0.98), fontsize=5.5, va="top",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3"))


def main():
    fig = plt.figure(figsize=(style.FULL_W, 4.7))
    gs = GridSpec(2, 6, figure=fig, hspace=0.68, wspace=2.0)
    axa = fig.add_subplot(gs[0, 0:2])
    axb1 = fig.add_subplot(gs[0, 2:4])
    axb2 = fig.add_subplot(gs[0, 4:6])
    axc = fig.add_subplot(gs[1, 0:3])
    axd = fig.add_subplot(gs[1, 3:6])
    panel_a(axa)
    panel_b(axb1, 1, "MANE exact given\nassembled", "Minimum length: accuracy")
    panel_b(axb2, 2, "no GENCODE CDS\noverlap (fraction)", "Minimum length: unannotated")
    axb2.legend(handles=[Line2D([], [], color="0.3", marker="o", ms=2.5, label="post-hoc filter")]
                + [Line2D([], [], ls="none", marker=m, ms=5 if m == "*" else 3.8, mfc="white",
                          mec="0.3", mew=0.9, label=lab) for _, m, lab in PREDEF],
                loc="upper right", bbox_to_anchor=(1.0, 0.86), handlelength=1.5, labelspacing=0.3)
    panel_cd(axc, "MANE_sens_exact_given_assembled", 100,
             "Δ MANE exact given assembled\n(percentage points)", "Homology evidence: accuracy")
    panel_cd(axd, "distinct_predicted_cds", 1e-3, "Δ distinct predicted CDSs\n(thousands)",
             "Homology evidence: predictions", annotate=True)
    axc.legend(handles=[Line2D([], [], color="0.3", marker="o", ms=3, label="Swiss-Prot"),
                        Line2D([], [], color="0.3", ls="--", marker="o", ms=3, mfc="white",
                               label="Swiss-Prot + Pfam"),
                        Line2D([], [], color="0.3", ls="none", marker="D", ms=3.5, label="Pfam only")],
               loc="lower left", handlelength=1.8, labelspacing=0.3)
    for ax, l in zip([axa, axb1, axb2, axc, axd], ["a", "b", "", "c", "d"]):
        if l:
            style.panel_label(ax, l, x=-0.36 if ax in (axa, axb1) else -0.17, y=1.07)
    fig.legend(handles=[Patch(color=style.ASM_COLOR[a], label=style.ASM_LABEL[a]) for a in ASM],
               loc="lower center", bbox_to_anchor=(0.5, 0.93), ncol=3, handlelength=1)
    style.save(fig, OUT, "fig7_parameters")


if __name__ == "__main__":
    main()
