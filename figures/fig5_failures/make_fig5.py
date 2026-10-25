#!/usr/bin/env python3
"""Figure 5: true CDSs that fail the coding test, and the two fallbacks (length, homology).

usage: make_fig5.py [homology_vs_abinitio.tsv] [out_dir]
  homology_vs_abinitio.tsv default: data/benchmark/homology_vs_abinitio.tsv. Regenerate it
  from a new summary.tsv with benchmark_v2/scripts/compare_homology.py and pass its path here.
  out_dir default: figures/fig5_failures

Panels
  a  failure rate of matched true CDSs by CDS length, per organism; pooled failure and pooled
     "still lost after the length rule"            data/long_orf_failures/summary_by_length.tsv
  b  exclusive attribution of failures >= 300 nt   data/long_orf_failures/explanation_flags.tsv
  c  D. melanogaster failures >= 300 nt by GC3 quintile and chromosome, with the share rescued by
     a GC3-matched model                           data/dmel_gc3_rescue.tsv (derive_dmel_gc3_rescue.py)
  d  homology minus ab initio, stop-level sensitivity, stranded complete vs fragments
  e  same, exact-CDS F1                            homology_vs_abinitio.tsv
  f  failing CDSs retained by the length rule, with or without a homology hit
                                                   data/long_orf_failures/length_vs_homology_fallback.tsv
Evidence: E12 (a-c, f), E07 (d, e).
"""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402

REPO = HERE.parents[1]
LOF = REPO / "data/long_orf_failures"
HOM = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "data/benchmark/homology_vs_abinitio.tsv"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE

BINS = ["0-299", "300-599", "600-999", "1000-1999", "2000-4999", ">=5000"]
BIN_LABEL = ["<0.3", "0.3–\n0.6", "0.6–\n1", "1–2", "2–5", "≥5"]
ATTR = [("atypical_gc3_composition", "atypical GC3", "#2c7fb8"),
        ("low_complexity_or_repetitive", "low complexity /\nrepetitive", "#d95f02"),
        ("no_homology_orphan", "no homology", "#7570b3"),
        ("stronger_out_of_frame_orf", "stronger ORF in\nanother frame", "#e6ab02"),
        ("unexplained", "unexplained", "#bdbdbd")]
SHORT = {o: style.ORG_LABEL[o].replace(". ", ".") for o in style.ORGS}
SHORT.update({"pfalciparum": "P. falc.", "tthermophila": "T. therm.", "spombe": "S. pombe",
              "athal": "A. thal.", "dmel": "D. mel.", "mouse": "M. mus.", "human": "H. sap.",
              "creinhardtii": "C. rein."})


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def italic(o):
    return r"$\it{" + SHORT[o].replace(" ", r"\ ") + "}$"


def panel_a(ax):
    r = rows(LOF / "summary_by_length.tsv")
    x = np.arange(len(BINS))
    for o in style.ORGS:
        y = [float(next(t for t in r if t["org"] == o and t["cds_len_bin"] == b)["fail_coding_pct"])
             for b in BINS]
        y = [v if v > 0 else np.nan for v in y]
        ax.plot(x, y, "o-", color=style.ORG_COLOR[o], lw=0.8, ms=2.5, alpha=0.85, label=style.ORG_LABEL[o])
    for col, ls, lab in [("fail_coding_pct", "-", "all, failing"), ("still_lost_pct", "--", "all, still lost\nafter length rule")]:
        y = [float(next(t for t in r if t["org"] == "ALL" and t["cds_len_bin"] == b)[col]) for b in BINS]
        y = [v if v > 0 else np.nan for v in y]
        ax.plot(x, y, ls, color="black", lw=1.8, label=lab, zorder=5)
    ax.set_yscale("log")
    ax.set_ylim(0.015, 150)
    ax.set_yticks([0.1, 1, 10, 100])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.set_xticks(x)
    ax.set_xticklabels(BIN_LABEL)
    ax.set_xlabel("annotated CDS length (kb)")
    ax.set_ylabel("true CDSs failing the coding test (%)")
    h, l = ax.get_legend_handles_labels()
    leg = ax.legend(h[-2:], l[-2:], loc="upper right", fontsize=6, handlelength=2.2)
    ax.add_artist(leg)
    ax.legend(h[:-2], [italic(o) for o in style.ORGS], loc="upper left", bbox_to_anchor=(-0.05, -0.3),
              ncol=2, fontsize=6, handlelength=1.4, columnspacing=0.8, labelspacing=0.3)


def panel_b(ax):
    r = rows(LOF / "explanation_flags.tsv")
    orgs = style.ORGS + ["ALL"]
    y = np.arange(len(orgs))[::-1]
    for yi, o in zip(y, orgs):
        left = 0
        for key, lab, col in ATTR:
            t = next(t for t in r if t["org"] == o and t["flag"] == "attr:" + key)
            v = float(t["fail_pct"])
            ax.barh(yi, v, left=left, color=col, height=0.72, edgecolor="white", linewidth=0.6,
                    label=lab if o == "ALL" else None)
            left += v
        n = next(t for t in r if t["org"] == o and t["flag"].startswith("attr:"))["n_fail"]
        ax.text(101, yi, f"{int(n):,}", va="center", fontsize=6, color="#444444")
    ax.set_yticks(y)
    ax.set_yticklabels([italic(o) if o != "ALL" else "all" for o in orgs])
    ax.set_xlim(0, 100)
    ax.set_xlabel("failures ≥300 nt (%)")
    ax.text(101, y[0] + 0.85, "n", fontsize=6, color="#444444")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.05, -0.2), ncol=2, fontsize=6, handlelength=1,
              columnspacing=0.8)


def panel_c(ax):
    r = {t["group"]: t for t in rows(HERE / "data/dmel_gc3_rescue.tsv")}
    groups = [f"gc3_q{i}" for i in range(1, 6)] + ["chr4", "chrY"]
    labels = ["Q1\n(low)", "Q2", "Q3", "Q4", "Q5\n(high)", "4", "Y"]
    x = np.array([0, 1, 2, 3, 4, 5.6, 6.6])
    fail = np.array([100 * int(r[g]["fail_score"]) / int(r[g]["n"]) for g in groups])
    resc = np.array([100 * int(r[g]["rescued_gc3_model"]) / int(r[g]["n"]) for g in groups])
    ax.bar(x, fail, width=0.72, color="#cfcfcf", label="fail")
    ax.bar(x, resc, width=0.72, color="#2c7fb8", label="rescued by\nGC3-matched model")
    for xi, g, f in zip(x, groups, fail):
        if g in ("gc3_q1", "chr4", "chrY"):
            ax.text(xi, f + 1, f"{r[g]['fail_score']}/{r[g]['n']}", ha="center", va="bottom", fontsize=5.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.text(1.8, -0.2, "GC3 quintile", transform=ax.get_xaxis_transform(), ha="center", va="top")
    ax.text(6.1, -0.2, "chromosome", transform=ax.get_xaxis_transform(), ha="center", va="top")
    ax.set_ylabel("CDSs ≥300 nt failing (%)")
    ax.set_title(r"$\it{D.\ melanogaster}$")
    ax.set_ylim(0, 58)
    ax.axvline(4.8, color="#bbbbbb", lw=0.5, ls=":")
    ax.legend(loc="center", bbox_to_anchor=(0.38, 0.62), fontsize=6, handlelength=1)


def homology(ax, metric, ylabel):
    r = rows(HOM)
    x = np.arange(len(style.ORGS))
    for v, col, lab, dx in [("complete", style.COMPLETE_C, "complete", -0.12),
                            ("fragments", style.FRAG_C, "fragments", 0.12)]:
        y = [100 * float(next(t for t in r if t["org"] == o and t["variant"] == v)[metric]) for o in style.ORGS]
        ax.bar(x + dx, y, width=0.24, color=col, label=lab)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([italic(o) for o in style.ORGS], rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel(ylabel)
    ax.set_ylim(-0.5, 10)


def panel_f(ax):
    r = {t["org"]: t for t in rows(LOF / "length_vs_homology_fallback.tsv")}
    orgs = [o for o in style.ORGS if o != "tthermophila"]
    y = np.arange(len(orgs))[::-1]
    hom = np.array([int(r[o]["also_homology_hit"]) for o in orgs])
    only = np.array([int(r[o]["length_only"]) for o in orgs])
    ax.barh(y, hom, color="#9ecae1", height=0.72, label="also a homology hit")
    ax.barh(y, only, left=hom, color="#08519c", height=0.72, label="length rule only")
    for yi, h, o in zip(y, hom, only):
        ax.text(h + o + 3, yi, f"{o}/{h + o}", va="center", fontsize=6, color="#444444")
    ax.set_yticks(y)
    ax.set_yticklabels([italic(o) for o in orgs])
    ax.set_xlim(0, 215)
    ax.set_xlabel("failing true CDSs kept by\nthe length rule")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    a = r["ALL"]
    ax.set_title(f"all: {a['length_only']} of {a['rescued_by_length']} length rule only", fontsize=6.5)
    ax.legend(loc="upper right", fontsize=6, handlelength=1)


def main():
    fig = plt.figure(figsize=(style.FULL_W, 140 * style.MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.0, 1.0], height_ratios=[1, 0.9],
                          wspace=0.55, hspace=0.9)
    axs = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]
    panel_a(axs[0])
    panel_b(axs[1])
    panel_c(axs[2])
    homology(axs[3], "d_sensitivity_stop", "homology − ab initio:\nstop-level sensitivity (points)")
    axs[3].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6, handlelength=1)
    homology(axs[4], "d_f1_exact", "homology − ab initio:\nexact-CDS F1 (points)")
    panel_f(axs[5])
    for ax, letter in zip(axs, "abcdef"):
        style.panel_label(ax, letter, x=-0.3 if letter in "bf" else -0.22)
    style.save(fig, OUT, "fig5_failures")


if __name__ == "__main__":
    main()
