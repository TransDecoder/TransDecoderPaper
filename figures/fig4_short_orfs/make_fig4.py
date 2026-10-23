#!/usr/bin/env python3
"""Figure 4: minimum ORF length and the short-ORF rule.

usage: make_fig4.py [results_dir] [out_dir]
  results_dir default: benchmark_v2/results_v3
  out_dir     default: figures/fig4_short_orfs

Panels
  a  human and mouse complete transcripts + lncRNA negatives: share of lncRNAs with a prediction vs
     stop-level sensitivity for true proteins of 50-99 aa; LongOrfs minimum 100..50 codons, and
     50 codons + the short-ORF rule at factors 150/200/250             short_orf_summary.tsv
  b  same, fragments + lncRNA negatives                                 short_orf_summary.tsv
  c  predictions per negative sequence at 100 codons, 50 codons and 50 codons + rule (factor 200):
     GC-matched random sequences (complete.wRand, 8 organisms) and lncRNAs (complete.wLnc)
                                                                        short_orf_summary.tsv
  d  share of human lncRNAs with a prediction by method, lncRNAs already in GENCODE v28 vs added
     later (complete.wLnc)                                              lnc_by_era.human.tsv
Evidence: E04 (a-c), E07 (d).
Negatives per input: one random sequence per truth transcript (Table 1 gene counts); lncRNAs 33,668
(human) and 31,467 (mouse).
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
RES = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "benchmark_v2/results_v3"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE

N_LNC = {"human": 33668, "mouse": 31467}
N_TRUTH = {"human": 19167, "mouse": 21409, "dmel": 13884, "athal": 26867, "spombe": 5100,
           "tthermophila": 26456, "pfalciparum": 5281, "creinhardtii": 17621}
MINS = ["m100", "m90", "m80", "m70", "m60", "m50"]
RULE = ["m50+g150", "m50+g200", "m50+g250"]
DEFAULT = "m50+g200"
SP_COLOR = {"human": "#1b9e77", "mouse": "#7570b3"}
CFG_C = {"m100": "#bdbdbd", "m50": "#f03b20", DEFAULT: style.TD_C}
CFG_LABEL = {"m100": "100 codons", "m50": "50 codons", DEFAULT: "50 + rule"}
TOOL_LABEL = {"pytd_dev_v2": "TransDecoder", "prodigal": "Prodigal", "gmst": "GeneMarkS-T",
              "borf": "Borf", "cpat": "CPAT", "rnasamba": "RNAsamba", "cpc2": "CPC2"}


def load():
    with open(RES / "short_orf_summary.tsv") as fh:
        return {(r["org"], r["v"], r["cfg"]): r for r in csv.DictReader(fh, delimiter="\t")}


def tradeoff(ax, S, variant, title):
    for sp in ("mouse", "human"):
        c = SP_COLOR[sp]
        a = 1.0 if sp == "human" else 0.55

        def xy(cfgs):
            return ([100 * int(S[(sp, variant, k)]["neg"]) / N_LNC[sp] for k in cfgs],
                    [100 * float(S[(sp, variant, k)]["sens_short"]) for k in cfgs])
        x, y = xy(MINS)
        ax.plot(x, y, "o-", color=c, alpha=a, ms=3, lw=1, label=None)
        x2, y2 = xy(["m50"] + RULE)
        ax.plot(x2, y2, "s--", color=c, alpha=a, ms=3, lw=0.9, mfc="white")
        xd, yd = xy([DEFAULT])
        ax.plot(xd, yd, "*", color=c, alpha=a, ms=9, mec="black", mew=0.5, zorder=6)
        if sp == "human":
            for k, xi, yi in zip(MINS, x, y):
                ax.annotate(k[1:], (xi, yi), xytext=(3, -7), textcoords="offset points", fontsize=5.5,
                            color="#444444")
            for k, xi, yi in zip(RULE, x2[1:], y2[1:]):
                ax.annotate("×" + k.split("g")[1], (xi, yi), xytext=(-4, 4), textcoords="offset points",
                            fontsize=5.5, color="#444444", ha="right")
    ax.set_xlabel("lncRNAs with a prediction (%)")
    ax.set_ylabel("sensitivity, true proteins of 50–99 aa (%)")
    ax.set_xlim(0, 70)
    ax.set_ylim(0, 100)
    ax.set_title(title)


def per_seq(ax, S):
    groups = [(o, "complete.wRand", N_TRUTH[o]) for o in style.ORGS] + \
             [(o, "complete.wLnc", N_LNC[o]) for o in ("human", "mouse")]
    x = np.r_[np.arange(len(style.ORGS)), len(style.ORGS) + 0.7 + np.arange(2)]
    for k, dx in [("m100", -0.25), ("m50", 0), (DEFAULT, 0.25)]:
        y = [int(S[(o, v, k)]["neg"]) / n for o, v, n in groups]
        y = [v if v > 0 else np.nan for v in y]
        ax.bar(x + dx, y, width=0.25, color=CFG_C[k], label=CFG_LABEL[k])
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 4)
    ax.set_xticks(x)
    labs = [r"$\it{" + style.ORG_LABEL[o].replace(" ", r"\ ") + "}$" for o in style.ORGS] + \
           [r"$\it{H.\ sapiens}$", r"$\it{M.\ musculus}$"]
    ax.set_xticklabels(labs, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel("predictions per negative sequence")
    ax.axvline(len(style.ORGS) - 0.3, color="#bbbbbb", lw=0.5, ls=":")
    tr = ax.get_xaxis_transform()
    ax.text(3.5, 1.01, "GC-matched random sequences", transform=tr, ha="center", va="bottom", fontsize=6.5)
    ax.text(len(style.ORGS) + 1.2, 1.01, "lncRNAs", transform=tr, ha="center", va="bottom", fontsize=6.5)
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.0), fontsize=6, handlelength=1, ncol=3,
              columnspacing=0.8)


def era(ax):
    rows, blk = [], None
    with open(RES / "lnc_by_era.human.tsv") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("=="):
                blk = line.split()[1]
            elif blk == "complete.wLnc" and f[0] != "tool":
                rows.append(dict(tool=f[0], era=f[1], n=int(f[2]), called=int(f[3])))
    tools = list(dict.fromkeys(r["tool"] for r in rows))
    tot = {t: sum(r["called"] for r in rows if r["tool"] == t) / sum(r["n"] for r in rows if r["tool"] == t)
           for t in tools}
    tools.sort(key=lambda t: -tot[t])
    x = np.arange(len(tools))
    for e, dx, col, lab in [("pre_v28", -0.18, "#636363", "in GENCODE v28"),
                            ("post_v28", 0.18, "#bdbdbd", "added after v28")]:
        y = [100 * next(r["called"] / r["n"] for r in rows if r["tool"] == t and r["era"] == e) for t in tools]
        ax.bar(x + dx, y, width=0.36, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([TOOL_LABEL[t] for t in tools], rotation=45, ha="right", rotation_mode="anchor")
    for lab in ax.get_xticklabels():
        if lab.get_text() == "TransDecoder":
            lab.set_fontweight("bold")
    ax.set_ylabel("lncRNAs with a prediction (%)")
    ax.set_title("human lncRNAs, by GENCODE release")
    ax.legend(loc="upper right", fontsize=6, handlelength=1)


def main():
    S = load()
    fig = plt.figure(figsize=(style.FULL_W, 128 * style.MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], wspace=0.32, hspace=0.62)
    axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axc, axd = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    tradeoff(axa, S, "complete.wLnc", "complete transcripts + lncRNAs")
    tradeoff(axb, S, "fragments.wLnc", "fragments + lncRNAs")
    L = plt.Line2D
    h = [L([], [], color=SP_COLOR["human"], lw=3), L([], [], color=SP_COLOR["mouse"], lw=3, alpha=0.55),
         L([], [], marker="o", ls="-", color="#666666", ms=3),
         L([], [], marker="s", ls="--", color="#666666", mfc="white", ms=3),
         L([], [], marker="*", ls="", color="#666666", mec="black", mew=0.5, ms=9)]
    l = [r"$\it{H.\ sapiens}$", r"$\it{M.\ musculus}$", "min. 100→50 codons",
         "50 + rule ×150–250", "default"]
    axb.legend(h, l, loc="lower right", fontsize=6, handlelength=1.8, borderaxespad=0.1)
    per_seq(axc, S)
    era(axd)
    for ax, letter in zip([axa, axb, axc, axd], "abcd"):
        style.panel_label(ax, letter, x=-0.14)
    style.save(fig, OUT, "fig4_short_orfs")


if __name__ == "__main__":
    main()
