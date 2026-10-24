#!/usr/bin/env python3
"""Sensitivity against calls on negatives for the short-ORF rule families (fresh full runs), with
other tools as points.

Input: score_negatives.py output. Labels '<rule>_c<C>' (rule: cliff, fixed, inverse, taper) are drawn
as one curve per rule over C; other labels given on the command line are drawn as single points.
Per variant: the median over organisms of each metric.

usage: plot_rules.py <scores.tsv> <out_prefix> [label,label,... to draw as points]
Writes <out_prefix>.medians.tsv and <out_prefix>.{pdf,png}.
"""
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

VARIANTS = [("complete.wRand", "complete + random"), ("fragments.wRand", "fragments + random"),
            ("complete.wLnc", "complete + lncRNA (human, mouse)"),
            ("fragments.wLnc", "fragments + lncRNA (human, mouse)")]
RULES = [("cliff", "cliff: C below 100 codons only (release)", "k"),
         ("fixed", "fixed: C at all lengths", "tab:orange"),
         ("inverse", "inverse: C x 100/L", "tab:green"),
         ("taper", "taper: C at 100 codons -> 0 at auto-accept", "tab:blue")]
RX = re.compile(r"^(cliff|fixed|inverse|taper)_c([0-9.]+)$")
MARK = "ovs^D<>pPhX*"


def main():
    df = pd.read_csv(sys.argv[1], sep="\t")
    out = sys.argv[2]
    points = sys.argv[3].split(",") if len(sys.argv) > 3 else []
    g = df.groupby(["variant", "label"])
    med = g[["sens_stop", "sens_exact", "precision_stop", "neg_rate"]].median().reset_index()
    med["n_orgs"] = g.size().values
    med.to_csv(out + ".medians.tsv", sep="\t", index=False, float_format="%.4f")
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.5))
    for j, (v, title) in enumerate(VARIANTS):
        m = med[med.variant == v]
        for i, metric in enumerate(["sens_stop", "sens_exact"]):
            ax = axes[i, j]
            for rule, name, col in RULES:
                f = m[m.label.str.match(rf"^{rule}_c[0-9.]+$")].copy()
                if not len(f):
                    continue
                f["c"] = f.label.map(lambda x: float(RX.match(x).group(2)))
                f = f.sort_values("c")
                ax.plot(100 * f.neg_rate, f[metric], "-o", color=col, ms=3, label=name)
                for _, r in f.iterrows():
                    ax.annotate(f"{r.c:g}", (100 * r.neg_rate, r[metric]), fontsize=6, color=col,
                                xytext=(3, 2), textcoords="offset points")
            for k, p in enumerate(points):
                q = m[m.label == p]
                if len(q):
                    ax.plot(100 * q.neg_rate, q[metric], MARK[k % len(MARK)], color=f"C{(k + 4) % 10}",
                            ms=6, label=p)
            ax.set_xscale("symlog", linthresh=1)
            ax.set_xlabel("negatives with a prediction (%)")
            ax.set_ylabel({"sens_stop": "sensitivity (stop)", "sens_exact": "sensitivity (exact CDS)"}[metric])
            if i == 0:
                ax.set_title(title, fontsize=10)
            ax.grid(alpha=0.3)
    h, lab = [], []
    for ax in axes.flat:  # union of legend entries
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    fig.legend(h, lab, loc="lower center", ncol=5, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=150)


if __name__ == "__main__":
    main()
