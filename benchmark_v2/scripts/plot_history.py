#!/usr/bin/env python3
"""Supplementary history panel: accuracy of key TransDecoder releases on the benchmark inputs.

usage: plot_history.py <history_scores.tsv> <out_prefix>
Labels (in release order): r2012, v2.0, v4.1.0, v5.0.0, v5.7.1, pytd_release (current rule; Perl devel
gives identical predictions), pytd_decoy_q0.05 (this paper's default). Median over organisms.
Writes <out_prefix>.{pdf,png}.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ORDER = [("r2012", "2012\n(SourceForge)"), ("v2.0", "v2.0\n(2015)"), ("v4.1.0", "v4.1\n(2017)"),
         ("v5.0.0", "v5.0\n(2017)"), ("v5.7.1", "v5.7.1\n(2023)"), ("pytd_release", "v7\nlength rule"),
         ("pytd_decoy_q0.05", "v7\ndecoy FDR")]
VARIANTS = [("complete.wRand", "complete + random"), ("fragments.wRand", "fragments + random"),
            ("complete.wLnc", "complete + lncRNA"), ("fragments.wLnc", "fragments + lncRNA")]
METRICS = [("sens_stop", "sensitivity (stop)"), ("sens_exact", "sensitivity (exact CDS)"),
           ("neg_rate", "negatives with a prediction (%)")]


def main():
    d = pd.read_csv(sys.argv[1], sep="\t")
    m = d.groupby(["variant", "label"])[[k for k, _ in METRICS]].median()
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    x = range(len(ORDER))
    for ax, (k, name) in zip(axes, METRICS):
        for j, (v, vname) in enumerate(VARIANTS):
            y = [m.loc[(v, lab), k] if (v, lab) in m.index else float("nan") for lab, _ in ORDER]
            if k == "neg_rate":
                y = [100 * t for t in y]
            ax.plot(x, y, "-o", color=f"C{j}", ms=4, label=vname)
        ax.set_xticks(list(x))
        ax.set_xticklabels([t for _, t in ORDER], fontsize=7)
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(sys.argv[2] + ".pdf")
    fig.savefig(sys.argv[2] + ".png", dpi=150)


if __name__ == "__main__":
    main()
