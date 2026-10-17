#!/usr/bin/env python3
"""Figures from score_all.py output.

usage: plot_summary.py <summary.tsv> <outdir>

Writes, for each input setting (complete / fragments, stranded / unstranded, with
random negatives), a grid of sensitivity / precision / F1 at stop and exact level by
tool and organism, plus a runtime figure.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ORG_ORDER = ["human", "mouse", "dmel", "athal", "spombe", "creinhardtii", "pfalciparum", "tthermophila"]
ORG_LABEL = {"human": "H. sapiens", "mouse": "M. musculus", "dmel": "D. melanogaster",
             "athal": "A. thaliana", "spombe": "S. pombe", "creinhardtii": "C. reinhardtii",
             "pfalciparum": "P. falciparum", "tthermophila": "T. thermophila"}
TOOL_ORDER = ["pytd", "pytd_homology", "pytd_norefine", "perl_td", "perl_td_v600", "perl_td_v571",
              "codan", "prodigal", "longest_orf", "all_orfs"]


def main():
    summary, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    d = pd.read_csv(summary, sep="\t")
    d["org_label"] = d["org"].map(ORG_LABEL)
    orgs = [ORG_LABEL[o] for o in ORG_ORDER if o in set(d["org"])]
    tools = [t for t in TOOL_ORDER if t in set(d["tool"])]

    for variant, dv in d.groupby("variant"):
        long = dv.melt(id_vars=["org_label", "tool"],
                       value_vars=["sensitivity_stop", "precision_stop", "f1_stop",
                                   "sensitivity_exact", "precision_exact", "f1_exact"],
                       var_name="metric", value_name="value")
        g = sns.catplot(data=long, x="org_label", y="value", hue="tool", col="metric", col_wrap=3,
                        kind="bar", order=orgs, hue_order=tools, height=3.2, aspect=1.6,
                        col_order=["sensitivity_stop", "precision_stop", "f1_stop",
                                   "sensitivity_exact", "precision_exact", "f1_exact"])
        g.set_xticklabels(rotation=45, ha="right", fontstyle="italic")
        g.set_axis_labels("", "")
        g.set(ylim=(0, 1))
        g.fig.suptitle(variant, y=1.02)
        g.savefig(os.path.join(outdir, f"accuracy.{variant}.pdf"), bbox_inches="tight")
        plt.close(g.fig)

    rt = d[d["runtime_sec"].notna() & (d["variant"] == "complete")]
    if len(rt):
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.barplot(data=rt, x="org_label", y="runtime_sec", hue="tool", order=orgs,
                    hue_order=[t for t in tools if t in set(rt["tool"])], ax=ax)
        ax.set_yscale("log")
        ax.set_xlabel("")
        ax.set_ylabel("run time (s)")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontstyle="italic")
        fig.savefig(os.path.join(outdir, "runtime.complete.pdf"), bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
