#!/usr/bin/env python3
"""Sensitivity against calls on negatives for the PyTD stringency sweep, with other tools as points.

Input: score_negatives.py output (org, variant, label, ..., sens_stop, sens_exact, neg_rate).
Labels starting with 's100_' or 'sall_' are sweep settings; other labels are drawn as single points.
Per variant: the median over organisms (wRand: 8 organisms; wLnc: human and mouse) of each metric.

usage: plot_stringency.py <scores.tsv> <out_prefix> [label,label,... to draw as points]
Writes <out_prefix>.medians.tsv and <out_prefix>.{pdf,png}.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

VARIANTS = [("complete.wRand", "complete + random"), ("fragments.wRand", "fragments + random"),
            ("complete.wLnc", "complete + lncRNA (human, mouse)"),
            ("fragments.wLnc", "fragments + lncRNA (human, mouse)")]
FAMILIES = [("s100_", "short ORFs (<100 codons) only", "tab:blue"),
            ("sall_", "all ORFs", "tab:orange")]
MARK = "ovs^D<>pPhX*"


def contrast(label):
    return float(label.split("_c")[1].split("_")[0])


def main():
    df = pd.read_csv(sys.argv[1], sep="\t")
    out = sys.argv[2]
    points = sys.argv[3].split(",") if len(sys.argv) > 3 else []
    med = (df.groupby(["variant", "label"])[["sens_stop", "sens_exact", "precision_stop", "neg_rate"]]
           .median().reset_index())
    med["n_orgs"] = df.groupby(["variant", "label"]).size().values
    med.to_csv(out + ".medians.tsv", sep="\t", index=False, float_format="%.4f")
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.5))
    for j, (v, title) in enumerate(VARIANTS):
        m = med[med.variant == v]
        for i, metric in enumerate(["sens_stop", "sens_exact"]):
            ax = axes[i, j]
            for pre, name, col in FAMILIES:
                f = m[m.label.str.startswith(pre) & ~m.label.str.endswith("_noauto")].copy()
                f["c"] = f.label.map(contrast)
                f = f.sort_values("c")
                ax.plot(100 * f.neg_rate, f[metric], "-o", color=col, ms=3, label=name)
                for _, r in f.iterrows():
                    ax.annotate(f"{r.c:g}", (100 * r.neg_rate, r[metric]), fontsize=6, color=col,
                                xytext=(3, 2), textcoords="offset points")
            na = m[m.label.str.endswith("_noauto")]
            ax.plot(100 * na.neg_rate, na[metric], "x", color="tab:orange", label="all ORFs, no length rule")
            d = m[m.label == "s100_c200"]
            ax.plot(100 * d.neg_rate, d[metric], "*", color="k", ms=12, label="final defaults")
            for k, p in enumerate(points):
                q = m[m.label == p]
                if len(q):
                    ax.plot(100 * q.neg_rate, q[metric], MARK[k % len(MARK)], color=f"C{(k + 2) % 10}",
                            ms=6, label=p)
            ax.set_xscale("symlog", linthresh=1)
            ax.set_xlabel("negatives with a prediction (%)")
            ax.set_ylabel({"sens_stop": "sensitivity (stop)", "sens_exact": "sensitivity (exact CDS)"}[metric])
            if i == 0:
                ax.set_title(title, fontsize=10)
            ax.grid(alpha=0.3)
    h, lab = axes[0, 0].get_legend_handles_labels()
    for ax in axes.flat:  # union of legend entries
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh)
                lab.append(ll)
    fig.legend(h, lab, loc="lower center", ncol=6, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png", dpi=150)


if __name__ == "__main__":
    main()
