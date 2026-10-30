#!/usr/bin/env python3
"""E36 figure: truncated share (keep-open is exact) by transcript length and by O/L, benchmark vs assemblies.

usage: python 10_figure_v1.py <oracle_by_bin.tsv> <out_prefix>
Panels: (a) stop-ending candidates by L; (b) 3'-open candidates by L; (c) all candidates by O/L.
Series in fixed categorical order (reference palette slots 1-5) with distinct markers (identity is not colour alone).
Writes <out_prefix>.png and .pdf.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

src, out = sys.argv[1:3]
D = pd.read_csv(src, sep="\t")
D = D[D.bin != "all"]
SER = [("bench_fragments", "benchmark fragments", "#2a78d6", "o"),
       ("bench_fragments_mixed", "benchmark mixed", "#eb6834", "s"),
       ("sim_human_d8", "sim. human Trinity d8", "#1baf7a", "^"),
       ("sim_athal_d8", "sim. Arabidopsis Trinity d8", "#eda100", "D"),
       ("gm12878_trinity", "GM12878 Trinity (MANE)", "#e87ba4", "v")]
panels = [("L", "stop", "(a) candidate ORF has a stop: by transcript length L"),
          ("L", "open3", "(b) candidate ORF lacks a stop: by L"),
          ("OL", "all", "(c) all candidates: by O/L")]
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
for ax, (feat, e3, title) in zip(axes, panels):
    x = D[(D.feature == feat) & (D.end3 == e3)]
    bins = list(dict.fromkeys(x.bin))
    for g, lab, col, mk in SER:
        s = x[(x.grp == g)].set_index("bin").reindex(bins)
        n = s.n_keep + s.n_trim
        y = s.trunc_share.where(n >= 30)
        if y.notna().sum() == 0:
            continue
        ax.plot(range(len(bins)), y.values, color=col, marker=mk, ms=8, lw=2, label=lab)
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels([b.replace("[", "").replace(")", "").replace(".0", "").replace(", ", "-").replace("-inf", "+")
                        for b in bins], rotation=40, ha="right", fontsize=9)
    ax.set_title(title, fontsize=10, loc="left")
    ax.grid(axis="y", color="#e5e5e2", lw=0.8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_ylim(-0.03, 1.03)
axes[0].set_ylabel("truncated share = keep-open exact / (keep + trim exact)")
axes[0].legend(fontsize=8, frameon=False, loc="lower left")
fig.text(0.5, -0.02, "bins with < 30 labelled candidates omitted; benchmark complete inputs have share 0 in every bin",
         ha="center", fontsize=8, color="#555555")
fig.tight_layout()
fig.savefig(out + ".png", dpi=150, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
