#!/usr/bin/env python3
"""Figure 2: accuracy on eight eukaryotes.

usage: make_fig2.py <summary.tsv> [out_dir]     (out_dir default: figures/fig2_benchmark)

summary.tsv is the output of benchmark_v2/scripts/score_all.py (E07). Panels, stranded inputs only:
  (a) stop-level F1, complete      (b) exact F1, complete
  (c) stop-level precision, complete + GC-matched random sequences (complete.wRand)
  (d) stop-level F1, fragments     (e) exact F1, fragments
  (f) run time on complete inputs, median over organisms (range as whiskers)
Missing tool x organism cells (e.g. CodAn runs still in progress) are drawn grey and marked
"pending"; cells of tools run with the standard code on T. thermophila are hatched.
Also writes fig2_values.tsv (every plotted value) and fig2_medians.tsv (medians over organisms,
for checking Table tab:bench).
"""
import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

SUMMARY = Path(sys.argv[1])
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE

# (summary tool name, label); pytd_dev_v2 == perl_td_dev_v2 on all 84 inputs (E08)
TOOLS = [("pytd_dev_v2", "TransDecoder"), ("perl_td_v600", "TransDecoder v6.0.0"),
         ("prodigal", "Prodigal"), ("gmst", "GeneMarkS-T"), ("codan", "CodAn"),
         ("cpat", "CPAT"), ("cpc2", "CPC2"), ("rnasamba", "RNAsamba"), ("borf", "Borf"),
         ("longest_orf", "longest ORF")]
GROUP_BREAKS = [2, 5, 9]          # lines after TransDecoder, gene finders, classifiers
WRONG_CODE = {"gmst", "codan", "cpc2", "rnasamba", "borf"}   # standard code on T. thermophila
RUNTIME_TOOLS = [("pytd_dev_v2", "PyTransDecoder"), ("perl_td_dev_v2", "TransDecoder (Perl)"),
                 ("prodigal", "Prodigal"), ("gmst", "GeneMarkS-T"), ("codan", "CodAn"),
                 ("cpat", "CPAT"), ("cpc2", "CPC2"), ("rnasamba", "RNAsamba"), ("borf", "Borf"),
                 ("angel", "ANGEL")]
ORG_SHORT = {"pfalciparum": "P. fal.", "tthermophila": "T. the.", "spombe": "S. pom.",
             "athal": "A. tha.", "dmel": "D. mel.", "mouse": "M. mus.", "human": "H. sap.",
             "creinhardtii": "C. rei."}
CMAP = LinearSegmentedColormap.from_list("f1", ["#f7fbff", "#c6dbef", "#6baed6", "#2171b5",
                                                "#08306b"])
VMIN, VMAX = 0.4, 1.0

rows = {(r["org"], r["variant"], r["tool"]): r
        for r in csv.DictReader(open(SUMMARY), delimiter="\t")}


def value(org, variant, tool, metric):
    r = rows.get((org, variant, tool))
    if r is None or r.get(metric, "") in ("", "NA", "nan"):
        return None
    return float(r[metric])


def heatmap(ax, variant, metric, title, show_rows):
    n_t, n_o = len(TOOLS), len(style.ORGS)
    mat = np.full((n_t, n_o), np.nan)
    for i, (tool, _) in enumerate(TOOLS):
        for j, org in enumerate(style.ORGS):
            v = value(org, variant, tool, metric)
            if v is not None:
                mat[i, j] = v
    ax.imshow(np.ma.masked_invalid(mat), cmap=CMAP, vmin=VMIN, vmax=VMAX, aspect="auto")
    for i, (tool, _) in enumerate(TOOLS):
        for j, org in enumerate(style.ORGS):
            v = mat[i, j]
            wrong = org == "tthermophila" and tool in WRONG_CODE
            if np.isnan(v):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc="#dddddd", ec="white",
                                           lw=0.5))
                ax.text(j, i, "pend.", ha="center", va="center", fontsize=4.3, color="#555555")
                continue
            if wrong:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc="white", ec="#aaaaaa",
                                           hatch="/////", lw=0))
                continue
            txt = "1.0" if round(v, 2) >= 1 else f"{v:.2f}"[1:]
            dark = (v - VMIN) / (VMAX - VMIN) > 0.55
            ax.text(j, i, txt, ha="center", va="center", fontsize=4.6,
                    color="white" if dark else "#222222")
    for b in GROUP_BREAKS:
        ax.axhline(b - 0.5, color="white", lw=1.6)
    ax.set_xticks(range(n_o))
    ax.set_xticklabels([ORG_SHORT[o] for o in style.ORGS], rotation=55, ha="right",
                       rotation_mode="anchor", fontsize=5.6, style="italic")
    ax.set_yticks(range(n_t))
    ax.set_yticklabels([t[1] for t in TOOLS] if show_rows else [], fontsize=6)
    ax.tick_params(length=0, pad=1.5)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, pad=3)
    return mat


def runtime(ax):
    labels, med, lo, hi, n = [], [], [], [], []
    for tool, label in RUNTIME_TOOLS:
        vals = [value(o, "complete", tool, "runtime_sec") for o in style.ORGS]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        labels.append(label)
        med.append(statistics.median(vals))
        lo.append(min(vals))
        hi.append(max(vals))
        n.append(len(vals))
    y = np.arange(len(labels))[::-1]
    for yi, m, a, b, k, lab in zip(y, med, lo, hi, n, labels):
        c = style.TD_C if "TransDecoder" in lab else "#6b6b6b"
        ax.plot([a, b], [yi, yi], color=c, lw=0.8, alpha=0.6)
        ax.plot(m, yi, "o", color=c, ms=3.5)
        note = f"{m:,.0f} s" + ("" if k == 8 else f"  (n = {k})")
        if lab == "ANGEL":
            note = f"{m/3600:.0f} h  (A. thaliana only)"
        elif m > 3600:
            note = f"{m/3600:.1f} h" + ("" if k == 8 else f"  (n = {k})")
        ax.text(max(b, m) * 1.35, yi, note, va="center", fontsize=5.4, color="#444444")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xscale("log")
    ax.set_xlim(5, 3e6)
    ax.set_xticks([10, 1e3, 1e5])
    ax.set_xticklabels(["10 s", "10$^3$ s", "10$^5$ s"])
    ax.set_xlabel("wall-clock time (log scale)")
    ax.tick_params(axis="y", length=0)
    ax.set_title("run time, complete inputs\n(median and range over organisms)", pad=3)
    return list(zip(labels, med, lo, hi, n))


def main():
    fig = plt.figure(figsize=(style.FULL_W, 128 * style.MM))
    w, h = 0.215, 0.33
    xs = [0.115, 0.36, 0.605]
    y_top, y_bot = 0.585, 0.105
    ax = {}
    ax["a"] = fig.add_axes([xs[0], y_top, w, h])
    ax["b"] = fig.add_axes([xs[1], y_top, w, h])
    ax["c"] = fig.add_axes([xs[2], y_top, w, h])
    ax["d"] = fig.add_axes([xs[0], y_bot, w, h])
    ax["e"] = fig.add_axes([xs[1], y_bot, w, h])
    ax["f"] = fig.add_axes([0.765, y_bot, 0.14, h])
    mats = {
        ("complete", "f1_stop"): heatmap(ax["a"], "complete", "f1_stop", "complete: stop-level F1", True),
        ("complete", "f1_exact"): heatmap(ax["b"], "complete", "f1_exact", "complete: exact F1", False),
        ("complete.wRand", "precision_stop"): heatmap(ax["c"], "complete.wRand", "precision_stop",
                                                      "complete + random: precision", False),
        ("fragments", "f1_stop"): heatmap(ax["d"], "fragments", "f1_stop", "fragments: stop-level F1", True),
        ("fragments", "f1_exact"): heatmap(ax["e"], "fragments", "f1_exact", "fragments: exact F1", False),
    }
    rt = runtime(ax["f"])
    # colour bar + hatch legend
    cax = fig.add_axes([0.845, y_top + 0.14, 0.012, h - 0.14])
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(VMIN, VMAX))
    cb = fig.colorbar(sm, cax=cax)
    cb.set_ticks([0.4, 0.6, 0.8, 1.0])
    cb.set_ticklabels(["$\\leq$0.4", "0.6", "0.8", "1.0"])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.6, length=2)
    cb.set_label("F1 or precision", fontsize=6)
    lax = fig.add_axes([0.845, y_top - 0.07, 0.12, 0.19])
    lax.axis("off")
    lax.add_patch(plt.Rectangle((0, 0.62), 0.1, 0.13, fc="white", ec="#aaaaaa", hatch="/////",
                                lw=0.5, transform=lax.transAxes))
    lax.text(0.15, 0.685, "run with the standard code\n(no ciliate code available)",
             fontsize=5.4, va="center", transform=lax.transAxes)
    lax.add_patch(plt.Rectangle((0, 0.25), 0.1, 0.13, fc="#dddddd", ec="white",
                                transform=lax.transAxes))
    lax.text(0.15, 0.315, "run pending", fontsize=5.4, va="center", transform=lax.transAxes)
    for k in "abcde":
        style.panel_label(ax[k], k, x=-0.08 if k not in "ad" else -0.55, y=1.04)
    style.panel_label(ax["f"], "f", x=-1.25, y=1.04)
    style.save(fig, OUT, "fig2_benchmark")

    # values and medians for checking
    with open(OUT / "fig2_values.tsv", "w") as fh:
        fh.write("variant\tmetric\ttool\torg\tvalue\n")
        for (variant, metric), mat in mats.items():
            for i, (tool, _) in enumerate(TOOLS):
                for j, org in enumerate(style.ORGS):
                    v = mat[i, j]
                    fh.write(f"{variant}\t{metric}\t{tool}\t{org}\t{'NA' if np.isnan(v) else f'{v:.4f}'}\n")
    with open(OUT / "fig2_medians.tsv", "w") as fh:
        fh.write("variant\tmetric\ttool\tn_orgs\tmedian\tmin\tmax\n")
        for (variant, metric), mat in mats.items():
            for i, (tool, _) in enumerate(TOOLS):
                v = mat[i][~np.isnan(mat[i])]
                if len(v):
                    fh.write(f"{variant}\t{metric}\t{tool}\t{len(v)}\t{np.median(v):.4f}\t"
                             f"{v.min():.4f}\t{v.max():.4f}\n")
        for label, m, a, b, k in rt:
            fh.write(f"complete\truntime_sec\t{label}\t{k}\t{m:.1f}\t{a:.1f}\t{b:.1f}\n")


if __name__ == "__main__":
    main()
