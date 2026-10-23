#!/usr/bin/env python3
"""Figure 3 (start-codon refinement) and Fig. S2 (alpha x w heatmaps), from committed tables only.

usage: make_fig3.py [out_dir]          (default: figures/fig3_start)

Inputs (repo-relative):
  analyses/start_param_sensitivity/results/{open_end_mix,grid_scores,optima,plateau,
      kozak_score_distribution,gm12878}.tsv
  benchmark_v2/results_v3/summary.tsv   (perl_td_v600 exact F1 only)
Outputs: fig3_start.{pdf,png}, figS2_start_grid.{pdf,png}, fig3_values.tsv (plotted values).
Evidence records: E13 (grid, open ends, GM12878, Kozak distributions), E07 (v6.0.0), E09.

The grid was run with the Kozak weight in nats (0, 5, 10, 25, 50, 100). The code now takes the
weight per bit; w_bits = w_nats * ln 2, so the default 25 per nat is 17.3 per bit (code: 17;
E13 addendum: at most 0.13% of predictions differ). Axes are labelled in bits.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import style  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "analyses/start_param_sensitivity/results"
SUMMARY = REPO / "benchmark_v2/results_v3/summary.tsv"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "figures/fig3_start"

ALPHAS = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.25, 1.5, 2, 3]
W_NATS = [0, 5, 10, 25, 50, 100]
DEFAULT = ("0.4", "25")


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bits(w_nats):
    b = w_nats * math.log(2)
    return f"{b:.0f}" if b >= 10 or b == 0 else f"{b:.1f}"


def cfg(a, w):
    return f"a{a:g}_w{w:g}"


GRID = {(r["org"], r["input"], r["config"]): float(r["f1_exact"]) for r in rows(RES / "grid_scores.tsv")}
NOREF = {(r["org"], r["input"]): float(r["norefine_f1"]) for r in rows(RES / "optima.tsv")}
V600 = {(r["org"], r["variant"]): float(r["f1_exact"]) for r in rows(SUMMARY)
        if r["tool"] == "perl_td_v600" and r["variant"] in ("complete", "fragments")}
PLATEAU = {(r["input"], r["config"]): float(r["mean_f1"]) for r in rows(RES / "plateau.tsv")}
VALUES = []   # (panel, key, value) for fig3_values.tsv


# ---------------------------------------------------------------------------------------------
def panel_open_ends(ax):
    """(a) open 5' ends among no-refine predictions matching a truth stop: truncated vs UTR."""
    mix = {(r["org"], r["input"]): r for r in rows(RES / "open_end_mix.tsv")}
    orgs = style.ORGS[::-1]
    y = np.arange(len(orgs))
    h = 0.38
    for off, inp, edge in ((h / 2 + 0.02, "complete", style.COMPLETE_C), (-h / 2 - 0.02, "fragments", style.FRAG_C)):
        for i, o in enumerate(orgs):
            r = mix[(o, inp)]
            tr, ut = int(r["truncated"]), int(r["utr_extended"])
            n = tr + ut
            ax.barh(i + off, tr / n, h, color="#555555", edgecolor=edge, lw=0.8)
            ax.barh(i + off, ut / n, h, left=tr / n, color="#dddddd", edgecolor=edge, lw=0.8)
            VALUES.append(("a", f"{o}.{inp}.truncated_frac", f"{tr / n:.3f}"))
    ax.set_yticks(y)
    ax.set_yticklabels([style.ORG_LABEL[o] for o in orgs], style="italic")
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of open 5$'$ ends")
    ax.set_title("Open 5$'$ ends (no refinement)", loc="left")
    ax.legend(handles=[Rectangle((0, 0), 1, 1, fc="#555555"), Rectangle((0, 0), 1, 1, fc="#dddddd"),
                       Rectangle((0, 0), 1, 1, fc="white", ec=style.COMPLETE_C, lw=1),
                       Rectangle((0, 0), 1, 1, fc="white", ec=style.FRAG_C, lw=1)],
              labels=["truncated CDS", "UTR extension", "complete (upper bar)", "fragments (lower bar)"],
              loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, handlelength=1.2,
              columnspacing=0.8)
    ax.tick_params(axis="y", length=0)


MARK = [("no refinement", "norefine", dict(marker="o", mfc="white", mec="#444444", ms=4)),
        ("v6.0.0", "v600", dict(marker="s", mfc="#bbbbbb", mec="#777777", ms=3.5)),
        ("$\\alpha$ 0.4, $w$ 0 (no Kozak term)", "w0", dict(marker="^", mfc="white", mec=style.TD_C, ms=4)),
        ("final ($\\alpha$ 0.4, $w$ 17 per bit)", "final", dict(marker="o", mfc=style.TD_C, mec=style.TD_C, ms=4.2))]


def panel_dumbbell(ax, inp, show_y):
    """(b) exact F1 per organism: no refinement, v6.0.0, w=0, final."""
    orgs = style.ORGS[::-1]
    for i, o in enumerate(orgs):
        v = {"norefine": NOREF[(o, inp)], "v600": V600[(o, inp)],
             "w0": GRID[(o, inp, "a0.4_w0")], "final": GRID[(o, inp, "a0.4_w25")]}
        ax.plot([min(v.values()), max(v.values())], [i, i], color="#cccccc", lw=0.8, zorder=1)
        ax.annotate("", xy=(v["final"], i), xytext=(v["norefine"], i),
                    arrowprops=dict(arrowstyle="-", color=style.TD_C, lw=1.4, alpha=0.35), zorder=1)
        for _, k, kw in MARK:
            ax.plot(v[k], i, ls="none", zorder=3, **kw)
            VALUES.append(("b", f"{o}.{inp}.{k}", f"{v[k]:.4f}"))
    ax.set_yticks(range(len(orgs)))
    ax.set_yticklabels([style.ORG_LABEL[o] for o in orgs] if show_y else [], style="italic")
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.40, 1.0) if inp == "complete" else ax.set_xlim(0.58, 0.88)
    ax.set_xlabel("exact-CDS F1")
    ax.set_title("complete" if inp == "complete" else "fragments", loc="left",
                 color=style.COMPLETE_C if inp == "complete" else style.FRAG_C)
    ax.grid(axis="x", color="#eeeeee", lw=0.5)
    ax.set_axisbelow(True)


def panel_alpha(ax):
    """(c) mean exact F1 over organisms vs alpha at w=17 bits (solid) and w=0 (dashed)."""
    for inp, col, lab in (("complete", style.COMPLETE_C, "complete"), ("fragments", style.FRAG_C, "fragments"),
                          ("both", "#222222", "mean of both")):
        for w, ls, lw in ((25, "-", 1.3), (0, "--", 0.8)):
            ys = [PLATEAU[(inp, cfg(a, w))] for a in ALPHAS]
            ax.plot(ALPHAS, ys, ls=ls, lw=lw, color=col, marker="o" if w == 25 else None, ms=2.2,
                    label=lab if w == 25 else None)
            for a, yv in zip(ALPHAS, ys):
                VALUES.append(("c", f"{inp}.w{w}.a{a:g}", f"{yv:.4f}"))
    ax.axvline(0.4, color="#888888", lw=0.6, ls=":")
    ax.set_xscale("function", functions=(np.sqrt, np.square))
    ax.set_xticks([0, 0.1, 0.4, 1, 2, 3])
    ax.set_xticklabels(["0", "0.1", "0.4", "1", "2", "3"])
    ax.set_xlim(0, 3.05)
    ax.set_xlabel("$\\alpha$ (UTR score per nt; square-root scale)")
    ax.set_ylabel("mean exact-CDS F1 (8 organisms)")
    ax.set_title("Trade-off set by $\\alpha$", loc="left")
    h, lab = ax.get_legend_handles_labels()
    h += [plt.Line2D([], [], color="#666666", lw=1.3), plt.Line2D([], [], color="#666666", lw=0.8, ls="--")]
    lab += ["$w$ = 17 per bit", "$w$ = 0"]
    ax.legend(h, lab, loc="lower left", ncol=2, handlelength=1.6, columnspacing=0.6, borderaxespad=0.2)
    ax.set_ylim(0.25, 0.95)
    ax.text(0.42, 0.945, "default", fontsize=6, color="#666666", va="top")


def panel_kozak(ax):
    """(d) Kozak score (bits) at true starts, internal ATGs and the first ATG of 5'-cut fragments."""
    d = {(r["org"], r["set"]): r for r in rows(RES / "kozak_score_distribution.tsv")}
    sets = [("true_start", "annotated start", "#1b7837"), ("internal_ATG", "internal ATG", "#999999"),
            ("frag_firstATG", "first ATG of fragment", style.FRAG_C)]
    stats, pos, cols = [], [], []
    for i, o in enumerate(style.ORGS):
        for j, (s, _, c) in enumerate(sets):
            r = d[(o, s)]
            stats.append(dict(whislo=float(r["p10"]), q1=float(r["p25"]), med=float(r["p50"]),
                              q3=float(r["p75"]), whishi=float(r["p90"]), fliers=[]))
            pos.append(i * 4 + j)
            cols.append(c)
            VALUES.append(("d", f"{o}.{s}.median", r["p50"]))
    bp = ax.bxp(stats, positions=pos, widths=0.75, patch_artist=True, showfliers=False,
                medianprops=dict(color="black", lw=0.8), whiskerprops=dict(lw=0.6),
                capprops=dict(lw=0.6), boxprops=dict(lw=0.5))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
    ax.axhline(0, color="#888888", lw=0.5, ls=":")
    ax.set_xticks([i * 4 + 1 for i in range(len(style.ORGS))])
    ax.set_xticklabels([style.ORG_LABEL[o] for o in style.ORGS], style="italic", rotation=35, ha="right")
    ax.set_ylabel("Kozak score $K$ (bits)")
    ax.set_title("Start context of ATGs (self-trained Kozak model)", loc="left")
    ax.legend(handles=[Rectangle((0, 0), 1, 1, fc=c, alpha=0.75) for _, _, c in sets],
              labels=[lab for _, lab, _ in sets], loc="upper left", ncol=3, handlelength=0.9, handletextpad=0.4,
              bbox_to_anchor=(0, 1.02), columnspacing=0.7)
    ax.set_ylim(-4, 6.2)


def fig3():
    fig = plt.figure(figsize=(style.FULL_W, 150 * style.MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 0.95], hspace=0.62)
    top = gs[0].subgridspec(1, 3, width_ratios=[1.05, 1, 1], wspace=0.12)
    bot = gs[1].subgridspec(1, 2, width_ratios=[0.8, 1.45], wspace=0.28)
    ax_a = fig.add_subplot(top[0])
    panel_open_ends(ax_a)
    ax_b1 = fig.add_subplot(top[1])
    ax_b2 = fig.add_subplot(top[2], sharey=ax_b1)
    panel_dumbbell(ax_b1, "complete", True)
    panel_dumbbell(ax_b2, "fragments", False)
    ax_b1.legend(handles=[plt.Line2D([], [], ls="none", **kw) for _, _, kw in MARK],
                 labels=[lab for lab, _, _ in MARK], loc="upper center", bbox_to_anchor=(1.0, -0.17),
                 ncol=2, handletextpad=0.3, columnspacing=1.0)
    ax_c = fig.add_subplot(bot[0])
    panel_alpha(ax_c)
    ax_d = fig.add_subplot(bot[1])
    panel_kozak(ax_d)
    style.panel_label(ax_a, "a", x=-0.42)
    style.panel_label(ax_b1, "b", x=-0.05)
    style.panel_label(ax_c, "c", x=-0.25)
    style.panel_label(ax_d, "d", x=-0.1)
    style.save(fig, OUT, "fig3_start")


# ---------------------------------------------------------------------------------------------
def heat(ax, mat, xlabels, ylabels, default_xy, title):
    # colour spans the best value minus 5 points, so the near-optimal region stands out
    hi = np.nanmax(mat)
    lo = max(np.nanmin(mat), hi - 5)
    im = ax.imshow(mat, cmap="Blues", origin="lower", aspect="auto", vmin=lo, vmax=hi)
    for (i, j), v in np.ndenumerate(mat):
        frac = (v - lo) / (hi - lo) if hi > lo else 0
        ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.5,
                color="white" if frac > 0.6 else "black",
                fontweight="bold" if v == hi else "normal")
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels)
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Kozak weight $w$ (per bit)")
    ax.set_ylabel("$\\alpha$")
    ax.set_title(title, loc="left")
    dx, dy = default_xy
    ax.add_patch(Rectangle((dx - 0.5, dy - 0.5), 1, 1, fill=False, ec="#e41a1c", lw=1.2))
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    return im


def figS2():
    fig = plt.figure(figsize=(style.FULL_W, 215 * style.MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.38)
    top = gs[0].subgridspec(1, 3, wspace=0.32)
    bot = gs[1].subgridspec(1, 4, wspace=0.45)
    wl = [bits(w) for w in W_NATS]
    titles = {"complete": "Benchmark, complete", "fragments": "Benchmark, fragments",
              "both": "Benchmark, mean of both"}
    axs = []
    for k, inp in enumerate(("complete", "fragments", "both")):
        ax = fig.add_subplot(top[k])
        mat = np.array([[100 * PLATEAU[(inp, cfg(a, w))] for w in W_NATS] for a in ALPHAS])
        heat(ax, mat, wl, [f"{a:g}" for a in ALPHAS], (W_NATS.index(25), ALPHAS.index(0.4)),
             titles[inp] + "\nmean exact-CDS F1 (×100)")
        axs.append(ax)
    gm = {(r["assembly"], r["config"]): r for r in rows(RES / "gm12878.tsv")}
    ga, gw = ["0.2", "0.4", "0.6", "1.0"], [0, 25, 50]
    metrics = [("MANE_sens_exact_given_assembled", "MANE exact | assembled"),
               ("pred_complete_start_precision", "exact-start share")]
    gax = []
    for k, (asm, (col, lab)) in enumerate([(a, m) for m in metrics for a in ("lraa", "stringtie")]):
        ax = fig.add_subplot(bot[k])
        mat = np.array([[100 * float(gm[(asm, f"a{a}_w{w}")][col]) for w in gw] for a in ga])
        nr = 100 * float(gm[(asm, "norefine")][col])
        name = "LRAA" if asm == "lraa" else "StringTie"
        heat(ax, mat, [bits(w) for w in gw], [f"{float(a):g}" for a in ga], (1, 1),
             f"{name}\nno refinement {nr:.1f}")
        for a in ga:
            for w in gw:
                VALUES.append(("S2", f"gm12878.{asm}.{col}.a{a}_w{w}", gm[(asm, f"a{a}_w{w}")][col]))
        gax.append(ax)
    # metric names centred over each pair of GM12878 panels
    for k, (_, lab) in enumerate(metrics):
        l, r = gax[2 * k].get_position(), gax[2 * k + 1].get_position()
        fig.text((l.x0 + r.x1) / 2, l.y1 + 0.045, f"GM12878: {lab} (×100)", ha="center", va="bottom",
                 fontweight="bold")
    style.panel_label(axs[0], "a", x=-0.3, y=1.08)
    style.panel_label(gax[0], "b", x=-0.55, y=1.2)
    style.save(fig, OUT, "figS2_start_grid")


def checks():
    """Recompute numbers the manuscript quotes and print them (see NOTES.md)."""
    out = []
    comp = [GRID[(o, "complete", "a0.4_w25")] - NOREF[(o, "complete")] for o in style.ORGS]
    frag = [GRID[(o, "fragments", "a0.4_w25")] - NOREF[(o, "fragments")] for o in style.ORGS]
    out.append(f"refine gain complete pts: min {100*min(comp):.1f} max {100*max(comp):.1f} "
               f"median {100*np.median(comp):.1f}")
    out.append(f"refine effect fragments pts: min {100*min(frag):.1f} max {100*max(frag):.1f}; "
               f"positive {sum(f > 0.0005 for f in frag)}")
    w0_below = sum(GRID[(o, "fragments", "a0.4_w0")] < NOREF[(o, "fragments")] for o in style.ORGS)
    kz_f = [GRID[(o, "fragments", "a0.4_w25")] - GRID[(o, "fragments", "a0.4_w0")] for o in style.ORGS]
    kz_c = [GRID[(o, "complete", "a0.4_w25")] - GRID[(o, "complete", "a0.4_w0")] for o in style.ORGS]
    out.append(f"w0 fragments below no-refine: {w0_below}/8")
    out.append(f"Kozak term fragments pts: {100*min(kz_f):.1f}..{100*max(kz_f):.1f}; complete "
               f"{100*min(kz_c):.1f}..{100*max(kz_c):.1f} mean {100*np.mean(kz_c):.1f}")
    v6c = [GRID[(o, "complete", "a0.4_w25")] - V600[(o, "complete")] for o in style.ORGS]
    v6f = [GRID[(o, "fragments", "a0.4_w25")] - V600[(o, "fragments")] for o in style.ORGS]
    out.append(f"vs v6.0.0 complete pts {100*min(v6c):.1f}..{100*max(v6c):.1f}; median v600 "
               f"{np.median([V600[(o, 'complete')] for o in style.ORGS]):.3f} final "
               f"{np.median([GRID[(o, 'complete', 'a0.4_w25')] for o in style.ORGS]):.3f}; fragments "
               f"{100*min(v6f):.1f}..{100*max(v6f):.1f}")
    best = PLATEAU[("both", "a0.4_w25")]
    plat = sorted(c for (i, c), v in PLATEAU.items() if i == "both" and v >= best - 0.005)
    out.append(f"'both' cells within 0.5 pt of default ({best:.4f}): {len(plat)}: {','.join(plat)}")
    return out


if __name__ == "__main__":
    fig3()
    figS2()
    with open(OUT / "fig3_values.tsv", "w") as fh:
        fh.write("panel\tkey\tvalue\n")
        for p, k, v in VALUES:
            fh.write(f"{p}\t{k}\t{v}\n")
    print("\n".join(checks()))
