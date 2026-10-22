#!/usr/bin/env python3
"""Summarise the start-refinement parameter grid (run_grid.py, run_gm12878.sh).

usage: summarise.py [--grid $P/work/start_param_grid] [--out results/]

Writes to --out:
  grid_scores.tsv       one row per org x input x config (score_predictions.py metrics)
  checks.tsv            sanity: a0.4_w25 preds.tsv vs the pytd_dev_v2 run; stop-level invariance
  optima.tsv            per org x input: default F1, best cell(s), gap, rank, plateau size
  plateau.tsv           per input x cell: mean F1, worst-case and mean regret over organisms,
                        number of organisms within 0.5 / 1.0 F1 points of their optimum
  loo.tsv               leave-one-organism-out choice of (alpha, w) and its held-out F1
  gm12878.tsv           GM12878 MANE / GENCODE metrics per assembly x config (if run)
  heatmap_f1_exact_<input>.{png,pdf}, heatmap_sens_exact_<input>.{png,pdf},
  heatmap_summary.{png,pdf}, gm12878_heatmap.{png,pdf}
"""
import argparse
import collections
import filecmp
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
ORGS = "athal creinhardtii dmel human mouse pfalciparum spombe tthermophila".split()
INPUTS = ["complete", "fragments"]
DEFAULT = (0.4, 25.0)
CFG_RE = re.compile(r"^a([0-9.]+)_w([0-9.]+)$")
# sequential blue ramp (light -> dark = low -> high), from the dataviz reference palette
SEQ = LinearSegmentedColormap.from_list("seq", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
INK, MUTED = "#1f1f1e", "#6b6b66"


def parse_cfg(cfg):
    m = CFG_RE.match(cfg)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


def read_kv_tsv(path):
    with open(path) as fh:
        hdr, vals = fh.read().strip().split("\n")
    return dict(zip(hdr.split("\t"), vals.split("\t")))


def collect(grid):
    rows = []
    for org in ORGS:
        for inp in INPUTS:
            d = f"{grid}/runs/{org}/{inp}"
            if not os.path.isdir(d):
                continue
            for cfg in sorted(os.listdir(d)):
                sf = f"{d}/{cfg}/score.tsv"
                if not os.path.exists(sf):
                    continue
                r = read_kv_tsv(sf)
                a, w = parse_cfg(cfg)
                rt = open(f"{d}/{cfg}/runtime_sec.txt").read().strip()
                rows.append(dict(org=org, input=inp, config=cfg, alpha=a, weight=w,
                                 n_truth=int(r["n_truth"]), n_pred=int(r["n_pred"]),
                                 tp_stop=int(r["tp_stop"]), tp_exact=int(r["tp_exact"]),
                                 sens_exact=float(r["sensitivity_exact"]),
                                 prec_exact=float(r["precision_exact"]), f1_exact=float(r["f1_exact"]),
                                 f1_stop=float(r["f1_stop"]), runtime_sec=float(rt)))
    return rows


def write_tsv(path, rows, cols=None):
    if not rows:
        return
    cols = cols or list(rows[0])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("NA" if r.get(c) is None else
                               (f"{r[c]:.4f}" if isinstance(r[c], float) and c not in ("alpha", "weight")
                                else (f"{r[c]:g}" if isinstance(r[c], float) else str(r[c])))
                               for c in cols) + "\n")


def f1_from_counts(tp, n_truth, n_pred):
    s, p = tp / n_truth, tp / n_pred
    return 2 * s * p / (s + p) if s + p else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default=f"{P}/work/start_param_grid")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "results"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows = collect(a.grid)
    write_tsv(f"{a.out}/grid_scores.tsv", rows,
              "org input config alpha weight n_truth n_pred tp_stop tp_exact sens_exact prec_exact f1_exact f1_stop runtime_sec".split())
    by = collections.defaultdict(dict)          # (org, input) -> config -> row
    for r in rows:
        by[(r["org"], r["input"])][r["config"]] = r
    alphas = sorted({r["alpha"] for r in rows if r["alpha"] is not None})
    weights = sorted({r["weight"] for r in rows if r["weight"] is not None})
    cells = [(al, w) for al in alphas for w in weights]
    cname = lambda c: f"a{c[0]:g}_w{c[1]:g}"  # noqa: E731
    dname = cname(DEFAULT)

    # ---- checks
    checks = []
    for (org, inp), cfgs in sorted(by.items()):
        ref = f"{P}/work/bench_v3/runs/{org}/{inp}/pytd_dev_v2/preds.tsv"
        mine = f"{a.grid}/runs/{org}/{inp}/{dname}/preds.tsv"
        same = os.path.exists(mine) and filecmp.cmp(ref, mine, shallow=False)
        n_pred = {r["n_pred"] for r in cfgs.values()}
        tp_stop = {r["tp_stop"] for r in cfgs.values()}
        n_pred_r = {r["n_pred"] for c, r in cfgs.items() if c != "norefine"}
        tp_stop_r = {r["tp_stop"] for c, r in cfgs.items() if c != "norefine"}
        checks.append(dict(org=org, input=inp, n_configs=len(cfgs),
                           default_preds_identical_to_pytd_dev_v2=same,
                           n_pred_values=",".join(map(str, sorted(n_pred))),
                           tp_stop_values=",".join(map(str, sorted(tp_stop))),
                           stop_invariant_over_refining_configs=len(n_pred_r) == 1 and len(tp_stop_r) == 1,
                           stop_invariant_incl_norefine=len(n_pred) == 1 and len(tp_stop) == 1))
    write_tsv(f"{a.out}/checks.tsv", checks)

    # ---- optima per org x input
    opt_rows = []
    for inp in INPUTS:
        for org in ORGS:
            cfgs = by.get((org, inp), {})
            grid_rows = [cfgs[cname(c)] for c in cells if cname(c) in cfgs]
            if len(grid_rows) < len(cells) or dname not in cfgs:
                continue
            best = max(r["f1_exact"] for r in grid_rows)
            bests = [r["config"] for r in grid_rows if r["f1_exact"] == best]
            dflt = cfgs[dname]["f1_exact"]
            rank = 1 + sum(r["f1_exact"] > dflt for r in grid_rows)
            opt_rows.append(dict(
                org=org, input=inp, default_f1=dflt, default_sens=cfgs[dname]["sens_exact"],
                default_tp_exact=cfgs[dname]["tp_exact"], best_f1=best, best_configs=",".join(bests),
                best_tp_exact=max(r["tp_exact"] for r in grid_rows),
                gap_f1_points=round(100 * (best - dflt), 2),
                gap_tp_exact=max(r["tp_exact"] for r in grid_rows) - cfgs[dname]["tp_exact"],
                default_rank=f"{rank}/{len(grid_rows)}",
                n_cells_within_0p2pt=sum(best - r["f1_exact"] <= 0.002 + 1e-9 for r in grid_rows),
                n_cells_within_0p5pt=sum(best - r["f1_exact"] <= 0.005 + 1e-9 for r in grid_rows),
                n_cells_within_1pt=sum(best - r["f1_exact"] <= 0.010 + 1e-9 for r in grid_rows),
                worst_grid_f1=min(r["f1_exact"] for r in grid_rows),
                worst_grid_config=min(grid_rows, key=lambda r: r["f1_exact"])["config"],
                norefine_f1=cfgs["norefine"]["f1_exact"] if "norefine" in cfgs else None,
                a0_w0_f1=cfgs[cname((0.0, 0.0))]["f1_exact"] if cname((0.0, 0.0)) in cfgs else None))
    write_tsv(f"{a.out}/optima.tsv", opt_rows)

    # ---- plateau summaries
    plat_rows = []
    full = {inp: [o for o in ORGS if all(cname(c) in by.get((o, inp), {}) for c in cells)] for inp in INPUTS}
    plat = {}
    full["both"] = [o for o in ORGS if all(o in full[i] for i in INPUTS)]

    def f1c(o, scope, c, lam=0.5):
        """exact F1 of cell c; scope 'both' = lam * complete + (1 - lam) * fragments"""
        if scope != "both":
            return by[(o, scope)][cname(c)]["f1_exact"]
        return lam * by[(o, "complete")][cname(c)]["f1_exact"] + (1 - lam) * by[(o, "fragments")][cname(c)]["f1_exact"]

    for inp in INPUTS + ["both"]:
        orgs = full[inp]
        if not orgs:
            continue
        best = {o: max(f1c(o, inp, c) for c in cells) for o in orgs}
        for c in cells:
            f = [f1c(o, inp, c) for o in orgs]
            reg = [best[o] - f1c(o, inp, c) for o in orgs]
            r = dict(input=inp, config=cname(c), alpha=c[0], weight=c[1], n_orgs=len(orgs),
                     mean_f1=sum(f) / len(f),
                     mean_regret_points=100 * sum(reg) / len(reg), max_regret_points=100 * max(reg),
                     n_orgs_within_0p5pt=sum(x <= 0.005 + 1e-9 for x in reg),
                     n_orgs_within_1pt=sum(x <= 0.010 + 1e-9 for x in reg))
            plat_rows.append(r)
            plat[(inp, c)] = r
    write_tsv(f"{a.out}/plateau.tsv", plat_rows)

    # ---- complete:fragments mixing weight -> best cell (mean over organisms)
    mix_rows = []
    if full["both"] and DEFAULT in cells:
        orgs = full["both"]
        for k in range(11):
            lam = k / 10
            sc = {c: sum(f1c(o, "both", c, lam) for o in orgs) / len(orgs) for c in cells}
            bc = max(cells, key=lambda c: (sc[c], -abs(c[0] - DEFAULT[0])))
            near = [cname(c) for c in cells if sc[bc] - sc[c] <= 0.0025 + 1e-9]
            mix_rows.append(dict(weight_complete=lam, best=cname(bc), best_mean_f1=sc[bc],
                                 default_mean_f1=sc[DEFAULT],
                                 default_gap_points=round(100 * (sc[bc] - sc[DEFAULT]), 2),
                                 default_rank=f"{1 + sum(sc[c] > sc[DEFAULT] for c in cells)}/{len(cells)}",
                                 cells_within_0p25pt=",".join(near)))
    write_tsv(f"{a.out}/mix.tsv", mix_rows)

    # ---- leave-one-organism-out
    loo_rows = []
    for scope in INPUTS + ["both"]:
        inps = INPUTS if scope == "both" else [scope]
        orgs = [o for o in ORGS if all(o in full[i] for i in inps)]

        def mean_f1(c, train):
            v = [by[(o, i)][cname(c)]["f1_exact"] for o in train for i in inps]
            return sum(v) / len(v)

        if len(orgs) < 2:
            continue
        for held in orgs:
            train = [o for o in orgs if o != held]
            scores = {c: mean_f1(c, train) for c in cells}
            top = max(scores.values())
            chosen = sorted([c for c in cells if scores[c] == top])
            c = chosen[0]
            held_f1 = mean_f1(c, [held])
            held_def = mean_f1(DEFAULT, [held])
            held_oracle = max(mean_f1(x, [held]) for x in cells)
            loo_rows.append(dict(scope=scope, held_out=held, chosen=cname(c),
                                 n_tied=len(chosen), train_mean_f1=top,
                                 train_mean_f1_default=mean_f1(DEFAULT, train),
                                 heldout_f1_chosen=held_f1, heldout_f1_default=held_def,
                                 heldout_f1_oracle=held_oracle,
                                 chosen_minus_default_points=round(100 * (held_f1 - held_def), 2),
                                 oracle_minus_chosen_points=round(100 * (held_oracle - held_f1), 2)))
    write_tsv(f"{a.out}/loo.tsv", loo_rows)

    # ---- figures
    def panel(ax, inp, org, metric, title):
        cfgs = by[(org, inp)]
        M = [[cfgs[cname((al, w))][metric] for w in weights] for al in alphas]
        best = max(max(r) for r in M)
        lo = best - 0.05
        ax.imshow([[max(v, lo) for v in r] for r in M], cmap=SEQ, vmin=lo, vmax=best, aspect="auto",
                  origin="lower")
        for i, al in enumerate(alphas):
            for j, w in enumerate(weights):
                v = M[i][j]
                dark = (max(v, lo) - lo) / (best - lo + 1e-12) > 0.55
                ax.text(j, i, f"{100 * v:.1f}", ha="center", va="center", fontsize=6.5,
                        color="white" if dark else INK,
                        fontweight="bold" if v == best else "normal")
        di, dj = alphas.index(DEFAULT[0]), weights.index(DEFAULT[1])
        ax.add_patch(Rectangle((dj - 0.5, di - 0.5), 1, 1, fill=False, ec="#e34948", lw=2))
        ax.set_xticks(range(len(weights)), [f"{w:g}" for w in weights], fontsize=7)
        ax.set_yticks(range(len(alphas)), [f"{al:g}" for al in alphas], fontsize=7)
        extra = ""
        if "norefine" in cfgs:
            extra = f"  (no refinement {100 * cfgs['norefine'][metric]:.1f})"
        ax.set_title(f"{title}{extra}", fontsize=8, color=INK)
        for s in ax.spines.values():
            s.set_visible(False)

    label = {"f1_exact": "exact-start F1", "sens_exact": "exact-start sensitivity"}
    for metric in ("f1_exact", "sens_exact"):
        for inp in INPUTS:
            orgs = full[inp]
            if not orgs:
                continue
            fig, axes = plt.subplots(2, 4, figsize=(15, 7.5), constrained_layout=True)
            for ax, org in zip(axes.flat, orgs):
                panel(ax, inp, org, metric, org)
            for ax in axes.flat[len(orgs):]:
                ax.axis("off")
            for ax in axes[:, 0]:
                ax.set_ylabel("contrast alpha", fontsize=8)
            for ax in axes[-1, :]:
                ax.set_xlabel("Kozak weight w", fontsize=8)
            fig.suptitle(f"{inp}: {label[metric]} (x100) over the (alpha, w) grid. Colour spans each "
                         f"organism's best down to 5 points below; bold = best cell; red box = default "
                         f"(0.4, 25)", fontsize=9, color=INK)
            for ext in ("png", "pdf"):
                fig.savefig(f"{a.out}/heatmap_{metric}_{inp}.{ext}", dpi=150)
            plt.close(fig)

    # summary: mean F1, max regret, organisms within 0.5 pt -- per input
    fig, axes = plt.subplots(3, 3, figsize=(14, 11), constrained_layout=True)
    for row, inp in enumerate(INPUTS + ["both"]):
        if not full[inp]:
            continue
        specs = [("mean_f1", "mean exact F1 over organisms (x100)", lambda v: f"{100 * v:.1f}", False),
                 ("max_regret_points", "worst-organism shortfall from its own best (F1 points)",
                  lambda v: f"{v:.1f}", True),
                 ("n_orgs_within_0p5pt", "organisms within 0.5 F1 points of their best",
                  lambda v: f"{v:d}", False)]
        for col, (key, ttl, fm, rev) in enumerate(specs):
            ax = axes[row, col]
            M = [[plat[(inp, (al, w))][key] for w in weights] for al in alphas]
            flat = [v for r in M for v in r]
            if key == "max_regret_points":
                vmin, vmax = 0, min(max(flat), 5)
                C = [[vmax - min(v, vmax) for v in r] for r in M]   # dark = small shortfall
                vmin, vmax = 0, vmax
            elif key == "mean_f1":
                vmax = max(flat); vmin = vmax - 0.05
                C = [[max(v, vmin) for v in r] for r in M]
            else:
                vmin, vmax = 0, len(full[inp])
                C = M
            ax.imshow(C, cmap=SEQ, vmin=vmin, vmax=vmax, aspect="auto", origin="lower")
            for i in range(len(alphas)):
                for j in range(len(weights)):
                    frac = (C[i][j] - vmin) / (vmax - vmin + 1e-12)
                    ax.text(j, i, fm(M[i][j]), ha="center", va="center", fontsize=6.5,
                            color="white" if frac > 0.55 else INK)
            di, dj = alphas.index(DEFAULT[0]), weights.index(DEFAULT[1])
            ax.add_patch(Rectangle((dj - 0.5, di - 0.5), 1, 1, fill=False, ec="#e34948", lw=2))
            ax.set_xticks(range(len(weights)), [f"{w:g}" for w in weights], fontsize=7)
            ax.set_yticks(range(len(alphas)), [f"{al:g}" for al in alphas], fontsize=7)
            ax.set_title(f"{inp}: {ttl}", fontsize=8, color=INK)
            ax.set_xlabel("Kozak weight w", fontsize=8)
            ax.set_ylabel("contrast alpha", fontsize=8)
            for s in ax.spines.values():
                s.set_visible(False)
    fig.suptitle("Across 8 organisms (darker = better; red box = default 0.4, 25; "
                 "'both' = per-organism mean of complete and fragments F1)", fontsize=9, color=INK)
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}/heatmap_summary.{ext}", dpi=150)
    plt.close(fig)

    # ---- GM12878
    gdir = f"{a.grid}/gm12878"
    if os.path.isdir(gdir):
        keys = ["MANE_assembled", "MANE_exact", "MANE_sens_exact_given_assembled", "MANE_compatible",
                "GENCODE_complete_exact", "GENCODE_complete_sens_exact_given_assembled",
                "distinct_predicted_complete", "distinct_predicted_5prime_partial",
                "pred_complete_MANE_precision_exact", "pred_all_GENCODE_exact_frac"]
        grows = []
        for asm in sorted(os.listdir(gdir)):
            for cfg in sorted(os.listdir(f"{gdir}/{asm}")):
                sf = f"{gdir}/{asm}/{cfg}/eval.summary.tsv"
                if not os.path.exists(sf):
                    continue
                kv = dict(line.rstrip("\n").split("\t") for line in open(sf))
                al, w = parse_cfg(cfg)
                r = dict(assembly=asm, config=cfg, alpha=al, weight=w)
                r.update({k: kv[k] for k in keys})
                # prediction side (eval.pred_status.tsv): complete predictions that are an in-frame
                # piece of a GENCODE CDS, split into exact (right start) and not exact (wrong start,
                # e.g. a truncated ORF trimmed to an internal ATG); partial predictions kept compatible
                pc = collections.Counter()
                with open(f"{gdir}/{asm}/{cfg}/eval.pred_status.tsv") as fh:
                    next(fh)
                    for line in fh:
                        x = line.rstrip("\n").split("\t")
                        if x[4] != "1":
                            continue
                        if x[1] == "complete":
                            pc["exact" if x[3] == "1" else "not_exact"] += 1
                        else:
                            pc["partial"] += 1
                r.update(pred_complete_gencode_exact=pc["exact"], pred_complete_compatible_not_exact=pc["not_exact"],
                         pred_partial_gencode_compatible=pc["partial"],
                         pred_complete_start_precision=f"{pc['exact'] / max(1, pc['exact'] + pc['not_exact']):.4f}")
                ref = f"{P}/work/gm12878_td/eval/{asm}_final_abinitio.summary.tsv"
                if cfg == dname and os.path.exists(ref):
                    rk = dict(line.rstrip("\n").split("\t") for line in open(ref))
                    r["matches_existing_final_abinitio"] = all(rk[k] == kv[k] for k in kv if k != "label")
                grows.append(r)
        if grows:
            cols = ["assembly", "config", "alpha", "weight"] + keys + [
                "pred_complete_gencode_exact", "pred_complete_compatible_not_exact",
                "pred_partial_gencode_compatible", "pred_complete_start_precision",
                "matches_existing_final_abinitio"]
            write_tsv(f"{a.out}/gm12878.tsv", grows, cols)
            asms = sorted({r["assembly"] for r in grows})
            gal = sorted({r["alpha"] for r in grows if r["alpha"] is not None})
            gw = sorted({r["weight"] for r in grows if r["weight"] is not None})
            gmet = [("MANE_sens_exact_given_assembled", "MANE exact | assembled"),
                    ("pred_complete_start_precision",
                     "exact start among GENCODE-\ncompatible complete calls")]
            fig, axes = plt.subplots(len(gmet), len(asms), figsize=(5.2 * len(asms), 3.9 * len(gmet)),
                                     constrained_layout=True, squeeze=False)
            for row, (key, ttl) in enumerate(gmet):
                for ax, asm in zip(axes[row], asms):
                    d = {(r["alpha"], r["weight"]): float(r[key])
                         for r in grows if r["assembly"] == asm and r["alpha"] is not None}
                    nr = [float(r[key]) for r in grows if r["assembly"] == asm and r["config"] == "norefine"]
                    M = [[d.get((al, w), float("nan")) for w in gw] for al in gal]
                    vmax = max(v for r in M for v in r); vmin = vmax - 0.03
                    ax.imshow([[max(v, vmin) for v in r] for r in M], cmap=SEQ, vmin=vmin, vmax=vmax,
                              aspect="auto", origin="lower")
                    for i in range(len(gal)):
                        for j in range(len(gw)):
                            v = M[i][j]
                            ax.text(j, i, f"{100 * v:.2f}", ha="center", va="center", fontsize=8,
                                    color="white" if (max(v, vmin) - vmin) / (vmax - vmin) > 0.55 else INK,
                                    fontweight="bold" if v == vmax else "normal")
                    if DEFAULT[0] in gal and DEFAULT[1] in gw:
                        di, dj = gal.index(DEFAULT[0]), gw.index(DEFAULT[1])
                        ax.add_patch(Rectangle((dj - 0.5, di - 0.5), 1, 1, fill=False, ec="#e34948", lw=2))
                    ax.set_xticks(range(len(gw)), [f"{w:g}" for w in gw], fontsize=8)
                    ax.set_yticks(range(len(gal)), [f"{al:g}" for al in gal], fontsize=8)
                    ax.set_xlabel("Kozak weight w", fontsize=8); ax.set_ylabel("contrast alpha", fontsize=8)
                    t = f"GM12878 {asm}: {ttl} (x100)"
                    if nr:
                        t += f"\nno refinement {100 * nr[0]:.2f}; colour spans best - 3 pts"
                    ax.set_title(t, fontsize=8, color=INK)
                    for s_ in ax.spines.values():
                        s_.set_visible(False)
            for ext in ("png", "pdf"):
                fig.savefig(f"{a.out}/gm12878_heatmap.{ext}", dpi=150)
            plt.close(fig)

    print(f"{len(rows)} grid rows; outputs in {a.out}")


if __name__ == "__main__":
    main()
