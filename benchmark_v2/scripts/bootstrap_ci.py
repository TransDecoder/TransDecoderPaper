#!/usr/bin/env python3
"""Bootstrap confidence intervals for the benchmark metrics.

For each organism and input, sequences of the input FASTA (coding transcripts and any random or
lncRNA negatives) are resampled with replacement, B times, and the same resample is applied to
every tool, so differences between tools are paired. Per-sequence counts (truth CDSs, predictions,
stop-level and exact TPs) are those of score_predictions.py, with which the point estimates agree
exactly. Medians over organisms are computed per replicate (organisms resampled independently).

Writes two TSVs:
  <out>.per_input.tsv  org, variant, tool, metric, estimate, lo, hi, and for each tool versus the
                       reference tool: diff (tool - reference), diff_lo, diff_hi
  <out>.median.tsv     the same for the median over organisms

usage: bootstrap_ci.py <inputs_dir> <runs_dir> <out_prefix> [--B 1000] [--variants complete,...]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

TOOLS = {  # label -> (run subdirectory, predictions file)
    "pytd_dev_v2": ("pytd_dev_v2", "preds.tsv"),
    "perl_td_v600": ("perl_td_v600", "preds.tsv"),
    "longest_orf": ("pytd", "longest_orf.preds.tsv"),
    "longest_orf_atg": ("pytd", "longest_orf_atg.preds.tsv"),
    "prodigal": ("prodigal", "preds.tsv"),
    "gmst": ("gmst", "preds.tsv"),
    "codan": ("codan", "preds.tsv"),
    "cpat": ("cpat", "preds.tsv"),
    "cpc2": ("cpc2", "preds.tsv"),
    "rnasamba": ("rnasamba", "preds.tsv"),
    "borf": ("borf", "preds.tsv"),
}
REF = "pytd_dev_v2"
ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
METRICS = ["f1_stop", "f1_exact", "precision_stop", "sensitivity_stop", "sensitivity_exact"]


def seq_ids(fasta):
    with open(fasta) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]


def per_seq_counts(ids, truth, preds):
    """Array (n_seq, 4): n_truth, n_pred, tp_stop, tp_exact, matching score_predictions.py."""
    a = np.zeros((len(ids), 4), dtype=np.int32)
    for i, tid in enumerate(ids):
        t = truth.get(tid)
        a[i, 0] = 1 if t else 0
        plist = preds.get(tid, [])
        a[i, 1] = len(plist)
        best = None
        for p in plist:
            if t and p[2] == t[2]:
                p5, p3 = ends(*p)
                t5, t3 = ends(*t[:3])
                if abs(p3 - t3) <= TOL:
                    best = "exact" if abs(p5 - t5) <= TOL or best == "exact" else (best or "stop")
        if best:
            a[i, 2] = 1
            a[i, 3] = 1 if best == "exact" else 0
    return a


def metrics(tot):
    """tot: (..., 4) summed counts -> dict of metric arrays."""
    n_truth, n_pred, tps, tpe = (tot[..., k].astype(float) for k in range(4))
    with np.errstate(divide="ignore", invalid="ignore"):
        out = {}
        for lvl, tp in (("stop", tps), ("exact", tpe)):
            sens = np.where(n_truth > 0, tp / n_truth, 0.0)
            prec = np.where(n_pred > 0, tp / n_pred, 0.0)
            f1 = np.where(sens + prec > 0, 2 * sens * prec / (sens + prec), 0.0)
            out[f"sensitivity_{lvl}"], out[f"precision_{lvl}"], out[f"f1_{lvl}"] = sens, prec, f1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs")
    ap.add_argument("runs")
    ap.add_argument("out")
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--variants", default="complete,fragments,complete.wRand,fragments.wRand")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    variants = args.variants.split(",")
    tools = list(TOOLS)
    # reps[variant][tool][metric] -> list over orgs of arrays (B,) (NaN where the run is missing)
    reps = {v: {t: {m: [] for m in METRICS} for t in tools} for v in variants}
    rows = []
    for v in variants:
        for org in ORGS:
            fa = os.path.join(args.inputs, org, v + ".fa")
            truth = read_truth(os.path.join(args.inputs, org, v + ".truth.tsv"))
            ids = seq_ids(fa)
            n = len(ids)
            counts, have = [], []
            for t in tools:
                sub, fn = TOOLS[t]
                pf = os.path.join(args.runs, org, v, sub, fn)
                ok = os.path.exists(pf) and os.path.getsize(pf) > 0
                have.append(ok)
                counts.append(per_seq_counts(ids, truth, read_preds(pf)) if ok else np.zeros((n, 4), np.int32))
            C = np.stack(counts, axis=1)  # (n, tools, 4)
            point = metrics(C.sum(axis=0))
            boot = {m: np.empty((args.B, len(tools))) for m in METRICS}
            for b0 in range(0, args.B, 50):
                nb = min(50, args.B - b0)
                W = rng.multinomial(n, np.full(n, 1.0 / n), size=nb).astype(np.float32)  # (nb, n)
                tot = np.einsum("bn,ntk->btk", W, C.astype(np.float32))
                mb = metrics(tot)
                for m in METRICS:
                    boot[m][b0:b0 + nb] = mb[m]
            ri = tools.index(REF)
            for j, t in enumerate(tools):
                for m in METRICS:
                    vals = boot[m][:, j] if have[j] else np.full(args.B, np.nan)
                    reps[v][t][m].append(vals)
                    if not have[j]:
                        continue
                    d = boot[m][:, j] - boot[m][:, ri]
                    rows.append((org, v, t, m, point[m][j], *np.percentile(vals, [2.5, 97.5]),
                                 point[m][j] - point[m][ri], *np.percentile(d, [2.5, 97.5])))
            print(f"done {org} {v} n={n}", file=sys.stderr)
    hdr = "org\tvariant\ttool\tmetric\testimate\tlo\thi\tdiff_vs_ref\tdiff_lo\tdiff_hi\n"
    with open(args.out + ".per_input.tsv", "w") as fh:
        fh.write(hdr)
        for r in rows:
            fh.write("\t".join(r[:4]) + "\t" + "\t".join(f"{x:.4f}" for x in r[4:]) + "\n")
    with open(args.out + ".median.tsv", "w") as fh:
        fh.write(hdr.replace("org\t", "n_orgs\t"))
        for v in variants:
            for t in tools:
                for m in METRICS:
                    A = np.stack(reps[v][t][m], axis=1)  # (B, orgs)
                    R = np.stack(reps[v][REF][m], axis=1)
                    ok = ~np.isnan(A[0])
                    if not ok.any():
                        continue
                    med = np.median(A[:, ok], axis=1)
                    dmed = med - np.median(R[:, ok], axis=1)
                    # point estimate: median of the point estimates, taken from per_input rows
                    pts = [r[4] for r in rows if r[1] == v and r[2] == t and r[3] == m]
                    rpts = [r[4] for r in rows if r[1] == v and r[2] == REF and r[3] == m
                            and r[0] in [x[0] for x in rows if x[1] == v and x[2] == t and x[3] == m]]
                    est = float(np.median(pts))
                    fh.write(f"{int(ok.sum())}\t{v}\t{t}\t{m}\t{est:.4f}\t"
                             + "\t".join(f"{x:.4f}" for x in np.percentile(med, [2.5, 97.5]))
                             + f"\t{est - float(np.median(rpts)):.4f}\t"
                             + "\t".join(f"{x:.4f}" for x in np.percentile(dmed, [2.5, 97.5])) + "\n")


if __name__ == "__main__":
    sys.exit(main())
