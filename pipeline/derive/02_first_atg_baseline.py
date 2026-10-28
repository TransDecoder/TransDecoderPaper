#!/usr/bin/env python3
"""First-ATG longest-ORF baseline at the release (E07 addendum, 2026-10-28).

The baseline is built exactly as E07's `longest_orf_atg`
(benchmark_v2/scripts/longest_orf_with_start.py, run unchanged):
  - 5'-open LongOrfs ORFs are trimmed to their first in-frame ATG, and dropped if they have none;
  - the longest resulting ORF per transcript is kept.
Here it is built from the LongOrfs output of the release runs
(<runs>/<org>/<variant>/pytd_v7/input.fa.transdecoder_dir/longest_orfs.cds; LongOrfs minimum
50 aa), not from the older `pytd` runs (minimum 100 aa) that E07 used.

Predictions are written to <work>/<org>/<variant>/longest_orf_atg.preds.tsv; the run directories
are only read.

Scoring:
  - Counts per sequence are those of score_predictions.py, via bootstrap_ci.per_seq_counts and
    bootstrap_ci.metrics (stop within 2 nt; exact = start also within 2 nt; one credit per truth CDS).
  - Medians are over the 8 organisms, as in Table 1 (notebook 02).
  - Mixed inputs: median over organisms within each trial, then the mean over the 3 trials, as in
    mixed_medians.tsv.
  - Negatives: share of rand_* sequences with a call.
  - Bootstrap: the scheme of bootstrap_ci.py. Sequences of each input are resampled B times (seed 1)
    with the same resample for both tools, and the median over organisms is taken per replicate.
    95% percentile intervals are given for each tool and for baseline minus pytd_v7.

usage: 02_first_atg_baseline.py [P] [--B 1000]
writes $T/data/benchmark/final_v7/first_atg_baseline.tsv (long format)
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("P", nargs="?", default="/broad/hptmp/bhaas/TRANSDECODER_PAPER")
ap.add_argument("--B", type=int, default=1000)
ap.add_argument("--seed", type=int, default=1)
args = ap.parse_args()

P = Path(args.P)
T = P / "TransDecoderPaper"
SCRIPTS = T / "benchmark_v2/scripts"
sys.path.insert(0, str(SCRIPTS))
from bootstrap_ci import metrics, per_seq_counts, seq_ids  # noqa: E402
from score_predictions import read_preds, read_truth  # noqa: E402

INPUTS = P / "data/bench_v3/inputs"
RUNS = P / "work/bench_v3/runs"
WORK = P / "work/bench_v3/first_atg_v7"
OUT = T / "data/benchmark/final_v7/first_atg_baseline.tsv"
ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
VARIANTS = ["complete", "fragments", "complete.wRand", "fragments.wRand",
            "fragments_mixed_t1", "fragments_mixed_t2", "fragments_mixed_t3"]
TOOLS = ["pytd_v7", "longest_orf_atg_v7"]
METRICS = ["sensitivity_stop", "precision_stop", "f1_stop", "sensitivity_exact", "precision_exact", "f1_exact"]
rng = np.random.default_rng(args.seed)

rows = []          # (scope, variant, org_or_trial, tool, metric, estimate, lo, hi, diff, diff_lo, diff_hi)
reps = {}          # (variant, tool, metric) -> list over orgs of (B,) arrays
points = {}        # (variant, org, tool, metric) -> point estimate
for v in VARIANTS:
    for org in ORGS:
        run = RUNS / org / v / "pytd_v7"
        cds = run / "input.fa.transdecoder_dir/longest_orfs.cds"
        out = WORK / org / v / "longest_orf_atg.preds.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as fh:
            subprocess.run([sys.executable, str(SCRIPTS / "longest_orf_with_start.py"), str(cds)],
                           stdout=fh, check=True)
        ids = seq_ids(INPUTS / org / f"{v}.fa")
        truth = read_truth(str(INPUTS / org / f"{v}.truth.tsv"))
        preds = {"pytd_v7": read_preds(str(run / "preds.tsv")), "longest_orf_atg_v7": read_preds(str(out))}
        C = np.stack([per_seq_counts(ids, truth, preds[t]) for t in TOOLS], axis=1)  # (n, tools, 4)
        pt = metrics(C.sum(axis=0))
        rand = [i for i in ids if i.startswith("rand_")]
        for j, t in enumerate(TOOLS):
            for m in METRICS:
                points[(v, org, t, m)] = float(pt[m][j])
            if rand:
                called = sum(1 for i in rand if preds[t].get(i))
                points[(v, org, t, "neg_with_pred_frac")] = called / len(rand)
        n = len(ids)
        boot = {m: np.empty((args.B, len(TOOLS))) for m in METRICS}
        for b0 in range(0, args.B, 50):
            nb = min(50, args.B - b0)
            W = rng.multinomial(n, np.full(n, 1.0 / n), size=nb).astype(np.float32)
            mb = metrics(np.einsum("bn,ntk->btk", W, C.astype(np.float32)))
            for m in METRICS:
                boot[m][b0:b0 + nb] = mb[m]
        for j, t in enumerate(TOOLS):
            for m in METRICS:
                reps.setdefault((v, t, m), []).append(boot[m][:, j])
        print(f"done {v} {org} n={n}", file=sys.stderr)


def f(x):
    return f"{x:.4f}"


# per-input point estimates
for (v, org, t, m), x in sorted(points.items()):
    rows.append(("per_input", v, org, t, m, f(x), "", "", "", "", ""))

# medians over organisms (non-mixed variants), with bootstrap intervals
for v in VARIANTS:
    if v.startswith("fragments_mixed"):
        continue
    for m in METRICS + ["neg_with_pred_frac"]:
        est = {t: float(np.median([points[(v, o, t, m)] for o in ORGS if (v, o, t, m) in points]))
               for t in TOOLS if any((v, o, t, m) in points for o in ORGS)}
        if not est:
            continue
        if m in METRICS:
            med = {t: np.median(np.stack(reps[(v, t, m)], axis=1), axis=1) for t in TOOLS}
            d = med["longest_orf_atg_v7"] - med["pytd_v7"]
            for t in TOOLS:
                lo, hi = np.percentile(med[t], [2.5, 97.5])
                dd = (f(est[t] - est["pytd_v7"]), *map(f, np.percentile(d, [2.5, 97.5]))) \
                    if t != "pytd_v7" else ("", "", "")
                rows.append(("median_8_orgs", v, "all", t, m, f(est[t]), f(lo), f(hi), *dd))
        else:
            for t in TOOLS:
                rows.append(("median_8_orgs", v, "all", t, m, f(est[t]), "", "", "", "", ""))

# mixed inputs: median over organisms per trial, mean over trials
for m in METRICS:
    for t in TOOLS:
        per_trial = []
        for k in (1, 2, 3):
            v = f"fragments_mixed_t{k}"
            med = float(np.median([points[(v, o, t, m)] for o in ORGS]))
            per_trial.append(med)
            rows.append(("mixed_trial_median", v, f"t{k}", t, m, f(med), "", "", "", "", ""))
        rows.append(("mixed_mean_of_trial_medians", "fragments_mixed", "t1-t3", t, m,
                     f(float(np.mean(per_trial))), f(min(per_trial)), f(max(per_trial)), "", "", ""))

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("scope\tvariant\torg_or_trial\ttool\tmetric\testimate\tlo\thi\tdiff_vs_pytd_v7\tdiff_lo\tdiff_hi\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
print(f"wrote {OUT} ({len(rows)} rows)", file=sys.stderr)
