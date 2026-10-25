#!/usr/bin/env python3
"""Score the UTR-prior runs against the release (evidence E29).

usage: python 03_utr_eval.py <out_dir>

Benchmark (stranded complete, complete_fullutr, fragments, fragments_mixed_t1..t3; 8 organisms):
  release  = $P/work/bench_v3/runs/<org>/<variant>/pytd_v7/preds.tsv (complete*, fragments) or the
             lambda = 0 run of 03_utr_runs.sh (mixed inputs; identical code path to the release)
  utr_l1   = $P/work/utr_prior/runs/<org>/<variant>/utr_l1/preds.tsv
Scoring as benchmark_v2/scripts/score_predictions.py (whole input: exact/stop F1) and
pipeline/derive/02_partial_types.py (per reference type: sensitivity, precision).
GM12878: eval_genome_cds.py summaries of the lambda 0 / 1 runs (LRAA, StringTie).
Writes utr_prior_eval.bench.tsv, utr_prior_eval.bench_types.tsv, utr_prior_eval.gm12878.tsv.
"""
import collections
import importlib.util
import os
import sys
from pathlib import Path

P = Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER")
R = P / "TransDecoderPaper"
sys.path.insert(0, str(R / "benchmark_v2/scripts"))
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

spec = importlib.util.spec_from_file_location("pt", R / "pipeline/derive/02_partial_types.py")
pt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pt)

OUT = Path(sys.argv[1])
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii"]
VARS = ["complete", "complete_fullutr", "fragments", "fragments_mixed_t1", "fragments_mixed_t2", "fragments_mixed_t3"]


def overall(truth, preds):
    """tp_stop, tp_exact, n_pred (one credit per reference, as score_predictions.py)."""
    n_pred = sum(len(v) for v in preds.values())
    tps = tpe = 0
    for tid, t in truth.items():
        best = None
        for p in preds.get(tid, []):
            if p[2] != t[2]:
                continue
            p5, p3 = ends(*p)
            t5, t3 = ends(*t[:3])
            if abs(p3 - t3) <= TOL:
                lvl = "exact" if abs(p5 - t5) <= TOL else "stop"
                if best != "exact":
                    best = lvl
        tps += best is not None
        tpe += best == "exact"
    return tps, tpe, n_pred


rows, trows = [], []
for org in ORGS:
    for v in VARS:
        truth = read_truth(P / f"data/bench_v3/inputs/{org}/{v}.truth.tsv")
        rel = P / f"work/bench_v3/runs/{org}/{v}/pytd_v7/preds.tsv"
        if v.startswith("fragments_mixed"):
            rel = P / f"work/utr_prior/runs/{org}/{v}/utr_l0/preds.tsv"
        runs = {"release": rel, "utr_l1": P / f"work/utr_prior/runs/{org}/{v}/utr_l1/preds.tsv"}
        for name, pf in runs.items():
            if not pf.exists():
                continue
            preds = read_preds(pf)
            tps, tpe, npred = overall(truth, preds)
            n = len(truth)
            rows.append(dict(org=org, input=v, run=name, n_truth=n, n_pred=npred, tp_stop=tps, tp_exact=tpe,
                             f1_stop=round(2 * tps / (n + npred), 4), f1_exact=round(2 * tpe / (n + npred), 4)))
            for typ, k in sorted(pt.score(truth, preds).items()):
                trows.append(dict(org=org, input=v, run=name, type=typ, n_truth=k["n_truth"], n_pred=k["n_pred"],
                                  tp_stop=k["tp_stop"], tp_exact=k["tp_exact"],
                                  sens_exact=round(k["tp_exact"] / k["n_truth"], 4),
                                  prec_exact=round(k["tp_exact"] / k["n_pred"], 4) if k["n_pred"] else ""))


def write(path, rs):
    with open(path, "w") as o:
        o.write("\t".join(rs[0]) + "\n")
        for r in rs:
            o.write("\t".join(str(r[c]) for c in rs[0]) + "\n")


write(OUT / "utr_prior_eval.bench.tsv", rows)
write(OUT / "utr_prior_eval.bench_types.tsv", trows)

KEYS = ["distinct_predicted_cds", "distinct_predicted_complete", "distinct_predicted_5prime_partial", "MANE_assembled",
        "MANE_exact", "MANE_sens_exact_given_assembled", "GENCODE_complete_assembled", "GENCODE_complete_exact",
        "GENCODE_complete_sens_exact_given_assembled", "pred_complete_GENCODE_exact_frac",
        "pred_all_GENCODE_compatible_frac"]
g = []
for asm in ["lraa", "stringtie"]:
    for lam in ["0", "1"]:
        f = P / f"work/utr_prior/gm12878_eval/{asm}_utr_l{lam}.summary.tsv"
        if not f.exists():
            continue
        d = dict(l.rstrip("\n").split("\t") for l in open(f))
        g.append(dict(assembly=asm, run="release" if lam == "0" else "utr_l1", **{k: d[k] for k in KEYS}))
if g:
    write(OUT / "utr_prior_eval.gm12878.tsv", g)
