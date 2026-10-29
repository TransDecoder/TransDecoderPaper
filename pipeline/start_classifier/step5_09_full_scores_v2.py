#!/usr/bin/env python3
"""E33 step 5, full homology mode: score Swiss-Prot vs UniRef50 runs (and ab initio) on the full-mode datasets.

usage: python step5_09_full_scores_v2.py <out_dir>

v2: adds the post hoc arm `hybrid` (Swiss-Prot hits + UniRef50 hits for ORFs without a Swiss-Prot hit,
step5_08_full_predict_v2.sh, step5_full/hybrid/<dataset>).

Arms per dataset:
  abinitio      release ab initio run (bench pytd_v7; athal depth-series A; human d8 pytd_v7)
  sprot_noassist  Swiss-Prot hits for selection only (bench pytd_v7_homology; athal C); not for human d8
  sprot         Swiss-Prot hits, assist on (bench pytd_v7h_homology; athal B; human d8 step5_full/sprot)
  uniref50      UniRef50-minus-taxon hits + the same Pfam hits, assist on (step5_full/uniref50)
Scoring as 16_score_adoption.py / score_predictions.py (stop within 2 nt, exact = start also within 2 nt; F1 =
2 tp / (n_truth + n_pred); per type with negatives = 2 tp_exact / (n_truth of type + predictions on
transcripts of that type)); simulated assemblies: predictions on excluded contigs removed (as E30/E32).
Stop sensitivity = tp_stop / n_truth, stop precision = tp_stop / n_pred.
Negatives: benchmark rand_* (wRand) and lnc_* (wLnc) transcripts; simulated assemblies: labels negatives.txt
(lncRNA-only contigs). Reported: predictions on negatives, negatives with >= 1 prediction.
Also parses the Predict log's homology start-refinement counts and the Predict run time.
Writes full_scores.tsv.
"""
import os
import re
import sys

import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
sys.path.insert(0, f"{P}/TransDecoderPaper/benchmark_v2/scripts")
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
F5 = f"{P}/work/start_classifier/step5_full"
R = f"{P}/work/bench_v3/runs"
AD = f"{P}/work/start_classifier/athal_depth"


def score(truth, preds):
    n_pred = sum(len(v) for v in preds.values())
    row = dict(n_truth=len(truth), n_pred=n_pred, tp_stop=0, tp_exact=0)
    per = {t: [0, 0, 0] for t in TYPES}
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
        row["tp_stop"] += best is not None
        row["tp_exact"] += best == "exact"
        a = per[t[3]]
        a[0] += 1
        a[1] += len(preds.get(tid, []))
        a[2] += best == "exact"
    row["f1_stop"] = 2 * row["tp_stop"] / (row["n_truth"] + n_pred)
    row["f1_exact"] = 2 * row["tp_exact"] / (row["n_truth"] + n_pred)
    row["sens_stop"] = row["tp_stop"] / row["n_truth"]
    row["prec_stop"] = row["tp_stop"] / n_pred if n_pred else float("nan")
    row["sens_exact"] = row["tp_exact"] / row["n_truth"]
    for t, (n, npr, tp) in per.items():
        row[f"n_{t}"], row[f"tp_exact_{t}"] = n, tp
        row[f"f1_exact_{t}"] = 2 * tp / (n + npr) if n + npr else float("nan")
    return row


def log_counts(run):
    d = {}
    f = f"{run}/predict.log"
    if not os.path.exists(f):
        return d
    for line in open(f):
        m = re.search(r"R1 \(BLASTP, keep open\) fired (\d+), R2 \(BLASTP, trim\) (\d+), R3 \(Pfam, keep open\) (\d+); "
                      r"homology decided (\d+) trim / (\d+) keep open, (\d+) conflicts.*?; (\d+) decisions differ", line)
        if m:
            d.update(zip(["R1", "R2", "R3", "hom_trim", "hom_keep", "conflicts", "changed_vs_abinitio"],
                         map(int, m.groups())))
        m = re.search(r"BLASTP (\d+) ORFs with hits, Pfam (\d+) ORFs with domains", line)
        if m:
            d["orfs_blast_hit"], d["orfs_pfam"] = int(m.group(1)), int(m.group(2))
    rt = f"{run}/runtime_sec.txt"
    if os.path.exists(rt):
        d["predict_sec"] = float(open(rt).read().strip())
    return d


datasets = []
for o in ("human", "athal", "dmel", "tthermophila"):
    for v in ("complete", "fragments", "complete.wRand", "fragments.wRand") + \
            (("complete.wLnc", "fragments.wLnc") if o == "human" else ()):
        datasets.append((f"bench.{o}.{v}", dict(
            truth=f"{P}/data/bench_v3/inputs/{o}/{v}.truth.tsv", excluded=None, negatives=None,
            arms=dict(abinitio=f"{R}/{o}/{v}/pytd_v7", sprot_noassist=f"{R}/{o}/{v}/pytd_v7_homology",
                      sprot=f"{R}/{o}/{v}/pytd_v7h_homology", uniref50=f"{F5}/uniref50/bench.{o}.{v}",
                      hybrid=f"{F5}/hybrid/bench.{o}.{v}"))))
for lv in ("d1", "d16"):
    L = f"{P}/data/simtrinity_athal/depth/{lv}"
    datasets.append((f"sim.athal.{lv}", dict(
        truth=f"{L}/inputs/fragments_simtrinity.truth.tsv", excluded=f"{L}/labels/sim_trinity.excluded.txt",
        negatives=f"{L}/labels/sim_trinity.negatives.txt",
        arms=dict(abinitio=f"{AD}/A/depth/{lv}/runs/pytd_v7", sprot_noassist=f"{AD}/C/depth/{lv}/runs/pytd_v7",
                  sprot=f"{AD}/B/depth/{lv}/runs/pytd_v7", uniref50=f"{F5}/uniref50/sim.athal.{lv}",
                  hybrid=f"{F5}/hybrid/sim.athal.{lv}"))))
L = f"{P}/data/simtrinity/depth/d8"
datasets.append(("sim.human.d8", dict(
    truth=f"{L}/inputs/fragments_simtrinity.truth.tsv", excluded=f"{L}/labels/sim_trinity.excluded.txt",
    negatives=f"{L}/labels/sim_trinity.negatives.txt",
    arms=dict(abinitio=f"{P}/work/simtrinity/depth/d8/runs/pytd_v7", sprot=f"{F5}/sprot/sim.human.d8",
              uniref50=f"{F5}/uniref50/sim.human.d8",
              hybrid=f"{F5}/hybrid/sim.human.d8"))))

rows = []
for ds, d in datasets:
    truth = read_truth(d["truth"])
    excl = {l.split("\t")[0] for l in open(d["excluded"])} if d["excluded"] else set()
    negs = {l.strip() for l in open(d["negatives"])} if d["negatives"] else None
    if negs is None:
        fa = os.path.realpath(f"{d['arms']['abinitio']}/input.fa")
        n_neg_total = sum(1 for l in open(fa) if l.startswith((">rand_", ">lnc_")))
    else:
        n_neg_total = len(negs - excl)
    for arm, run in d["arms"].items():
        pf = f"{run}/preds.tsv"
        if not os.path.exists(pf):
            print("missing", ds, arm, file=sys.stderr)
            continue
        preds = {k: v for k, v in read_preds(pf).items() if k not in excl}
        if negs is None:
            neg = {k: v for k, v in preds.items() if k.startswith(("rand_", "lnc_"))}
        else:
            neg = {k: v for k, v in preds.items() if k in negs}
        rows.append(dict(dataset=ds, arm=arm, run=run, **score(truth, preds),
                         n_pred_on_negatives=sum(len(v) for v in neg.values()), n_negatives_with_pred=len(neg), n_negatives=n_neg_total,
                         **log_counts(run)))
        print(ds, arm, "scored", file=sys.stderr, flush=True)
pd.DataFrame(rows).to_csv(f"{out}/full_scores.tsv", sep="\t", index=False, float_format="%.5f")
