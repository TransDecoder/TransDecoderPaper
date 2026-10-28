#!/usr/bin/env python3
"""E33 step 2: score the self-trained decisions against the release rule (truth used only here).

usage: python 05_evaluate.py <features_dir> <selftrain_dir> <manifest.tsv> <out_dir>

manifest columns: dataset, run_dir, truth, excluded ('-' for none).
For each dataset and method, the run's final predictions (preds.tsv = refined best_candidates) are
re-derived with the method's keep/trim decision on every 5'-open candidate with an in-frame ATG
(stops never change). Methods: release (the run), every posterior column of
<dataset>.decisions.tsv.gz (keep when >= 0.5), and 'perfect' (oracle: keep truncated, trim UTR
extensions whose start is the first ATG, release decision otherwise; the ceiling for start decisions).
Scoring as score_predictions.py / 13_depth_collect.py: stop within 2 nt, exact = start also within
2 nt, one credit per truth CDS; predictions on excluded contigs (simulated assemblies) removed.
Per type: F1 = 2 tp_exact / (n_truth + predictions on transcripts of that type); overall exact F1 =
2 tp_exact / (n_truth + all scored predictions).
Writes decision_metrics.tsv (accuracy / AUC on labelled candidates), exact_f1.tsv (overall and per
type), pi_estimates.tsv.
"""
import collections
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

fdir, sdir, manifest, odir = sys.argv[1:5]
os.makedirs(odir, exist_ok=True)
TOL = 2
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]


def ends(s, e, strand):
    return (s, e) if strand == "+" else (e, s)


def read_preds(path, excl):
    preds = collections.defaultdict(list)
    with open(path) as fh:
        fh.readline()
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if x[0] in excl:
                continue
            preds[x[0]].append((int(x[1]), int(x[2]), x[3]))
    return preds


def score(truth, preds):
    acc = {t: [0, 0, 0] for t in TYPES}  # n_truth, n_pred_on_type, tp_exact
    for tid, (s, e, strand, typ) in truth.items():
        a = acc[typ]
        a[0] += 1
        pl = preds.get(tid, [])
        a[1] += len(pl)
        t5, t3 = ends(s, e, strand)
        ex = False
        for q in pl:
            if q[2] != strand:
                continue
            q5, q3 = ends(*q)
            if abs(q3 - t3) <= TOL and abs(q5 - t5) <= TOL:
                ex = True
        a[2] += ex
    n_pred = sum(len(v) for v in preds.values())
    return acc, n_pred


man = pd.read_csv(manifest, sep="\t")
dm_rows, f1_rows, pi_rows = [], [], []
for _, m in man.iterrows():
    ds = m.dataset
    F = pd.read_csv(f"{fdir}/{ds}.tsv.gz", sep="\t", low_memory=False)
    D = pd.read_csv(f"{sdir}/{ds}.decisions.tsv.gz", sep="\t")
    F = F.merge(D, on="orf_id", how="left", validate="one_to_one")
    methods = [c for c in D.columns if c != "orf_id"]
    F["perfect"] = np.where(F.label == "truncated", 1.0,
                            np.where((F.label == "utr") & (F.true_atg_rank == 1), 0.0, 1.0 - F.trim_run))
    F["release"] = 1.0 - F.trim_run
    allm = ["release", "perfect"] + methods
    # decision metrics on labelled candidates
    lab = F[(F.label == "truncated") | ((F.label == "utr") & (F.true_atg_rank >= 1))]
    y = (lab.label == "truncated").astype(int).values
    for meth in allm:
        p = lab[meth].values
        sc = (lab.coding - lab.utr_score).values if meth == "release" else p
        dm_rows.append(dict(dataset=ds, method=meth, n=len(y), pi_true=y.mean(),
                            acc=((p >= 0.5) == (y == 1)).mean(),
                            auc=roc_auc_score(y, sc) if 0 < y.sum() < len(y) else np.nan,
                            keep_rate=(F[meth] >= 0.5).mean()))
    # truth and predictions
    excl = set()
    if m.excluded != "-":
        excl = {l.split("\t")[0] for l in open(m.excluded)}
    truth = {}
    with open(m.truth) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            truth[r["transcript_id"]] = (int(r["cds_start"]), int(r["cds_end"]), r["strand"], r["type"])
    base = read_preds(os.path.join(m.run_dir, "preds.tsv"), excl)
    # coordinates of each candidate when kept open / trimmed
    L, o5, ol, sg = F.tx_len.values, F.orf5.values, F.orf_len.values, F.seg.values
    plus = (F.strand == "+").values
    open_c = np.where(plus, o5 + 1, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5)
    trim_c = np.where(plus, o5 + 1 + sg, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5 - sg)
    tids, strands, trim_run = F.transcript_id.values, F.strand.values, F.trim_run.values
    for meth in allm:
        keep = F[meth].values >= 0.5
        preds = {k: list(v) for k, v in base.items()}
        n_changed = 0
        for i in np.where(keep == (trim_run == 1))[0]:  # decision differs from the release
            tid = tids[i]
            if tid in excl:
                continue
            old = (int(trim_c[0][i]), int(trim_c[1][i])) if trim_run[i] else (int(open_c[0][i]), int(open_c[1][i]))
            new = (int(open_c[0][i]), int(open_c[1][i])) if keep[i] else (int(trim_c[0][i]), int(trim_c[1][i]))
            lst = preds[tid]
            j = lst.index((old[0], old[1], strands[i]))
            lst[j] = (new[0], new[1], strands[i])
            n_changed += 1
        acc, n_pred = score(truth, preds)
        n_truth = sum(a[0] for a in acc.values())
        tpe = sum(a[2] for a in acc.values())
        row = dict(dataset=ds, method=meth, n_changed=n_changed, n_truth=n_truth, n_pred=n_pred, tp_exact=tpe,
                   f1_exact=2 * tpe / (n_truth + n_pred))
        for t in TYPES:
            n, npr, tp = acc[t]
            row[f"n_{t}"], row[f"tp_{t}"] = n, tp
            row[f"f1_{t}"] = 2 * tp / (n + npr) if n + npr else np.nan
        f1_rows.append(row)
    lg = dict(l.rstrip("\n").split("\t") for l in open(f"{sdir}/{ds}.log.tsv"))
    lg["pi_true_labelled"] = f"{y.mean():.4f}"
    pi_rows.append(lg)
    print(ds, "done", file=sys.stderr)

pd.DataFrame(dm_rows).to_csv(f"{odir}/decision_metrics.tsv", sep="\t", index=False, float_format="%.4f")
pd.DataFrame(f1_rows).to_csv(f"{odir}/exact_f1.tsv", sep="\t", index=False, float_format="%.5f")
pd.DataFrame(pi_rows).to_csv(f"{odir}/selftrain_log.tsv", sep="\t", index=False)
