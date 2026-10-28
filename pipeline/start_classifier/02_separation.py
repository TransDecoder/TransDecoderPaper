#!/usr/bin/env python3
"""E33 step 1: how well can truncated 5'-open ORFs be told from UTR extensions (oracle labels)?

usage: python 02_separation.py <features_dir> <out_dir>

Reads <features_dir>/*.tsv.gz (01_features.py). Uses labelled candidates (label truncated / utr;
utr rows with true_atg_rank 0 are dropped and counted). y = 1 for truncated (correct decision: keep
5'-open), 0 for utr (correct decision: trim; exact only when true_atg_rank == 1).

Models (p = P(truncated)); decision keep-open when p >= 0.5 unless stated:
  rule          release rule, score = C(u)/mu - (0.4|u| + 17K) (keep when >= 0); fixed, not fitted
  rule_lr       logistic regression on |u|, C(u)/mu, K (best linear reweighting of the rule)
  lr_intr       logistic regression, intrinsic features of u and its ATG only
  lr_full       logistic regression, intrinsic + transcript/ORF context features
  hgb_intr      gradient-boosted trees (sklearn HistGradientBoosting), intrinsic
  hgb_full      gradient-boosted trees, intrinsic + context
Schemes:
  cv_bench      benchmark pooled (48 inputs), 5-fold GroupKFold by organism+transcript
  looo          benchmark leave-one-organism-out
  bench2sim     train on all benchmark inputs, test on each simulated assembly; also with the
                logit shifted to the test set's true truncated share (oracle prior)
  cv_sim        5-fold CV within each simulated assembly (grouped by contig)
Metrics: AUC; accuracy at the model's decision; best accuracy over thresholds on the test set
(oracle threshold, a ceiling); exact = correct decisions that give the exact ORF (truncated kept, or
utr trimmed whose true start is the first ATG).
"""
import glob
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

fdir, odir = sys.argv[1], sys.argv[2]
os.makedirs(odir, exist_ok=True)
SEED = 20261028

df = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in sorted(glob.glob(f"{fdir}/*.tsv.gz"))],
               ignore_index=True)
parts = df["dataset"].str.split(".", expand=True)
df["source"], df["org"], df["input"] = parts[0], parts[1], parts[2]
df["itype"] = np.where(df["source"] == "sim", "sim_" + df["org"],
                       df["input"].str.replace(r"_t\d$", "", regex=True))

# ---- candidate counts per dataset
rows = []
for ds, g in df.groupby("dataset", sort=False):
    lab = g[g.label != "nostop"]
    rows.append(dict(dataset=ds, open5_with_atg=len(g), stop_matched=len(lab),
                     truncated=(lab.label == "truncated").sum(),
                     utr_first_atg=((lab.label == "utr") & (lab.true_atg_rank == 1)).sum(),
                     utr_later_atg=((lab.label == "utr") & (lab.true_atg_rank >= 2)).sum(),
                     utr_inconsistent=((lab.label == "utr") & (lab.true_atg_rank == 0)).sum(),
                     trim_mismatch_vs_run=(g.trim != g.trim_run).sum()))
cnt = pd.DataFrame(rows)
cnt["pi_truncated"] = cnt.truncated / (cnt.truncated + cnt.utr_first_atg + cnt.utr_later_atg)
cnt.to_csv(f"{odir}/candidates_by_dataset.tsv", sep="\t", index=False, float_format="%.4f")

d = df[(df.label == "truncated") | ((df.label == "utr") & (df.true_atg_rank >= 1))].copy()
d["y"] = (d.label == "truncated").astype(int)
d["exact_if_trim"] = ((d.label == "utr") & (d.true_atg_rank == 1)).astype(int)
d["rule_score"] = d.coding - d.utr_score
d["log_seg"] = np.log(d.seg)
d["coding_z"] = d.coding / np.sqrt(d.seg)
d["log_tx"] = np.log(d.tx_len)
d["log_orf"] = np.log(d.orf_len)
d["log_utr3"] = np.log1p(d.utr3)
d["log_orf_rel"] = np.log(d.orf_len_rel)
d["grp"] = d.org + ":" + d.transcript_id.astype(str)
SEGBINS = [0, 20, 50, 100, 200, 500, 10 ** 9]
SEGLAB = ["0-19", "20-49", "50-99", "100-199", "200-499", ">=500"]
d["seg_bin"] = pd.cut(d.seg, SEGBINS, right=False, labels=SEGLAB)

F_RULE = ["seg", "coding", "kozak"]
F_INTR = ["seg", "log_seg", "coding", "coding_pb", "coding_z", "kozak", "kozak_off", "seg_gc"]
F_FULL = F_INTR + ["log_tx", "log_orf", "log_orf_rel", "internal", "log_utr3"]
MODELS = {"rule_lr": ("lr", F_RULE), "lr_intr": ("lr", F_INTR), "lr_full": ("lr", F_FULL),
          "hgb_intr": ("hgb", F_INTR), "hgb_full": ("hgb", F_FULL)}


def fit(kind, X, y):
    if kind == "lr":
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
    else:
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                           l2_regularization=1.0, random_state=SEED)
    return m.fit(X, y)


def logit(p):
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return np.log(p / (1 - p))


def metrics(y, score, keep, exact_if_trim):
    """score: higher = more truncated-like; keep: boolean decision keep-open."""
    y = np.asarray(y); score = np.asarray(score, float); keep = np.asarray(keep, bool)
    out = dict(n=len(y), pi=y.mean() if len(y) else np.nan)
    out["auc"] = roc_auc_score(y, score) if 0 < y.sum() < len(y) else np.nan
    out["acc"] = (keep == (y == 1)).mean() if len(y) else np.nan
    # oracle threshold accuracy
    if len(y):
        o = np.argsort(-score, kind="mergesort")
        ys, ss = y[o], score[o]
        tp = np.cumsum(ys); fp = np.cumsum(1 - ys)
        # keep the top k (k = 0..n); correct = tp_k + (negatives - fp_k)
        neg = len(y) - y.sum()
        corr = np.concatenate([[neg], tp + (neg - fp)])
        # only cut between distinct scores
        valid = np.concatenate([[True], np.r_[ss[1:] != ss[:-1], True]])
        out["acc_best"] = corr[valid].max() / len(y)
    out["exact"] = int(((y == 1) & keep).sum() + ((~keep) & (np.asarray(exact_if_trim) == 1)).sum())
    out["exact_max"] = int(y.sum() + np.asarray(exact_if_trim).sum())
    return out


def evaluate(test, scores, keeps, tag):
    res = []
    for name in scores:
        m = metrics(test.y.values, scores[name], keeps[name], test.exact_if_trim.values)
        res.append(dict(tag, model=name, bin="all", **m))
        for b in SEGLAB:
            sel = (test.seg_bin == b).values
            if sel.sum() == 0:
                continue
            m = metrics(test.y.values[sel], scores[name][sel], keeps[name][sel], test.exact_if_trim.values[sel])
            res.append(dict(tag, model=name, bin=b, **m))
    return res


def predict_all(train, test, prior_shift=0.0):
    scores, keeps = {"rule": test.rule_score.values}, {"rule": test.rule_score.values >= 0}
    for name, (kind, feats) in MODELS.items():
        m = fit(kind, train[feats].values, train.y.values)
        lg = logit(m.predict_proba(test[feats].values)[:, 1]) + prior_shift
        scores[name], keeps[name] = lg, lg >= 0
    return scores, keeps


results = []
bench = d[d.source == "bench"].reset_index(drop=True)
sims = sorted(d[d.source == "sim"].dataset.unique())

# cv_bench: out-of-fold predictions, then metrics per input type and pooled
oof = {k: np.zeros(len(bench)) for k in ["rule"] + list(MODELS)}
for tr, te in GroupKFold(n_splits=5).split(bench, bench.y, bench.grp):
    s, _ = predict_all(bench.iloc[tr], bench.iloc[te])
    for k in s:
        oof[k][te] = s[k]
for it in ["ALL"] + sorted(bench.itype.unique()):
    sel = np.ones(len(bench), bool) if it == "ALL" else (bench.itype == it).values
    results += evaluate(bench[sel], {k: v[sel] for k, v in oof.items()}, {k: v[sel] >= 0 for k, v in oof.items()},
                        dict(scheme="cv_bench", test=it))
bench["oof_lr_full"], bench["oof_hgb_full"], bench["oof_lr_intr"] = oof["lr_full"], oof["hgb_full"], oof["lr_intr"]

# looo
for org in sorted(bench.org.unique()):
    tr, te = bench[bench.org != org], bench[bench.org == org]
    s, k = predict_all(tr, te)
    results += evaluate(te, s, k, dict(scheme="looo", test=org))
    for it in sorted(te.itype.unique()):
        sel = (te.itype == it).values
        results += [r for r in evaluate(te[sel], {a: b[sel] for a, b in s.items()}, {a: b[sel] for a, b in k.items()},
                                        dict(scheme="looo_type", test=f"{org}:{it}")) if r["bin"] == "all"]

# bench2sim (plain and oracle prior shift)
pi_tr = bench.y.mean()
for ds in sims:
    te = d[d.dataset == ds]
    s, k = predict_all(bench, te)
    results += evaluate(te, s, k, dict(scheme="bench2sim", test=ds))
    shift = logit(np.array([te.y.mean()]))[0] - logit(np.array([pi_tr]))[0]
    s2, k2 = predict_all(bench, te, prior_shift=shift)
    s2.pop("rule"); k2.pop("rule")
    results += evaluate(te, s2, k2, dict(scheme="bench2sim_oracle_prior", test=ds))

# cv_sim
for ds in sims:
    te_all = d[d.dataset == ds].reset_index(drop=True)
    o2 = {k: np.zeros(len(te_all)) for k in ["rule"] + list(MODELS)}
    for tr, te in GroupKFold(n_splits=5).split(te_all, te_all.y, te_all.transcript_id):
        s, _ = predict_all(te_all.iloc[tr], te_all.iloc[te])
        for k in s:
            o2[k][te] = s[k]
    results += evaluate(te_all, o2, {k: v >= 0 for k, v in o2.items()}, dict(scheme="cv_sim", test=ds))

res = pd.DataFrame(results)
res.to_csv(f"{odir}/separation.tsv", sep="\t", index=False, float_format="%.4f")

# ---- logistic regression coefficients on the pooled benchmark (standardised features)
coef_rows = []
for name in ["rule_lr", "lr_intr", "lr_full"]:
    kind, feats = MODELS[name]
    for setname, data in [("bench", bench)] + [(ds, d[d.dataset == ds]) for ds in sims]:
        m = fit(kind, data[feats].values, data.y.values)
        lr = m[-1]
        for f, c in zip(feats, lr.coef_[0]):
            coef_rows.append(dict(model=name, data=setname, feature=f, coef_std=c))
        coef_rows.append(dict(model=name, data=setname, feature="(intercept)", coef_std=lr.intercept_[0]))
pd.DataFrame(coef_rows).to_csv(f"{odir}/lr_coefficients.tsv", sep="\t", index=False, float_format="%.4f")

# ---- coding-score scaling with |u| by class, and the decoy null
grp_src = np.where(d.source == "bench", "bench", d.dataset)
d["grp_src"] = grp_src
sc = []
for (gs, b), g in d.groupby(["grp_src", "seg_bin"], observed=True):
    for cls, h in [("truncated", g[g.y == 1]), ("utr", g[g.y == 0])]:
        if len(h) < 5:
            continue
        dec = h.dec_coding.dropna()
        row = dict(data=gs, seg_bin=b, cls=cls, n=len(h), seg_median=h.seg.median(),
                   coding_median=h.coding.median(), coding_q25=h.coding.quantile(.25), coding_q75=h.coding.quantile(.75),
                   coding_mean=h.coding.mean(), coding_sd=h.coding.std(),
                   coding_pb_median=h.coding_pb.median(), kozak_median=h.kozak.median(),
                   decoy_n=len(dec), decoy_median=dec.median() if len(dec) else np.nan,
                   decoy_q25=dec.quantile(.25) if len(dec) else np.nan, decoy_q75=dec.quantile(.75) if len(dec) else np.nan,
                   decoy_mean=dec.mean() if len(dec) else np.nan, decoy_sd=dec.std() if len(dec) else np.nan)
        if len(dec) >= 5:
            ks = ks_2samp(h.coding.values, dec.values)
            row.update(ks_vs_decoy=ks.statistic, ks_p=ks.pvalue)
        sc.append(row)
pd.DataFrame(sc).to_csv(f"{odir}/coding_score_by_seg.tsv", sep="\t", index=False, float_format="%.4g")

# per-class fits: mean coding = a*|u|, sd = b*|u|^c (log-log), pooled over bins with >= 30 rows
fits = []
for gs, g in d.groupby("grp_src"):
    for cls, h in [("truncated", g[g.y == 1]), ("utr", g[g.y == 0]), ("decoy", None)]:
        if cls == "decoy":
            h = g[g.y == 0].dropna(subset=["dec_coding"]).assign(coding=lambda x: x.dec_coding)
        if len(h) < 100:
            continue
        a = (h.coding * h.seg).sum() / (h.seg ** 2).sum()
        bins = pd.qcut(h.seg, 10, duplicates="drop")
        agg = h.groupby(bins, observed=True).agg(s=("seg", "median"), sd=("coding", "std"), m=("coding", "mean"))
        agg = agg[(agg.sd > 0)]
        c, logb = np.polyfit(np.log(agg.s), np.log(agg.sd), 1)
        fits.append(dict(data=gs, cls=cls, n=len(h), mean_per_base=a, sd_exponent=c, sd_coef=math.exp(logb)))
pd.DataFrame(fits).to_csv(f"{odir}/coding_score_scaling_fit.tsv", sep="\t", index=False, float_format="%.4f")

# ---- compact summary to stdout
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
s = res[res.bin == "all"].pivot_table(index=["scheme", "test"], columns="model", values="auc")
print("AUC\n", s.round(3))
s = res[res.bin == "all"].pivot_table(index=["scheme", "test"], columns="model", values="acc")
print("ACC\n", s.round(3))
s = res[res.bin == "all"].pivot_table(index=["scheme", "test"], columns="model", values="exact")
print("EXACT\n", s)
