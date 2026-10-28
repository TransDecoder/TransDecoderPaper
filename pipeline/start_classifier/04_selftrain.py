#!/usr/bin/env python3
"""E33 step 2: self-trained keep-open / trim decisions for one input (no truth labels used).

usage: PYTHONPATH=<PyTD tree> python 04_selftrain.py <run_dir> <genetic_code> <features.tsv.gz> <dataset> <out_prefix>

Uses only the finished PyTD run (models and selected ORFs) and the input's transcripts. The label
columns of <features.tsv.gz> (01_features.py) are NOT read; only the candidates' features are.

(a) Synthetic self-training. Confident ORFs = the run's selected ORFs (best_candidates.gff3, i.e.
    before refinement) that begin with ATG (complete or 3' partial; the benchmark 'fragments' inputs
    have no complete CDSs) and are at least MIN_CONF_CODONS long, at most
    MAX_CONF (random subset, seeded). From each confident transcript:
    - truncated examples: N_CUT cuts at uniform positions inside the CDS (after the ATG, leaving at
      least the input's minimum candidate ORF length); the 5' piece is removed, so u = real CDS from the
      first whole codon to the first in-frame ATG downstream of the cut;
    - UTR examples: N_CUT cuts at uniform positions in the open 5' UTR segment (after the last in-frame
      stop upstream of the ATG, >= 3 nt before the ATG), so u = real UTR and the ATG is the true start.
    - with probability q (the share of the input's real candidates whose ORF lacks a stop) the 3' end
      is also cut inside the CDS (ORF becomes internal / 3' partial; no 3' UTR).
    If fewer than MIN_UTR UTR examples result, decoy UTR examples are added (as many as truncated
    examples / 2): a stop-free, ATG-free in-frame stretch of a PyTD dinucleotide shuffle of a confident
    transcript (seed '42:<tid>:<k>', core/decoys.py) of length |u| - 6 drawn from the real candidates'
    |u| < 200, followed by the confident ORF from 6 nt upstream of its ATG.
    Every synthetic candidate is scored with PyTD's own code exactly like a real one.
    Weights: 'nat' = 1 (classes balanced); 'match' = per class, reweighted so the joint histogram of
    (|u| decile, remaining transcript length quintile) equals that of the real candidates (clipped at 20;
    classes balanced).
    Models fitted on synthetic data (balanced, so training prior 0.5): logistic regression on
      u        |u|, log|u|, C/mu, C/mu/|u|, C/mu/sqrt|u|, K, Kozak bases off transcript, GC(u)
      uc       u + log 3' UTR length, ORF lacks stop, log ORF length
      uct      uc + log transcript length
    and gradient-boosted trees on uct (comparison only). The truncated share of the real candidates is
    estimated by EM prior adjustment (Saerens et al. 2002) and posteriors adjusted to it.
(b) Unsupervised mixture on the real candidates. Intrinsic log-odds L0 = the release rule score
    (C/mu - 0.4|u| - 17K) calibrated by logistic regression on the synthetic examples (b_rule*), or the
    (a) 'u' model's log-odds (b_au*). EM over the real candidates: pi and class-conditional histograms
    of context features (3' state = ORF lacks stop, else 3' UTR length quintile; ORF length quintile;
    for *_tx also transcript length quintile; +1 pseudo-counts), posterior logit = L0 + logit(pi) +
    sum of the features' log-likelihood ratios (naive Bayes). Initialised from the calibrated rule's
    (or 'u' model's) posteriors. b_rule_pi = EM over pi only (prior-adapted calibrated rule).
Decision: keep 5'-open when posterior >= 0.5.

Writes <out_prefix>.decisions.tsv.gz (one row per real candidate: orf_id and posterior per method),
<out_prefix>.synthetic.tsv.gz (synthetic examples), <out_prefix>.log.tsv (counts, pi estimates).
"""
import gzip
import math
import random
import re
import statistics
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from pytransdecoder.predict import TransDecoderPredict, DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT
from pytransdecoder.core import kozak as kozak_mod
from pytransdecoder.core import decoys
from pytransdecoder.core.prng import Xoshiro128ss
from pytransdecoder.core.sequence import reverse_complement

run, gcode, feat_file, dataset, outp = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
ALPHA, W = DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT
MIN_CONF_CODONS, MAX_CONF, N_CUT, MIN_UTR = 100, 4000, 6, 300
SEED = 20261028
rng = random.Random(f"E33s2:{dataset}")

tdp = TransDecoderPredict(transcripts_file=run / "input.fa", output_dir=run, genetic_code=gcode)
wd = tdp.workdir
hexfile = wd / "hexamer.scores"
scores = tdp._load_hexamer_scores(hexfile)
mu = tdp._coding_contrast(wd / "__checkpoints_TDpredict" / "top_training_orfs.cds", hexfile)
kfile = wd / "start_refinement.kozak.tsv"
kmodel = kozak_mod.read_model(kfile) if kfile.exists() else None
tx = tdp._load_transcripts()
stops = set(tdp.stop_codons)

# ---- real candidates (features only)
real = pd.read_csv(feat_file, sep="\t", usecols=["orf_id", "transcript_id", "tx_len", "orf_len", "seg", "coding",
                                                 "coding_pb", "kozak", "kozak_off", "utr_score", "internal",
                                                 "utr3", "seg_gc"])
q_int = real.internal.mean()
min_len = int(np.percentile(real.orf_len, 5))
min_len -= min_len % 3

# ---- confident ORFs
conf = []
for line in open(wd / "best_candidates.gff3"):
    p = line.rstrip("\n").split("\t")
    if line.startswith("#") or len(p) != 9 or p[2] != "mRNA":
        continue
    tid, s, e, strand = p[0], int(p[3]), int(p[4]), p[6]
    seq = tx[tid]
    if strand == "-":
        ori, o5 = reverse_complement(seq), len(seq) - e
    else:
        ori, o5 = seq, s - 1
    orf = ori[o5:o5 + e - s + 1]
    if orf[:3] == "ATG" and len(orf) >= 3 * MIN_CONF_CODONS:  # complete or 3' partial
        conf.append((tid, ori, o5, o5 + len(orf)))
rng.shuffle(conf)
conf = conf[:MAX_CONF]


def candidate(t, end, kind):
    """Score transcript t (ORF orientation) whose ORF begins at the first whole codon of the frame of
    'end' (exclusive ORF end) and runs to 'end'. Returns a feature dict or None."""
    o5 = end % 3
    orf = t[o5:end]
    atg = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"), None)
    if atg is None or len(orf) < min_len:
        return None
    u = orf[:atg]
    c = tdp._frame_contrast(u, scores) / mu
    k = kozak_mod.score(kmodel, t, o5 + atg) if kmodel else 0.0
    koff = sum(1 for b in kozak_mod.window(t, o5 + atg) if b is None)
    return dict(kind=kind, tx_len=len(t), orf_len=len(orf), seg=atg, coding=c, coding_pb=c / atg, kozak=k,
                kozak_off=koff, utr_score=ALPHA * atg + W * k, internal=int(orf[-3:] not in stops),
                utr3=len(t) - end, seg_gc=(u.count("G") + u.count("C")) / len(u))


def cut3(t, atg_abs, end):
    """With probability q_int cut the 3' end inside the CDS (in frame), else keep."""
    if rng.random() >= q_int:
        return t, end
    lo = atg_abs + 3 + min_len
    if lo >= end - 3:
        return t, end
    new_end = lo + 3 * rng.randrange((end - 3 - lo) // 3 + 1)
    return t[:new_end], new_end


syn = []
n_utr_src = 0
for tid, ori, s, e in conf:
    # truncated: cut inside the CDS
    hi = e - min_len
    if hi > s + 3:
        for _ in range(N_CUT):
            c = rng.randrange(s + 3, hi)
            t, end = ori[c:], e - c
            first_atg = next((i for i in range(end % 3 + 3, end - 3, 3) if t[i:i + 3] == "ATG"), None)
            if first_atg is None:
                continue
            t, end = cut3(t, first_atg, end)
            r = candidate(t, end, "truncated")
            if r:
                syn.append(r)
    # UTR: open 5' UTR after the last upstream in-frame stop
    p = s - 3
    while p >= 0 and ori[p:p + 3] not in stops:
        p -= 3
    p += 3  # first base after the stop (or the frame start)
    if s - p >= 3:
        n_utr_src += 1
        for _ in range(N_CUT):
            c = rng.randrange(p, s - 2)
            t, end = ori[c:], e - c
            t, end = cut3(t, s - c, end)
            r = candidate(t, end, "utr")
            if r and r["seg"] == s - c - (end % 3):
                syn.append(r)

n_trunc = sum(r["kind"] == "truncated" for r in syn)
n_utr_real = sum(r["kind"] == "utr" for r in syn)
n_decoy = 0
if n_utr_real < MIN_UTR:
    target = n_trunc // 2
    lens = [x for x in real.seg if 9 <= x < 200]
    tries = 0
    while n_decoy < target and tries < 20 * target and conf and lens:
        tries += 1
        tid, ori, s, e = conf[rng.randrange(len(conf))]
        L = rng.choice(lens) - 6
        if s < 6 or L < 3:
            continue
        ctx = ori[s - 6:s]
        if any(ctx[i:i + 3] in stops or ctx[i:i + 3] == "ATG" for i in (0, 3)):
            continue
        shuf = decoys.dinucleotide_shuffle(ori, Xoshiro128ss(f"42:{tid}:{tries}"))
        f = rng.randrange(3)
        cod = [shuf[i:i + 3] for i in range(f, len(shuf) - 2, 3)]
        ok, run_start = [], 0
        for j, cd in enumerate(cod + ["TAA"]):
            if cd in stops or cd == "ATG":
                if j - run_start >= L // 3:
                    ok.append((run_start, j))
                run_start = j + 1
        if not ok:
            continue
        a, b = ok[rng.randrange(len(ok))]
        st = a + rng.randrange(b - a - L // 3 + 1)
        dseq = "".join(cod[st:st + L // 3])
        t = dseq + ori[s - 6:]
        end = len(dseq) + 6 + (e - s)
        t, end = cut3(t, len(dseq) + 6, end)
        r = candidate(t, end, "utr")
        if r and r["seg"] == len(dseq) + 6:
            r["kind"] = "utr_decoy"
            syn.append(r)
            n_decoy += 1

S = pd.DataFrame(syn)
S["y"] = (S.kind == "truncated").astype(int)


def add_feats(d):
    d = d.copy()
    d["log_seg"] = np.log(d.seg)
    d["coding_z"] = d.coding / np.sqrt(d.seg)
    d["log_tx"] = np.log(d.tx_len)
    d["log_orf"] = np.log(d.orf_len)
    d["log_utr3"] = np.log1p(d.utr3)
    d["rule_score"] = d.coding - d.utr_score
    return d


S, R = add_feats(S), add_feats(real)

# weights
S["w_nat"] = np.where(S.y == 1, 0.5 / max((S.y == 1).sum(), 1), 0.5 / max((S.y == 0).sum(), 1)) * len(S)
segq = np.unique(np.quantile(R.seg, np.linspace(0, 1, 11)[1:-1]))
txq = np.unique(np.quantile(R.tx_len, np.linspace(0, 1, 6)[1:-1]))


def cell(d):
    return np.searchsorted(segq, d.seg.values, side="right") * 10 + np.searchsorted(txq, d.tx_len.values, side="right")


R_cells = pd.Series(cell(R)).value_counts(normalize=True)
S["cell"] = cell(S)
wm = np.zeros(len(S))
for k in (0, 1):
    m = (S.y == k).values
    sc = S.loc[m, "cell"].value_counts(normalize=True)
    w = S.loc[m, "cell"].map(lambda c: R_cells.get(c, 0.0) / sc[c]).clip(upper=20).values
    wm[m] = w / w.sum() * 0.5 * len(S)
S["w_match"] = wm

F_U = ["seg", "log_seg", "coding", "coding_pb", "coding_z", "kozak", "kozak_off", "seg_gc"]
F_UC = F_U + ["log_utr3", "internal", "log_orf"]
F_UCT = F_UC + ["log_tx"]


def logit(p):
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return np.log(p / (1 - p))


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def saerens(lo, n_iter=200):
    """EM prior adjustment from a training prior of 0.5; returns (pi, adjusted log-odds)."""
    pi = 0.5
    for _ in range(n_iter):
        new = sigmoid(lo + logit(np.array([pi]))[0]).mean()
        if abs(new - pi) < 1e-6:
            pi = new
            break
        pi = new
    return pi, lo + logit(np.array([pi]))[0]


def fit_lr(feats, wcol):
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000))
    m.fit(S[feats].values, S.y.values, logisticregression__sample_weight=S[wcol].values)
    return m


post, info = {}, {}
lo_u = None
for name, feats, wcol in [("a_u", F_U, "w_nat"), ("a_uc", F_UC, "w_nat"), ("a_uct", F_UCT, "w_nat"),
                          ("a_u_match", F_U, "w_match"), ("a_uc_match", F_UC, "w_match"),
                          ("a_uct_match", F_UCT, "w_match")]:
    m = fit_lr(feats, wcol)
    lo = logit(m.predict_proba(R[feats].values)[:, 1])
    if name == "a_u":
        lo_u = lo
        lo_u_syn = logit(m.predict_proba(S[feats].values)[:, 1])
    pi, adj = saerens(lo)
    post[name] = sigmoid(adj)
    post[name + "_noem"] = sigmoid(lo)
    info[f"pi_{name}"] = pi
hgb = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, l2_regularization=1.0,
                                     random_state=SEED)
hgb.fit(S[F_UCT].values, S.y.values, sample_weight=S.w_nat.values)
pi, adj = saerens(logit(hgb.predict_proba(R[F_UCT].values)[:, 1]))
post["a_hgb_uct"], info["pi_a_hgb_uct"] = sigmoid(adj), pi

# (b) mixture over real candidates
cal = LogisticRegression(max_iter=1000).fit(S[["rule_score"]].values, S.y.values, sample_weight=S.w_nat.values)
lo_rule = cal.decision_function(R[["rule_score"]].values)
info["rule_platt_slope"], info["rule_platt_icpt"] = cal.coef_[0][0], cal.intercept_[0]


def qbins(x, n=5):
    edges = np.unique(np.quantile(x, np.linspace(0, 1, n + 1)[1:-1]))
    return np.searchsorted(edges, x, side="right")


ctx = {"s3": np.where(R.internal.values == 1, -1, qbins(R.utr3.values)) + 1, "orf": qbins(R.orf_len.values),
       "tx": qbins(R.tx_len.values)}


def mixture(L0, feats, n_iter=300):
    p = sigmoid(L0)
    pi = p.mean()
    for _ in range(n_iter):
        llr = np.zeros(len(L0))
        for f in feats:
            b = ctx[f]
            nb = b.max() + 1
            hT = np.bincount(b, weights=p, minlength=nb) + 1
            hU = np.bincount(b, weights=1 - p, minlength=nb) + 1
            llr += np.log(hT / hT.sum())[b] - np.log(hU / hU.sum())[b]
        newp = sigmoid(L0 + logit(np.array([pi]))[0] + llr)
        pi_new = newp.mean()
        done = np.abs(newp - p).max() < 1e-6
        p, pi = newp, pi_new
        if done:
            break
    return pi, p


for name, L0, feats in [("b_rule_pi", lo_rule, []), ("b_rule_ctx", lo_rule, ["s3", "orf"]),
                        ("b_rule_ctx_tx", lo_rule, ["s3", "orf", "tx"]), ("b_au_ctx", lo_u, ["s3", "orf"]),
                        ("b_au_ctx_tx", lo_u, ["s3", "orf", "tx"])]:
    pi, p = mixture(L0, feats)
    post[name], info[f"pi_{name}"] = p, pi

out = pd.DataFrame({"orf_id": R.orf_id})
for k, v in post.items():
    out[k] = np.round(v, 5)
out.to_csv(f"{outp}.decisions.tsv.gz", sep="\t", index=False)
S.drop(columns=["cell"]).to_csv(f"{outp}.synthetic.tsv.gz", sep="\t", index=False, float_format="%.4f")
info.update(dataset=dataset, n_real=len(R), q_internal=round(q_int, 4), min_orf_len=min_len, n_conf=len(conf),
            n_conf_with_open_utr=n_utr_src, n_syn_truncated=n_trunc, n_syn_utr=n_utr_real, n_syn_utr_decoy=n_decoy,
            mu=mu)
with open(f"{outp}.log.tsv", "w") as o:
    for k, v in info.items():
        o.write(f"{k}\t{v:.4f}\n" if isinstance(v, float) else f"{k}\t{v}\n")
print(dataset, " ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}" for k, v in info.items()),
      file=sys.stderr)
