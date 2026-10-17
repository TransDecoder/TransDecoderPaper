#!/usr/bin/env python3
"""Compare start-codon refinement strategies on benchmark truth sets.

Takes a finished PyTransDecoder run made with --no-refine-starts (its
input.fa.transdecoder_dir holds best_candidates.gff3 and the training ORFs), re-trains
the start PWM under a given strategy, revises the starts of non-ATG (5'-open) ORFs,
and scores the result against the truth set.

Strategies
  legacy        PyTransDecoder's current procedure (random 30% noise-depletion seed,
                5 random 75/25 rounds scoring <=1000 test features, 10-threshold ROC,
                best-AUC window, best-F1 threshold); --seed sets its RNG.
  kfold         deterministic noise filter + k-fold CV over all features, exact
                (Mann-Whitney) AUC, threshold chosen on pooled out-of-fold scores
                (--thr-objective); --seed only sets the fold assignment.
  fixed         as kfold but the window is fixed at the full 20 up / 10 down.
  none          no refinement.

  --rule best   choose the highest-scoring candidate ATG at/above threshold (current)
  --rule first  choose the most upstream candidate ATG at/above threshold
  --gate coding only move a start if the sequence between the ORF's 5' end and the candidate
                ATG scores <= 0 under the ORF's own Markov coding model (i.e. looks like UTR,
                not like truncated coding sequence)

Output: one TSV row with the strategy, the selected window and threshold, the number
of revised starts, and stop-level / exact-level TPs against the truth.
"""
import argparse
import os
import random
import shutil
import sys
import tempfile

import numpy as np
from scipy.stats import rankdata

BASES = "GATC"
ATG_POS, LEFT, RIGHT = 20, 20, 10
FLEN = LEFT + 3 + RIGHT


def encode(features):
    idx = {b: i for i, b in enumerate(BASES)}
    return np.array([[idx[c] for c in f] for f in features], dtype=np.int8)


def pwm_probs(X, pseudo=0.1):
    counts = np.zeros((FLEN, 4))
    for j in range(4):
        counts[:, j] = (X == j).sum(axis=0)
    return (counts + pseudo) / (counts.sum(axis=1, keepdims=True) + 4 * pseudo)


def contributions(X, lo):
    """Per-position log-odds contributions, shape (n, FLEN)."""
    return lo[np.arange(FLEN)[None, :], X]


def window_scores(C, up, down):
    a, b = ATG_POS - up, ATG_POS + 3 + down
    return C[:, a:b].sum(axis=1) / (up + down)


def mw_auc(pos, neg):
    r = rankdata(np.concatenate([pos, neg]))
    n1, n2 = len(pos), len(neg)
    return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2)


def choose_threshold(pos, neg, objective):
    cands = np.unique(np.quantile(np.concatenate([pos, neg]), np.linspace(0, 1, 201)))
    best, best_t = -1, None
    w = len(pos) / len(neg)  # balanced: weight negatives to the positive count
    for t in cands:
        tp = (pos >= t).sum(); fn = len(pos) - tp
        fp = (neg >= t).sum(); tn = len(neg) - fp
        if objective == "f1_balanced":
            v = 2 * tp / (2 * tp + w * fp + fn)
        elif objective == "f1_raw":
            v = 2 * tp / (2 * tp + fp + fn)
        else:  # youden
            v = tp / len(pos) - fp / len(neg)
        if v > best:
            best, best_t = v, t
    return float(best_t)


def train_stable(P, N, strategy, seed, k, objective):
    # deterministic noise filter: drop positives that do not score above 0 under the
    # full-window log-odds of PWM(all positives) vs PWM(negatives)
    lo_all = np.log(pwm_probs(P) / pwm_probs(N))
    keep = contributions(P, lo_all).sum(axis=1) > 0
    P = P[keep]
    rng = np.random.default_rng(seed)
    fp, fn = rng.permutation(len(P)) % k, rng.permutation(len(N)) % k
    windows = [(u, d) for u in range(1, LEFT + 1) for d in range(1, RIGHT + 1)]
    if strategy == "fixed":
        windows = [(LEFT, RIGHT)]
    oof_p = {w: [] for w in windows}
    oof_n = {w: [] for w in windows}
    for f in range(k):
        lo = np.log(pwm_probs(P[fp != f]) / pwm_probs(N[fn != f]))
        Cp, Cn = contributions(P[fp == f], lo), contributions(N[fn == f], lo)
        for w in windows:
            oof_p[w].append(window_scores(Cp, *w))
            oof_n[w].append(window_scores(Cn, *w))
    aucs = {w: mw_auc(np.concatenate(oof_p[w]), np.concatenate(oof_n[w])) for w in windows}
    best_w = max(windows, key=lambda w: (aucs[w], -(w[0] + w[1])))
    thr = choose_threshold(np.concatenate(oof_p[best_w]), np.concatenate(oof_n[best_w]), objective)
    lo_final = np.log(pwm_probs(P) / pwm_probs(N))
    return lo_final, best_w, thr, aucs[best_w], int(keep.sum())


def revise(gff_in, transcripts, lo, window, thr, rule, revcomp, gate=None):
    """Return {(seqid, orf_id): (start, end, strand)} after start revision, and #revised."""
    out, n_rev = [], 0
    for line in open(gff_in):
        x = line.rstrip("\n").split("\t")
        if len(x) != 9 or x[2] != "mRNA":
            continue
        seqid, s, e, strand = x[0], int(x[3]), int(x[4]), x[6]
        seq = transcripts[seqid]
        o = revcomp(seq) if strand == "-" else seq
        sp = (len(o) - e + 1) if strand == "-" else s
        si = sp - 1
        if lo is None or o[si:si + 3] == "ATG":
            out.append((seqid, s, e, strand))
            continue
        orf_len = e - s + 1
        max_pos = max(si + 30, si + int(15 * orf_len / 100))
        best = None
        pos = o.find("ATG", 0)
        while pos != -1 and pos <= max_pos:
            if pos > si and (pos - si) % 3 == 0:
                a = pos - ATG_POS
                if a >= 0 and a + FLEN <= len(o):
                    feat = o[a:a + FLEN]
                    if set(feat) <= set(BASES):
                        sc = round(float(window_scores(contributions(encode([feat]), lo), *window)[0]), 3)
                        if sc >= thr and gate is not None and gate(o[si:pos]) > 0:
                            pass  # upstream segment looks coding: the ORF is likely truncated
                        elif sc >= thr and (best is None or (rule == "best" and sc > best[1])):
                            best = (pos, sc)
                            if rule == "first":
                                break
            pos = o.find("ATG", pos + 1)
        if best is None:
            out.append((seqid, s, e, strand))
            continue
        n_rev += 1
        ns = best[0] + 1
        if strand == "-":
            out.append((seqid, s, len(o) - ns + 1, strand))
        else:
            out.append((seqid, ns, e, strand))
    return out, n_rev


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="pytd_norefine run directory")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--pytd", required=True, help="PyTransDecoder source tree")
    ap.add_argument("--genetic-code", default="universal")
    ap.add_argument("--strategy", choices=["legacy", "kfold", "fixed", "none"], required=True)
    ap.add_argument("--rule", choices=["best", "first"], default="best")
    ap.add_argument("--gate", choices=["none", "coding"], default="none")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--thr-objective", choices=["f1_balanced", "f1_raw", "youden"], default="f1_balanced")
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    sys.path.insert(0, args.pytd)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../benchmark_v2/scripts")
    import pytransdecoder.predict as pred_mod
    from pytransdecoder.core.sequence import reverse_complement
    from score_predictions import read_truth, ends, TOL

    run_dir = os.path.abspath(args.run_dir)
    src_wd = os.path.join(run_dir, "input.fa.transdecoder_dir")
    gff_in = os.path.join(src_wd, "best_candidates.gff3")
    training = os.path.join(src_wd, "__checkpoints_TDpredict", "top_training_orfs.cds")

    tmp = tempfile.mkdtemp(prefix="startexp_")
    try:
        p = pred_mod.TransDecoderPredict(transcripts_file=os.path.join(run_dir, "input.fa"),
                                         output_dir=tmp, genetic_code=args.genetic_code)
        p.workdir.mkdir(parents=True, exist_ok=True)
        transcripts = p._load_transcripts()

        lo = window = thr = auc = None
        n_pos = n_kept = 0
        if args.strategy == "legacy":
            seed = args.seed
            orig = pred_mod.random.Random
            pred_mod.random.Random = lambda *_a, **_k: orig(seed)
            try:
                res = p._train_start_pwm(training)
            finally:
                pred_mod.random.Random = orig
            pwm_plus, pwm_minus, window, thr = res
            pp = np.array([[pwm_plus.pos_probs[i][b] for b in BASES] for i in range(FLEN)])
            pm = np.array([[pwm_minus.pos_probs[i][b] for b in BASES] for i in range(FLEN)])
            lo = np.log(pp / pm)
        elif args.strategy in ("kfold", "fixed"):
            posf, negf = p._extract_pwm_features(training, transcripts, LEFT, RIGHT)
            P, N = encode(posf), encode(negf)
            n_pos = len(P)
            lo, window, thr, auc, n_kept = train_stable(P, N, args.strategy, args.seed, args.k, args.thr_objective)

        gate = None
        if args.gate == "coding":
            hex_scores = {}
            for line in open(os.path.join(src_wd, "hexamer.scores")):
                if not line.startswith("#"):
                    x = line.rstrip("\n").split("\t")
                    if len(x) >= 4:
                        hex_scores[x[0]] = float(x[3])
            gate = lambda seg: p._score_cds_via_markov(seg, 5, hex_scores)
        preds, n_rev = revise(gff_in, transcripts, lo, window, thr if thr is not None else 0, args.rule,
                              reverse_complement, gate)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    truth = read_truth(args.truth)
    tp_stop = tp_exact = 0
    for seqid, s, e, strand in preds:
        t = truth.get(seqid)
        if not t or t[2] != strand:
            continue
        p5, p3 = ends(s, e, strand)
        t5, t3 = ends(*t[:3])
        if abs(p3 - t3) <= TOL:
            tp_stop += 1
            tp_exact += abs(p5 - t5) <= TOL
    cols = ["label", "strategy", "rule", "gate", "seed", "thr_objective", "window", "threshold", "cv_auc",
            "n_pos_features", "n_pos_kept", "n_preds", "n_revised", "tp_stop", "tp_exact", "n_truth"]
    vals = [args.label, args.strategy, args.rule, args.gate, args.seed, args.thr_objective,
            "NA" if window is None else f"{window[0]},{window[1]}",
            "NA" if thr is None else f"{thr:.3f}", "NA" if auc is None else f"{auc:.4f}",
            n_pos, n_kept, len(preds), n_rev, tp_stop, tp_exact, len(truth)]
    print("\t".join(cols))
    print("\t".join(map(str, vals)))


if __name__ == "__main__":
    sys.exit(main())
