#!/usr/bin/env python
"""Q5: offline simulation of start-choice rules on the actual PyTD pre-refinement ORFs.

Re-implements predict._revise_start_sites (frame-contrast gate + candidate in-frame ATGs) on
best_candidates.gff3 of the pytd_dev_pwm run, then applies alternative choice rules and scores
exact-start accuracy among stop-correct predictions against the truth.
usage: simulate_rules.py <org> <variant>
"""
import sys, os, re, json, collections, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pwmlib import *

org, variant = sys.argv[1], sys.argv[2]
OUT = f"{ROOT}/work/analysis/pwm_audit/out"
rund = f"{RUNS}/{org}/{variant}/pytd_dev_pwm/input.fa.transdecoder_dir"
seqs = read_fasta(f"{INPUTS}/{org}/{variant}.fa")
truth = read_truth(f"{INPUTS}/{org}/{variant}.truth.tsv")
table = int(next(iter(truth.values()))["transl_table"])
stops = STOPS[table]
ALPHA = 0.4

# hexamer scores
scores = {}
for l in open(f"{rund}/hexamer.scores"):
    if l[0] == "#":
        continue
    x = l.rstrip().split("\t")
    if len(x) >= 4:
        scores[x[0]] = float(x[3])


def per_pos_scores(seq):
    n = len(seq)
    out = np.zeros(n)
    for i in range(n):
        k = min(i, 5)
        out[i] = scores.get(f"{seq[i-k:i+1]}-{i%3}", 0.0)
    return out


def cumulative_contrast(orf):
    """C[i] = frame_contrast(orf[:i]) for i multiple of 3 (i>=6), via cumsums."""
    c0 = np.concatenate([[0.0], np.cumsum(per_pos_scores(orf))])
    c1 = np.concatenate([[0.0], np.cumsum(per_pos_scores(orf[1:]))])
    c2 = np.concatenate([[0.0], np.cumsum(per_pos_scores(orf[2:]))])

    def fc(i):
        # _score_cds_via_markov returns 0 for sequences shorter than 6
        s0 = c0[i] if i >= 6 else 0.0
        s1 = c1[i - 1] if i - 1 >= 6 else 0.0
        s2 = c2[i - 2] if i - 2 >= 6 else 0.0
        return s0 - max(s1, s2)
    return fc


# coding contrast median from training ORFs (as in _refine_start_codons)
pb = []
for l in open(f"{rund}/__checkpoints_TDpredict/top_training_orfs.cds"):
    pass
from Bio import SeqIO
for rec in SeqIO.parse(f"{rund}/__checkpoints_TDpredict/top_training_orfs.cds", "fasta"):
    s = str(rec.seq).upper()
    if len(s) >= 30:
        fc = cumulative_contrast(s)
        pb.append(fc(len(s)) / len(s))
pb.sort()
mid = len(pb) // 2
coding_contrast = pb[mid] if len(pb) % 2 else (pb[mid - 1] + pb[mid]) / 2

# PWMs
Pp = load_pwm_file(f"{rund}/start_refinement.enhanced.+.pwm")
Pm = load_pwm_file(f"{rund}/start_refinement.-.pwm")
Praw = load_pwm_file(f"{rund}/start_refinement.+.pwm")
(up, down), thr, _, _ = select_range_threshold(
    f"{rund}/start_refinement.enhanced.feature.scores.roc.auc", f"{rund}/start_refinement.enhanced.feature.scores.roc")

# truth-trained PWMs, 2-fold by transcript hash (annotated 5' UTRs only for training)
def fold_of(tid):
    return int(hashlib.md5(tid.encode()).hexdigest(), 16) % 2
sets, meta = truth_sets(seqs, truth, annotated_only=True)
truthP = {}
for f in (0, 1):
    tr = 1 - f
    pos = [w for t, w in sets["pos"] if fold_of(t) == tr]
    neg = [w for k in ("in_frame_down", "in_frame_up", "utr_oof") for t, w in sets[k] if fold_of(t) == tr]
    truthP[f] = (pwm_from_feats(pos), pwm_from_feats(neg))
# truth PWM with TD-style negatives
truthP_td = {}
for f in (0, 1):
    tr = 1 - f
    pos = [w for t, w in sets["pos"] if fold_of(t) == tr]
    neg = [w for t, w in sets["down_any"] if fold_of(t) == tr]
    truthP_td[f] = (pwm_from_feats(pos), pwm_from_feats(neg))


def lo(w, P, N, u=20, d=10):
    X = encode([w])
    return float(score_range(X, P, N, u, d, norm=False)[0])


# ---- parse pre-refinement ORFs ----
orfs = []  # (tid, start, end, strand)
for l in open(f"{rund}/best_candidates.gff3"):
    x = l.rstrip("\n").split("\t")
    if l[0] == "#" or len(x) != 9 or x[2] != "mRNA":
        continue
    orfs.append((x[0], int(x[3]), int(x[4]), x[6]))

# ---- per ORF, compute candidate list with features ----
rules = ["first", "current_pwm", "self_raw_argmax", "self_raw_margin1", "self_raw_margin2", "self_raw_margin3",
         "self_enh_chosen_argmax", "truth_argmax", "truth_margin1", "truth_margin2", "truth_margin3", "truth_td_argmax",
         "truth_td_margin2", "oracle"]
results = {m: collections.Counter() for m in ("open", "closed_ext")}
details = collections.Counter()


def choose(cands, feats, rule, fold):
    """cands: list of i (orf offsets); feats: dict i -> window or None. returns chosen i."""
    first = cands[0]
    if rule == "first":
        return first
    if rule == "current_pwm":
        for i in cands:
            w = feats[i]
            if w is None:
                continue
            sc = lo(w, Pp, Pm, up, down) / (up + down)
            if sc >= thr:
                return i
        return first
    if rule.startswith("self_raw"):
        P, N, u, d = Praw, Pm, 20, 10
    elif rule.startswith("self_enh_chosen"):
        P, N, u, d = Pp, Pm, up, down
    elif rule.startswith("truth_td"):
        (P, N), u, d = truthP_td[fold], 20, 10
    elif rule.startswith("truth"):
        (P, N), u, d = truthP[fold], 20, 10
    else:
        raise ValueError(rule)
    sc = {i: (lo(feats[i], P, N, u, d) if feats[i] is not None else -1e9) for i in cands}
    m = re.search(r"margin(\d+)", rule)
    if m:
        margin = float(m.group(1))
        best = first
        for i in cands[1:]:
            if sc[i] > sc[best] + margin:
                best = i
        return best
    return max(cands, key=lambda i: sc[i])


for tid, start, end, strand in orfs:
    r = truth.get(tid)
    if r is None or r["strand"] != strand:
        continue
    ts, te = int(r["cds_start"]), int(r["cds_end"])
    if abs(end - te) > 2:
        continue  # stop-incorrect: start choice irrelevant for the exact metric
    s = seqs[tid]
    oriented = s if strand == "+" else rc(s)
    orf_start = start - 1 if strand == "+" else len(s) - end
    orf = oriented[orf_start:orf_start + (end - start + 1)]
    fc = cumulative_contrast(orf)
    closed = orf[:3] == "ATG"
    mode = "closed_ext" if closed else "open"
    fold = fold_of(tid)
    true_i = (ts - 1) - orf_start  # true start as ORF offset
    cands = [0] if closed else []
    for i in range(3, len(orf) - 3, 3):
        if orf[i:i + 3] != "ATG":
            continue
        c = fc(i)
        if c >= ALPHA * coding_contrast * i:
            break
        cands.append(i)
    if not cands:
        # open end kept as 5' partial: exact only if the truth start is at the transcript end (never with ATG)
        results[mode]["n"] += 1
        results[mode]["kept_partial"] += 1
        for rule in rules:
            if abs(true_i - 0) <= 2:
                results[mode][rule] += 1
        continue
    results[mode]["n"] += 1
    results[mode]["ncand_total"] += len(cands)
    if len(cands) > 1:
        results[mode]["multi_cand"] += 1
    truth_in = true_i in cands
    results[mode]["truth_in_cands"] += truth_in
    if truth_in:
        results[mode]["truth_is_first"] += (cands[0] == true_i)
    feats = {i: window(oriented, orf_start + i) for i in cands}
    for rule in rules:
        if rule == "oracle":
            ch = true_i if truth_in else cands[0]
        else:
            ch = choose(cands, feats, rule, fold)
        if abs(ch - true_i) <= 2:
            results[mode][rule] += 1
        if rule == "current_pwm" and mode == "open":
            details["current_moved"] += (ch != cands[0])

out = {"org": org, "variant": variant, "coding_contrast": coding_contrast, "enh_range": [up, down], "enh_thr": thr,
       "stop_correct_orfs": sum(results[m]["n"] for m in results), "rules": rules}
for m in results:
    out[m] = dict(results[m])
out["note"] = ("'open' = ORFs with an open 5' end (what PyTD refines now, closed ORFs excluded). "
               "'closed_ext' = ORFs starting with ATG, where the same candidate loop is applied from the first ATG "
               "(offset 0 is always a candidate), a possible extension to handle N-terminal over-extension.")
json.dump(out, open(f"{OUT}/{org}.{variant}.sim.json", "w"), indent=1)
print(json.dumps(out, indent=1))
