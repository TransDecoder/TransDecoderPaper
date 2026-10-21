#!/usr/bin/env python
"""Q1-Q3: training-set composition, PWM parameters, real discriminative power.
usage: audit_signal.py <org> <variant>   -> writes <org>.<variant>.signal.json + text
"""
import sys, os, re, json, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pwmlib import *

org, variant = sys.argv[1], sys.argv[2]
OUT = f"{ROOT}/work/analysis/pwm_audit/out"
os.makedirs(OUT, exist_ok=True)
rund = f"{RUNS}/{org}/{variant}/pytd_dev_pwm/input.fa.transdecoder_dir"
seqs = read_fasta(f"{INPUTS}/{org}/{variant}.fa")
truth = read_truth(f"{INPUTS}/{org}/{variant}.truth.tsv")
res = {"org": org, "variant": variant}

# ---------------- Q1: training ORFs / positives -----------------
pat = re.compile(r"(\S+):(\d+)-(\d+)\(([+-])\)")
train = []
for l in open(f"{rund}/__checkpoints_TDpredict/top_training_orfs.cds"):
    if l[0] == ">":
        m = pat.search(l)
        train.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), "5prime_partial" in l or "internal" in l))
c = collections.Counter()
pos_feats = []
pos_class = []
for tid, a, b, strand, partial5 in train:
    s = seqs[tid]
    if strand == "-":
        s = rc(s)
        st = len(s) - max(a, b) + 1
    else:
        st = min(a, b)
    si = st - 1
    c["train_orfs"] += 1
    if s[si:si + 3] != "ATG":
        c["no_ATG_skipped"] += 1
        continue
    w = window(s, si)
    if w is None:
        c["too_little_upstream_skipped"] += 1
        continue
    c["positives"] += 1
    pos_feats.append(w)
    r = truth.get(tid)
    if r is None or strand != r["strand"]:
        cls = "no_truth_or_wrong_strand"
    else:
        ts = int(r["cds_start"]) - 1
        te = int(r["cds_end"])
        if si == ts:
            cls = "true_start"
        elif (max(a, b) == te or abs(max(a, b) - te) <= 2) and si < ts and (ts - si) % 3 == 0:
            cls = "upstream_inframe_ATG_same_stop"
        elif si > ts and (si - ts) % 3 == 0 and abs(max(a, b) - te) <= 2:
            cls = "downstream_inframe_ATG_same_stop"
        else:
            cls = "other_orf"
        if r["utr5_source"] != "annotated":
            cls += "|utr5_added"
    c["pos_" + cls] += 1
    pos_class.append(cls)
# compare with the features file written by the run
feat_file = [l.split()[0] for l in open(f"{rund}/start_refinement.+.features") if l.strip()]
c["features_file_n"] = len(feat_file)
c["features_match_reextraction"] = int(feat_file == pos_feats)
c["features_ATG_at_20"] = sum(1 for f in feat_file if f[20:23] == "ATG")
neg_file = [l.split()[0] for l in open(f"{rund}/start_refinement.-.features") if l.strip()]
c["neg_n"] = len(neg_file)
c["neg_ATG_at_20"] = sum(1 for f in neg_file if f[20:23] == "ATG")
enh = [l.split()[0] for l in open(f"{rund}/start_refinement.enhanced.+.features") if l.strip()]
c["enhanced_n"] = len(enh)
# which positives survived enhancement, by class
enh_set = set(enh)
cc = collections.Counter()
for f, cl in zip(pos_feats, pos_class):
    cc[(cl.split("|")[0], f in enh_set)] += 1
res["q1_counts"] = dict(c)
res["q1_enhanced_by_class"] = {f"{k[0]}|kept={k[1]}": v for k, v in cc.items()}

# negatives composition: classify negative windows from training ORFs by frame/region relative to the training ORF
nc = collections.Counter()
for tid, a, b, strand, partial5 in train:
    s = seqs[tid]
    if strand == "-":
        s = rc(s)
        st, en = len(s) - max(a, b) + 1, len(s) - min(a, b) + 1
    else:
        st, en = min(a, b), max(a, b)
    si = st - 1
    if s[si:si + 3] != "ATG":
        continue
    for i in find_all(s[st + 1:], "ATG"):
        p = st + 1 + i
        if window(s, p) is None:
            continue
        if p < en:
            nc["in_CDS_inframe" if (p - si) % 3 == 0 else "in_CDS_oof"] += 1
        else:
            nc["3UTR"] += 1
res["q1_negatives_composition"] = dict(nc)

# ---------------- Q2: PWM params ----------------
(up, down), thr, best_auc, best_f1 = select_range_threshold(
    f"{rund}/start_refinement.enhanced.feature.scores.roc.auc", f"{rund}/start_refinement.enhanced.feature.scores.roc")
(up0, down0), thr0, auc0, f10 = select_range_threshold(
    f"{rund}/start_refinement.feature.scores.roc.auc", f"{rund}/start_refinement.feature.scores.roc")
# threshold grid granularity
scores_all = [float(x.split("\t")[2]) for x in open(f"{rund}/start_refinement.enhanced.feature.scores") if not x.rstrip().endswith("NA")]
res["q2"] = dict(enh_range=[up, down], enh_threshold=thr, enh_auc_reported=best_auc, enh_bestF1=best_f1,
                 raw_range=[up0, down0], raw_threshold=thr0, raw_auc_reported=auc0, raw_bestF1=f10,
                 threshold_grid_step=(max(scores_all) - min(scores_all)) / 10, score_min=min(scores_all), score_max=max(scores_all))
# how many training positives (raw) pass the threshold with the enhanced pwm + chosen range
Pp = load_pwm_file(f"{rund}/start_refinement.enhanced.+.pwm")
Pm = load_pwm_file(f"{rund}/start_refinement.-.pwm")
Praw = load_pwm_file(f"{rund}/start_refinement.+.pwm")
Xpos = encode(pos_feats)
sp = score_range(Xpos, Pp, Pm, up, down)
res["q2"]["frac_raw_positives_pass_threshold"] = float((sp >= thr).mean())
Xneg = encode(neg_file)
sn = score_range(Xneg, Pp, Pm, up, down)
res["q2"]["frac_training_negatives_pass_threshold"] = float((sn >= thr).mean())
res["q2"]["auc_enh_on_raw_pos_vs_train_neg_chosen_range"] = auc(sp, sn)
res["q2"]["auc_enh_on_raw_pos_vs_train_neg_full"] = auc(score_range(Xpos, Pp, Pm, 20, 10), score_range(Xneg, Pp, Pm, 20, 10))

# ---------------- Q3: real signal using truth ----------------
sets, meta = truth_sets(seqs, truth, annotated_only=True)
setsA, _ = truth_sets(seqs, truth, annotated_only=False)
res["q3_set_sizes"] = {k: len(v) for k, v in sets.items()}
res["q3_set_sizes_incl_added_utr"] = {k: len(v) for k, v in setsA.items()}
X = {k: encode([w for _, w in v]) for k, v in sets.items()}
tids = {k: np.array([t for t, _ in v]) for k, v in sets.items()}

def auc_table(scoref):
    out = {}
    sp = scoref(X["pos"])
    for k in ("in_frame_down", "in_frame_up", "utr_oof", "down_any"):
        if len(X[k]):
            out[k] = round(auc(sp, scoref(X[k])), 4)
    return out

res["q3_auc_self_enhanced_chosen_range"] = auc_table(lambda A: score_range(A, Pp, Pm, up, down))
res["q3_auc_self_enhanced_full"] = auc_table(lambda A: score_range(A, Pp, Pm, 20, 10))
res["q3_auc_self_raw_full"] = auc_table(lambda A: score_range(A, Praw, Pm, 20, 10))
res["q3_auc_self_raw_kozak_6_3"] = auc_table(lambda A: score_range(A, Praw, Pm, 6, 3))
res["q3_auc_self_enh_kozak_6_3"] = auc_table(lambda A: score_range(A, Pp, Pm, 6, 3))

# truth-trained, 2-fold CV by transcript; negatives = several choices
rng = np.random.RandomState(0)
all_t = sorted(set(tids["pos"]))
fold = {t: rng.randint(2) for t in all_t}
fmask = {k: np.array([fold.get(t, 0) for t in tids[k]]) for k in X}

def cv_auc(neg_train_keys, ranges):
    """train P on pos[f], N on union(neg_train_keys)[f]; test on other fold; average AUC per contrast per range."""
    out = {r: collections.defaultdict(list) for r in ranges}
    kls = []
    for f in (0, 1):
        tr = 1 - f
        Ppos = pwm_from_feats([w for (t, w) in sets["pos"] if fold[t] == tr])
        negf = [w for k in neg_train_keys for (t, w) in sets[k] if fold.get(t, 0) == tr]
        Pneg = pwm_from_feats(negf)
        kls.append(kl_per_pos(Ppos, Pneg))
        for r in ranges:
            u, d = r
            sp = score_range(X["pos"][fmask["pos"] == f], Ppos, Pneg, u, d)
            for k in ("in_frame_down", "in_frame_up", "utr_oof", "down_any"):
                if (fmask[k] == f).sum():
                    out[r][k].append(auc(sp, score_range(X[k][fmask[k] == f], Ppos, Pneg, u, d)))
    return {f"{r[0]},{r[1]}": {k: round(float(np.mean(v)), 4) for k, v in d.items()} for r, d in out.items()}, np.mean(kls, 0)

ranges = [(20, 10), (6, 3), (3, 1), (10, 5), (6, 1), (10, 1), (3, 3)]
res["q3_auc_truth_TDneg_down_any"], kl_td = cv_auc(["down_any"], ranges)
res["q3_auc_truth_neg_inframe_down"], kl_ifd = cv_auc(["in_frame_down"], ranges)
res["q3_auc_truth_neg_inframe_up"], kl_ifu = cv_auc(["in_frame_up"], ranges) if len(sets["in_frame_up"]) > 30 else ({}, np.zeros(WLEN))
res["q3_auc_truth_neg_mixed_up_oof"], kl_mix = cv_auc(["in_frame_up", "utr_oof", "in_frame_down"], ranges)
res["q3_kl_by_pos"] = {
    "pos_vs_down_any": [round(float(v), 4) for v in kl_td],
    "pos_vs_inframe_down": [round(float(v), 4) for v in kl_ifd],
    "pos_vs_inframe_up": [round(float(v), 4) for v in kl_ifu],
}
# base composition at Kozak positions: -3 (idx 17), +4 (idx 23), -1 (19), -2 (18)
Ppos_all = pwm_from_feats([w for _, w in sets["pos"]])
res["q3_truth_pos_pwm_key_positions"] = {
    f"{lab}": {b: round(float(Ppos_all[i, j]), 3) for j, b in enumerate(BASES)}
    for lab, i in (("-3", 17), ("-2", 18), ("-1", 19), ("+4", 23), ("+5", 24))}
Pd = pwm_from_feats([w for _, w in sets["in_frame_down"]])
res["q3_inframe_down_pwm_key_positions"] = {
    f"{lab}": {b: round(float(Pd[i, j]), 3) for j, b in enumerate(BASES)}
    for lab, i in (("-3", 17), ("-2", 18), ("-1", 19), ("+4", 23), ("+5", 24))}
if len(sets["in_frame_up"]) > 30:
    Pu = pwm_from_feats([w for _, w in sets["in_frame_up"]])
    res["q3_inframe_up_pwm_key_positions"] = {
        f"{lab}": {b: round(float(Pu[i, j]), 3) for j, b in enumerate(BASES)}
        for lab, i in (("-3", 17), ("-2", 18), ("-1", 19), ("+4", 23), ("+5", 24))}
# self-trained PWMs at key positions
res["q2_self_enh_pwm_key_positions"] = {
    f"{lab}": {b: round(float(Pp[i, j]), 3) for j, b in enumerate(BASES)} for lab, i in (("-3", 17), ("-1", 19), ("+4", 23))}
res["q2_self_minus_pwm_key_positions"] = {
    f"{lab}": {b: round(float(Pm[i, j]), 3) for j, b in enumerate(BASES)} for lab, i in (("-3", 17), ("-1", 19), ("+4", 23))}
res["q2_self_raw_pwm_key_positions"] = {
    f"{lab}": {b: round(float(Praw[i, j]), 3) for j, b in enumerate(BASES)} for lab, i in (("-3", 17), ("-1", 19), ("+4", 23))}

# ---------------- Q4: candidate structure from truth ----------------
q4 = collections.Counter()
for tid, m in meta.items():
    q4["n"] += 1
    key = ("closed" if m["closed"] else "open") + ("_first" if m["n_up_inframe"] == 0 else "_later")
    q4[key] += 1
res["q4_truth_candidate_structure_annotated_utr"] = dict(q4)
q4a = collections.Counter()
_, metaA = truth_sets(seqs, truth, annotated_only=False)
for tid, m in metaA.items():
    q4a["n"] += 1
    key = ("closed" if m["closed"] else "open") + ("_first" if m["n_up_inframe"] == 0 else "_later")
    q4a[key] += 1
res["q4_truth_candidate_structure_all"] = dict(q4a)

json.dump(res, open(f"{OUT}/{org}.{variant}.signal.json", "w"), indent=1)
print(json.dumps(res, indent=1))
