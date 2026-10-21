"""Shared helpers for the start-codon PWM audit (read-only w.r.t. the code trees)."""
import os, re, math, collections
import numpy as np
from scipy.stats import rankdata

ROOT = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
INPUTS = f"{ROOT}/data/bench_v3/inputs"
RUNS = f"{ROOT}/work/bench_v3/runs"
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii"]
STOPS = {1: {"TAA", "TAG", "TGA"}, 6: {"TGA"}}

BASES = "GATC"
BIDX = {b: i for i, b in enumerate(BASES)}
PSEUDO = 0.1
LEFT, RIGHT = 20, 10
WLEN = LEFT + 3 + RIGHT
ATG = LEFT

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l[0] == ">":
            if n:
                d[n] = "".join(b).upper()
            n, b = l[1:].split()[0], []
        else:
            b.append(l.strip())
    if n:
        d[n] = "".join(b).upper()
    return d


def read_truth(p):
    t = {}
    hdr = None
    for l in open(p):
        x = l.rstrip("\n").split("\t")
        if hdr is None:
            hdr = x
            continue
        r = dict(zip(hdr, x))
        t[r["transcript_id"]] = r
    return t


def window(seq, atg_idx):
    """31-nt window with ATG at index 20, or None."""
    b = atg_idx - LEFT
    e = b + WLEN
    if b < 0 or e > len(seq):
        return None
    w = seq[b:e]
    if any(c not in BASES for c in w):
        return None
    return w


def encode(feats):
    """list of equal-length strings -> int array (n, L)."""
    arr = np.frombuffer("".join(feats).encode(), dtype=np.uint8).reshape(len(feats), -1)
    lut = np.full(256, -1, dtype=np.int8)
    for b, i in BIDX.items():
        lut[ord(b)] = i
    return lut[arr]


def pwm_from_feats(feats, pseudo=PSEUDO):
    X = encode(feats)
    L = X.shape[1]
    counts = np.zeros((L, 4))
    for j in range(4):
        counts[:, j] = (X == j).sum(0)
    return (counts + pseudo) / (counts.sum(1, keepdims=True) + 4 * pseudo)


def load_pwm_file(p):
    rows = []
    for l in open(p):
        x = l.rstrip().split("\t")
        if x[0] == "pos":
            continue
        rows.append([float(v) for v in x[1:5]])
    return np.array(rows)  # columns G A T C


def score_matrix(X, pwm_plus, pwm_minus, lo=0, hi=None, norm=True):
    """sum_{i in [lo,hi]} log(p+/p-) ; per-feature. X int-encoded (n,L). hi inclusive."""
    L = X.shape[1]
    if hi is None:
        hi = L - 1
    lo_mat = np.log(pwm_plus / pwm_minus)  # (L,4)
    idx = np.arange(lo, hi + 1)
    s = lo_mat[idx[:, None], X[:, idx].T].sum(0)
    if norm:
        s = s / (hi - lo + 1 - 3)
    return s


def score_range(X, pwm_plus, pwm_minus, up, down, norm=True):
    return score_matrix(X, pwm_plus, pwm_minus, ATG - up, ATG + 2 + down, norm=norm)


def auc(pos, neg):
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    r = rankdata(np.concatenate([pos, neg]))
    return (r[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def select_range_threshold(auc_file, roc_file):
    """Replicates predict._select_pwm_range_and_threshold."""
    best, best_auc = None, -1e9
    for l in open(auc_file):
        if not l.strip():
            continue
        r, a = l.split("\t")
        if float(a) > best_auc:
            best_auc, best = float(a), r
    thr, bf1 = None, -1e9
    with open(roc_file) as fh:
        next(fh)
        for l in fh:
            x = l.rstrip().split("\t")
            if x[0] != best:
                continue
            if float(x[8]) > bf1:
                bf1, thr = float(x[8]), float(x[1])
    up, down = [int(v) for v in best.split(",")]
    return (up, down), thr, best_auc, bf1


def kl_per_pos(pwm_a, pwm_b):
    return (pwm_a * np.log(pwm_a / pwm_b)).sum(1)


def find_all(seq, motif):
    i = 0
    out = []
    while True:
        i = seq.find(motif, i)
        if i < 0:
            return out
        out.append(i)
        i += 1


def truth_sets(seqs, truth, annotated_only=True, max_down_inframe=None):
    """For stranded truth (+): returns dict of feature lists keyed by class, each list of (tid, window).
    pos: true start.  in_frame_down: in-frame ATGs inside CDS after start.
    in_frame_up: in-frame ATGs upstream of the true start with no in-frame stop in between (N-ext decoys).
    utr_oof: out-of-frame ATGs in the 5' UTR.  down_any: all ATGs (any frame) downstream of start+2 (TD-style).
    """
    sets = collections.defaultdict(list)
    meta = {}
    for tid, r in truth.items():
        if annotated_only and r["utr5_source"] != "annotated":
            continue
        if r["strand"] != "+":
            continue
        s = seqs[tid]
        cs, ce = int(r["cds_start"]) - 1, int(r["cds_end"])  # 0-based start idx, end exclusive
        stops = STOPS[int(r["transl_table"])]
        if s[cs:cs + 3] != "ATG":
            continue
        w = window(s, cs)
        if w is None:
            continue
        sets["pos"].append((tid, w))
        # in-frame downstream within CDS
        nd = 0
        for i in range(cs + 3, ce - 3, 3):
            if s[i:i + 3] == "ATG":
                w2 = window(s, i)
                if w2:
                    sets["in_frame_down"].append((tid, w2))
                    nd += 1
                    if max_down_inframe and nd >= max_down_inframe:
                        break
        # in-frame upstream ATGs with no stop between
        i = cs - 3
        ups = []
        while i >= 0:
            c = s[i:i + 3]
            if c in stops:
                break
            if c == "ATG":
                ups.append(i)
            i -= 3
        closed = i >= 0
        for i in ups:
            w2 = window(s, i)
            if w2:
                sets["in_frame_up"].append((tid, w2))
        # out-of-frame upstream ATGs in 5' UTR
        for i in find_all(s[:cs], "ATG"):
            if (cs - i) % 3 == 0:
                continue
            w2 = window(s, i)
            if w2:
                sets["utr_oof"].append((tid, w2))
        # TD-style: all ATGs any frame downstream (from cs+2 on, as in the code)
        for i in find_all(s[cs + 2:], "ATG"):
            w2 = window(s, cs + 2 + i)
            if w2:
                sets["down_any"].append((tid, w2))
        meta[tid] = dict(closed=closed, n_up_inframe=len(ups), n_down_inframe=nd)
    return sets, meta
