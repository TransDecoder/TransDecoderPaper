#!/usr/bin/env python3
"""Offline simulation: Kozak-window PWM to move a predicted ATG start to a later in-frame ATG.

usage: sim_kozak.py <inputs_dir> <runs_dir> <org> <variant> <tool>
Prints TSV rows: org variant window margin gate n_stop_tp exact_base exact_new moved
Training (self, no truth): the 500 longest ATG-start predictions; positives = their start
windows; negatives = in-frame ATGs inside their CDS. Window = U nt upstream + D nt after the ATG
(ATG itself excluded). For each prediction starting with ATG, candidates are the next in-frame
ATGs (up to 5, leaving >=300 nt); a candidate wins if its log-odds exceeds the current start's
by > margin and (gate) the skipped segment's frame contrast < gate * mu * len.
"""
import sys, math, statistics, collections
inputs, runs, org, variant, tool = sys.argv[1:6]
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l[0] == ">":
            if n: d[n] = "".join(b).upper()
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = "".join(b).upper()
    return d
seqs = fasta(f"{inputs}/{org}/{variant}.fa")
tdd = f"{runs}/{org}/{variant}/{tool}/input.fa.transdecoder_dir"
hx = {}
for l in open(f"{tdd}/hexamer.scores"):
    if l[0] == "#": continue
    x = l.split("\t"); hx[x[0]] = float(x[3])
def S(seq): return sum(hx.get(f"{seq[max(0,i-5):i+1]}-{i%3}", 0.0) for i in range(len(seq)))
def contrast(seq): return S(seq) - max(S(seq[1:]), S(seq[2:]))
truth = {}
for l in open(f"{inputs}/{org}/{variant}.truth.tsv"):
    x = l.rstrip("\n").split("\t")
    if x[0] == "transcript_id": continue
    truth[x[0]] = (int(x[1]), int(x[2]), x[3])
preds = []
for l in open(f"{runs}/{org}/{variant}/{tool}/preds.tsv"):
    x = l.rstrip("\n").split("\t")
    if x[0] == "transcript_id": continue
    tid, a, b, st = x[0], int(x[1]), int(x[2]), x[3]
    s = seqs[tid]
    if st == "-":
        s = rc(s); a, b = len(s) - b + 1, len(s) - a + 1
    preds.append((tid, st, s, a - 1, b))          # 0-based start, end exclusive (oriented)
# truth in oriented coords of the predicted strand
def truth_oriented(tid, st, L):
    if tid not in truth: return None
    a, b, tst = truth[tid]
    if tst != st: return None
    return (a - 1, b) if st == "+" else (L - b, L - a + 1)
train = sorted([p for p in preds if p[2][p[3]:p[3]+3] == "ATG"], key=lambda p: p[3] - p[4])[:500]
mu = statistics.median(contrast(p[2][p[3]:p[4]]) / (p[4] - p[3]) for p in train)
WINDOWS = [(6, 2), (9, 3), (4, 1), (12, 4)]
def window(s, i, U, D):
    # positions i-U..i-1 and i+3..i+3+D-1; None where outside the sequence
    return [s[j] if 0 <= j < len(s) else None for j in list(range(i - U, i)) + list(range(i + 3, i + 3 + D))]
def build(feats, n):
    c = [collections.Counter() for _ in range(n)]
    for f in feats:
        for k, b in enumerate(f):
            if b and b in "ACGT": c[k][b] += 1
    return [{b: (cc[b] + 1) / (sum(cc.values()) + 4) for b in "ACGT"} for cc in c]
def inframe_atgs(s, a, b):
    return [i for i in range(a + 3, b - 3, 3) if s[i:i+3] == "ATG"]
results = []
for U, D in WINDOWS:
    pos, neg = [], []
    for tid, st, s, a, b in train:
        pos.append(window(s, a, U, D))
        neg.extend(window(s, i, U, D) for i in inframe_atgs(s, a, b))
    P, N = build(pos, U + D), build(neg, U + D)
    def llr(s, i):
        return sum(math.log(P[k][c] / N[k][c]) for k, c in enumerate(window(s, i, U, D)) if c and c in "ACGT")
    cache = []
    for tid, st, s, a, b in preds:
        t = truth_oriented(tid, st, len(s))
        stop_ok = t is not None and abs(t[1] - b) <= 2
        if not stop_ok: continue
        base_exact = abs(t[0] - a) <= 2
        cands = []
        if s[a:a+3] == "ATG":
            l0 = llr(s, a)
            for i in inframe_atgs(s, a, b)[:5]:
                if b - i < 300: break
                cands.append((i, llr(s, i) - l0, contrast(s[a:i]) / (mu * (i - a))))
        cache.append((t, a, base_exact, cands))
    for margin in [0, 1, 2, 3, 4, 6]:
        for gate in [None, 0.4, 1.0]:
            ex = mv = 0
            for t, a, base_exact, cands in cache:
                new = a
                ok = [c for c in cands if c[1] > margin and (gate is None or c[2] < gate)]
                if ok:
                    new = max(ok, key=lambda c: c[1])[0]; mv += 1
                ex += abs(t[0] - new) <= 2
            base = sum(c[2] for c in cache)
            results.append((f"{U},{D}", margin, gate, len(cache), base, ex, mv))
for w, m, g, n, base, ex, mv in results:
    print(f"{org}\t{variant}\t{w}\t{m}\t{g}\t{n}\t{base}\t{ex}\t{mv}")
