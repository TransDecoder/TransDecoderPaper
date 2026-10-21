#!/usr/bin/env python3
"""Keep-vs-trim for open 5' ends using frame contrast plus Kozak log-odds of the first ATG.

usage: trim_kozak.py <org>   (uses useg.{complete,fragments}.tsv rows, labels keep / atg1)
Features per row: r = contrast / (mu * dist)  (as in the shipped rule: trim if r < alpha)
                  K = Kozak (6,2) log-odds of the first in-frame ATG, self-trained per run.
Rule: trim if r < alpha + beta * K. Prints per (alpha, beta) correct counts.
"""
import sys, math, collections, statistics
import pandas as pd
org = sys.argv[1]
I, R = "data/bench_v3/inputs", "work/bench_v3/runs"
COMP = str.maketrans("ACGT", "TGCA")
def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l[0] == ">":
            if n: d[n] = "".join(b).upper()
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = "".join(b).upper()
    return d
U, D = 6, 2
def window(s, i): return [s[j] if 0 <= j < len(s) else None for j in list(range(i - U, i)) + list(range(i + 3, i + 3 + D))]
def build(feats):
    c = [collections.Counter() for _ in range(U + D)]
    for f in feats:
        for k, b in enumerate(f):
            if b: c[k][b] += 1
    return [{b: (cc[b] + 1) / (sum(cc.values()) + 4) for b in "ACGT"} for cc in c]
rows = []
for variant in ["complete", "fragments"]:
    d = pd.read_csv(f"work/analysis/start_errors/useg.{variant}.tsv", sep="\t")
    d = d[(d.org == org) & d.label.isin(["keep", "atg1"]) & d.atg1_dist.notna()]
    seqs = fasta(f"{I}/{org}/{variant}.fa")
    tdd = f"{R}/{org}/{variant}/pytd/input.fa.transdecoder_dir"
    hx = {}
    for l in open(f"{tdd}/hexamer.scores"):
        if l[0] != "#": x = l.split("\t"); hx[x[0]] = float(x[3])
    def S(q): return sum(hx.get(f"{q[max(0,i-5):i+1]}-{i%3}", 0.0) for i in range(len(q)))
    def contrast(q): return S(q) - max(S(q[1:]), S(q[2:]))
    preds = pd.read_csv(f"{R}/{org}/{variant}/pytd/preds.tsv", sep="\t")
    preds = preds[preds.strand == "+"].copy(); preds["L"] = preds.end - preds.start
    train = []
    for r in preds.sort_values("L", ascending=False).itertuples():
        s = seqs[r.transcript_id]
        if s[r.start-1:r.start+2] == "ATG": train.append((s, r.start - 1, r.end))
        if len(train) >= 500: break
    mu = statistics.median(contrast(s[a:b]) / (b - a) for s, a, b in train)
    pos, neg = [], []
    for s, a, b in train:
        pos.append(window(s, a)); neg.extend(window(s, i) for i in range(a + 3, b - 3, 3) if s[i:i+3] == "ATG")
    P, N = build(pos), build(neg)
    pn = pd.read_csv(f"{R}/{org}/{variant}/pytd_norefine/preds.tsv", sep="\t").drop_duplicates("transcript_id").set_index("transcript_id")
    for r in d.itertuples():
        s = seqs[r.transcript_id]; st = int(pn.loc[r.transcript_id, "start"]) - 1; i = st + int(r.atg1_dist)
        seg = s[st:i]
        K = sum(math.log(P[k][c] / N[k][c]) for k, c in enumerate(window(s, i)) if c)
        rows.append((variant, r.label, contrast(seg) / (mu * len(seg)), K))
df = pd.DataFrame(rows, columns=["variant", "label", "r", "K"])
base = ((df.r < 0.4) == (df.label == "atg1")).sum()
print(f"{org}\tn={len(df)}\ttrim_frac={(df.label=='atg1').mean():.2f}\tcurrent(alpha=0.4)={base}")
for a in [0.2, 0.3, 0.4, 0.5]:
    for b in [0, 0.05, 0.1, 0.15, 0.2, 0.3]:
        c = ((df.r < a + b * df.K) == (df.label == "atg1")).sum()
        print(f"{org}\t{a}\t{b}\t{c}\t{c-base:+d}")
