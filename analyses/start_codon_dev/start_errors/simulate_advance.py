#!/usr/bin/env python3
"""Move ORF starts downstream across in-frame ATGs while the skipped segment lacks coding
frame contrast. Applies to every prediction (closed ORFs too).

usage: simulate_advance.py <input.fa> <transdecoder_dir> <preds.tsv> <alpha> <mode> > new.tsv
mode: once (at most one move) | repeat (keep advancing) ; alpha scales the training ORFs'
per-base contrast (segment kept as coding if contrast >= alpha*mu*len).
Only ATG-start ORFs are advanced (open 5' ends were already handled by start refinement).
"""
import sys, statistics
fa, tdd, pf, alpha, mode = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
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
scores = {}
for l in open(f"{tdd}/hexamer.scores"):
    if l[0] == "#": continue
    x = l.split("\t"); scores[x[0]] = float(x[3])
def S(seq):
    return sum(scores.get(f"{seq[max(0,i-5):i+1]}-{i%3}", 0.0) for i in range(len(seq)))
def contrast(seq):
    return S(seq) - max(S(seq[1:]), S(seq[2:]))
seqs = fasta(fa)
rows = [l.rstrip("\n").split("\t") for l in open(pf)]
hdr, rows = rows[0], rows[1:]
def oriented(r):
    s = seqs[r[0]]; a, b = int(r[1]), int(r[2])
    if r[3] == "-":
        s = s.translate(COMP)[::-1]; a, b = len(s) - b + 1, len(s) - a + 1
    return s, a, b
cand = []
for r in rows:
    s, a, b = oriented(r)
    if s[a-1:a+2] == "ATG": cand.append((b - a, s[a-1:b]))
cand.sort(reverse=True)
mu = statistics.median(contrast(c) / len(c) for _, c in cand[:500])
print("\t".join(hdr)); nmove = 0
for r in rows:
    s, a, b = oriented(r)
    if s[a-1:a+2] == "ATG":
        orf = s[a-1:b]; off = 0
        while True:
            nxt = next((i for i in range(off + 3, len(orf) - 3, 3) if orf[i:i+3] == "ATG"), None)
            if nxt is None or len(orf) - nxt < 300: break
            seg = orf[off:nxt]
            if contrast(seg) < alpha * mu * len(seg):
                off = nxt; nmove += 1
                if mode == "once": break
            else:
                break
        if off:
            a += off
            if r[3] == "-":
                L = len(s); r[1], r[2] = str(L - b + 1), str(L - a + 1)
            else:
                r[1] = str(a)
    print("\t".join(r))
print(f"mu={mu:.4f} moved={nmove}", file=sys.stderr)
