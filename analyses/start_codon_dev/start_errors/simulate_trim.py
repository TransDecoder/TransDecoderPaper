#!/usr/bin/env python3
"""Post-process a run's predictions with the frame-contrast trim rule and write new preds.

usage: simulate_trim.py <input.fa> <transdecoder_dir> <preds.tsv> <alpha> [tetrahymena] > new_preds.tsv
For each + or - strand prediction whose start is not ATG (open 5' end), compute the frame
contrast C = S0 - max(S1, S2) of the segment from the open end to the first in-frame ATG (S_k:
TransDecoder Markov score with frame labels shifted by k). mu = median per-nt contrast over the
longest 500 ATG-start predictions (training-like). Trim to the ATG if C < alpha * mu * L.
"""
import sys, statistics
fa, tdd, pf, alpha = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
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
def S(seq, shift):
    return sum(scores.get(f"{seq[max(0,i-5):i+1]}-{(i+shift)%3}", 0.0) for i in range(len(seq)))
def contrast(seq):
    return S(seq, 0) - max(S(seq, 1), S(seq, 2))
seqs = fasta(fa)
rows = [l.rstrip("\n").split("\t") for l in open(pf)]
hdr, rows = rows[0], rows[1:]
def oriented(r):
    s = seqs[r[0]]; a, b = int(r[1]), int(r[2])
    if r[3] == "-":
        s = s.translate(COMP)[::-1]; a, b = len(s) - b + 1, len(s) - a + 1
    return s, a, b
# mu from training-like ORFs
cand = []
for r in rows:
    s, a, b = oriented(r)
    if s[a-1:a+2] == "ATG": cand.append((b - a, s[a-1:b]))
cand.sort(reverse=True)
mu = statistics.median(contrast(c) / len(c) for _, c in cand[:500])
print(f"mu={mu:.4f}", file=sys.stderr)
print("\t".join(hdr))
ntrim = 0
for r in rows:
    s, a, b = oriented(r)
    if s[a-1:a+2] != "ATG":
        orf = s[a-1:b]
        i = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i+3] == "ATG"), None)
        if i is not None and contrast(orf[:i]) < alpha * mu * i:
            a += i; ntrim += 1
            if r[3] == "-":
                L = len(s); r[1], r[2] = str(L - b + 1), str(L - a + 1)
            else:
                r[1] = str(a)
    print("\t".join(r))
print(f"trimmed={ntrim}", file=sys.stderr)
