#!/usr/bin/env python3
"""For truth CDSs whose stop a tool got right, classify the predicted start.

usage: classify_starts.py <inputs_dir> <runs_dir> <variant> <tool> [<tool> ...]
Prints, per organism and tool, counts of: exact; upstream to a non-ATG at the transcript
5' end (5'-partial call); upstream ATG; downstream ATG; other.
Stranded variants only (all truth on +).
"""
import sys, os, collections
inputs, runs, variant, tools = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]

def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l[0] == ">":
            if n: d[n] = "".join(b)
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = "".join(b)
    return d

print("org\ttool\tstop_tp\texact\tup_partial5\tup_ATG\tdown_ATG\tother\tfrac_exact\tfrac_up_partial5")
for org in sorted(os.listdir(inputs)):
    tf = f"{inputs}/{org}/{variant}.truth.tsv"
    if not os.path.exists(tf): continue
    seqs = fasta(f"{inputs}/{org}/{variant}.fa")
    truth = {}
    for l in open(tf):
        x = l.rstrip("\n").split("\t")
        if x[0] == "transcript_id": continue
        truth[x[0]] = (int(x[1]), int(x[2]), x[3])
    for t in tools:
        pf = f"{runs}/{org}/{variant}/{t}/preds.tsv"
        if not os.path.exists(pf): continue
        c = collections.Counter()
        for l in open(pf):
            x = l.split("\t")
            if x[0] not in truth or x[1] == "start": continue
            ts, te, st = truth[x[0]]; ps, pe, pst = int(x[1]), int(x[2]), x[3]
            if pst != st or abs(pe - te) > 2: continue
            if c.get(("seen", x[0])): continue
            c[("seen", x[0])] = 1
            c["stop_tp"] += 1
            s = seqs[x[0]]
            if abs(ps - ts) <= 2: c["exact"] += 1
            elif ps < ts and s[ps-1:ps+2] != "ATG": c["up_partial5"] += 1
            elif ps < ts: c["up_ATG"] += 1
            elif s[ps-1:ps+2] == "ATG": c["down_ATG"] += 1
            else: c["other"] += 1
        n = c["stop_tp"] or 1
        print(f"{org}\t{t}\t{c['stop_tp']}\t{c['exact']}\t{c['up_partial5']}\t{c['up_ATG']}\t{c['down_ATG']}\t{c['other']}\t{c['exact']/n:.3f}\t{c['up_partial5']/n:.3f}")
