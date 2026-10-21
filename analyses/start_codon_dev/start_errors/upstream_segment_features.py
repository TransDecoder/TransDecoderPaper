#!/usr/bin/env python3
"""Features of the segment between an open 5' ORF end and its in-frame ATGs.

For each 5'-partial call (no ATG at the predicted start) whose stop matches the truth, in a
stranded benchmark variant, write: the truth label (keep = truth start within 2 nt of the
open end; atg<k> = truth start is the k-th in-frame ATG; other), and for the first few
in-frame ATGs the length and Markov log-likelihood ratio (the run's own hexamer.scores, scored
as TransDecoder does) of the segment from the open end up to that ATG, plus the LLR of the
remainder of the ORF.

usage: upstream_segment_features.py <inputs_dir> <runs_dir> <variant> <tool> > out.tsv
"""
import sys, os
inputs, runs, variant, tool = sys.argv[1:5]
STOPS = {"TAA", "TAG", "TGA"}
MAXK = 5

def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l[0] == ">":
            if n: d[n] = "".join(b).upper()
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = "".join(b).upper()
    return d

def llr_profile(seq, scores, order=5, shift=0):
    out = []
    for i in range(len(seq)):
        mu = min(i, order)
        out.append(scores.get(f"{seq[i-mu:i+1]}-{(i+shift)%3}", 0.0))
    return out

hdr = ["org", "variant", "transcript_id", "label", "orf_len", "truth_type"]
for k in range(1, MAXK + 1):
    hdr += [f"atg{k}_dist", f"atg{k}_seg_llr", f"atg{k}_rest_llr", f"atg{k}_seg_s1", f"atg{k}_seg_s2"]
print("\t".join(hdr))
for org in sorted(os.listdir(inputs)):
    tf = f"{inputs}/{org}/{variant}.truth.tsv"
    rd = f"{runs}/{org}/{variant}/{tool}"
    pf = f"{rd}/preds.tsv"
    hx = f"{runs}/{org}/{variant}/pytd/input.fa.transdecoder_dir/hexamer.scores"
    if not (os.path.exists(tf) and os.path.exists(pf) and os.path.exists(hx)):
        continue
    stops = {"TGA"} if org == "tthermophila" else STOPS
    scores = {}
    for l in open(hx):
        if l[0] == "#": continue
        x = l.split("\t"); scores[x[0]] = float(x[3])
    seqs = fasta(f"{inputs}/{org}/{variant}.fa")
    truth = {}
    for l in open(tf):
        x = l.rstrip("\n").split("\t")
        if x[0] == "transcript_id": continue
        truth[x[0]] = (int(x[1]), int(x[2]), x[3], x[7])
    seen = set()
    for l in open(pf):
        x = l.rstrip("\n").split("\t")
        tid = x[0]
        if tid not in truth or tid in seen or x[1] == "start": continue
        ts, te, tst, ttype = truth[tid]; ps, pe = int(x[1]), int(x[2])
        if x[3] != "+" or tst != "+" or abs(pe - te) > 2: continue
        s = seqs[tid]
        if s[ps-1:ps+2] == "ATG": continue
        seen.add(tid)
        orf = s[ps-1:pe]
        prof = llr_profile(orf, scores)
        prof1 = llr_profile(orf, scores, shift=1)
        prof2 = llr_profile(orf, scores, shift=2)
        total = sum(prof)
        atgs = [i for i in range(3, len(orf) - 3, 3) if orf[i:i+3] == "ATG"][:MAXK]
        if abs(ts - ps) <= 2: label = "keep"
        else:
            label = "other"
            for k, i in enumerate(atgs, 1):
                if abs(ps + i - ts) <= 2: label = f"atg{k}"; break
        row = [org, variant, tid, label, str(len(orf)), ttype]
        for k in range(MAXK):
            if k < len(atgs):
                i = atgs[k]; seg = sum(prof[:i])
                row += [str(i), f"{seg:.2f}", f"{total-seg:.2f}", f"{sum(prof1[:i]):.2f}", f"{sum(prof2[:i]):.2f}"]
            else:
                row += ["NA"] * 5
        print("\t".join(row))
