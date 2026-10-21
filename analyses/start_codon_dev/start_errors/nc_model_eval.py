#!/usr/bin/env python3
"""Coding-vs-UTR log-likelihood of the open-5'-end segment, using a 5th-order non-coding Markov
model trained on 3'-UTR sequence (downstream of the stop) of the run's training ORFs.

usage: nc_model_eval.py <inputs_dir> <runs_dir> <useg.tsv ...>  (adds column atg1_nc_llr; prints TSV)
"""
import sys, os, math, collections
import pandas as pd
inputs, runs, files = sys.argv[1], sys.argv[2], sys.argv[3:]
ORDER = 5

def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l[0] == ">":
            if n: d[n] = "".join(b).upper()
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = "".join(b).upper()
    return d

def train_nc(seqs):
    k, km1 = collections.Counter(), collections.Counter()
    for s in seqs:
        for i in range(len(s)):
            for m in range(0, min(i, ORDER) + 1):
                k[s[i-m:i+1]] += 1; km1[s[i-m:i]] += 1
    return k, km1

def nc_logp(seq, model):
    k, km1 = model; tot = 0.0
    for i in range(len(seq)):
        m = min(i, ORDER)
        ctx = seq[i-m:i]; kmer = seq[i-m:i+1]
        tot += math.log((k[kmer] + 1) / (km1[ctx] + 4))
    return tot

def coding_logp(seq, scores, bg):
    tot = 0.0
    for i in range(len(seq)):
        m = min(i, ORDER)
        tot += scores.get(f"{seq[i-m:i+1]}-{i%3}", 0.0) + math.log(bg.get(seq[i], 0.25))
    return tot

d = pd.concat([pd.read_csv(f, sep="\t") for f in files])
out = []
for (org, variant), g in d.groupby(["org", "variant"]):
    rd = f"{runs}/{org}/{variant}/pytd/input.fa.transdecoder_dir"
    scores = {}
    for l in open(f"{rd}/hexamer.scores"):
        if l[0] == "#": continue
        x = l.split("\t"); scores[x[0]] = float(x[3])
    bg = {}
    for l in open(f"{rd}/base_freqs.dat"):
        x = l.split(); bg[x[0]] = float(x[2])
    seqs = fasta(f"{inputs}/{org}/{variant}.fa")
    # training-like ORFs: the 500 longest complete predictions of the pytd run; UTR = after the stop
    preds = pd.read_csv(f"{runs}/{org}/{variant}/pytd/preds.tsv", sep="\t")
    preds = preds[preds.strand == "+"].copy(); preds["len"] = preds.end - preds.start
    utrs = []
    for r in preds.sort_values("len", ascending=False).itertuples():
        s = seqs[r.transcript_id]
        if s[r.start-1:r.start+2] != "ATG": continue
        u = s[r.end:]
        if len(u) >= 30: utrs.append(u)
        if len(utrs) >= 500: break
    model = train_nc(utrs)
    vals = []
    for r in g.itertuples():
        s = seqs[r.transcript_id]
        # the ORF begins at the predicted (open) start; reconstruct: segment = first atg1_dist nt of the ORF
        pr = preds[preds.transcript_id == r.transcript_id] if False else None
        vals.append(None)
    g = g.copy()
    # segment sequences need the ORF start; recompute from truth-independent preds of pytd_norefine
    pn = pd.read_csv(f"{runs}/{org}/{variant}/pytd_norefine/preds.tsv", sep="\t").drop_duplicates("transcript_id").set_index("transcript_id")
    res = []
    for r in g.itertuples():
        st = int(pn.loc[r.transcript_id, "start"]) if r.transcript_id in pn.index else None
        if st is None or pd.isna(r.atg1_dist):
            res.append(float("nan")); continue
        seg = seqs[r.transcript_id][st-1:st-1+int(r.atg1_dist)]
        res.append(coding_logp(seg, scores, bg) - nc_logp(seg, model))
    g["atg1_nc_llr"] = res
    g["n_utr_train"] = len(utrs)
    out.append(g)
pd.concat(out).to_csv(sys.stdout, sep="\t", index=False)
