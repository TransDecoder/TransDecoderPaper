#!/usr/bin/env python3
"""ORFs retained with one homology database but not another, with their best Swiss-Prot hits.

usage: homology_series_retention.py <pred_status A> <pred_status B> <blastp outfmt6 used for A>
Prints the class counts of distinct CDSs present in A but not in B, and for the novel ones (no GENCODE
CDS overlap) the species and protein-name mnemonics of their best hit, and the hit regions on the
most frequent subject.
"""
import collections
import sys


def status(p):
    d = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if f[0] == "model_id":
            continue
        c = "exact" if f[3] == "1" else "compatible" if f[4] == "1" else "novel" if f[5] == "1" else "other"
        d[f[0]] = (c, int(f[2]))
    return d


a, b = status(sys.argv[1]), status(sys.argv[2])
hits = {}
for line in open(sys.argv[3]):
    f = line.split("\t")
    hits.setdefault(f[0], f)
only = [m for m in a if m not in b]
print("in A only:", len(only), dict(collections.Counter(a[m][0] for m in only)),
      "short(<100aa):", sum(a[m][1] < 100 for m in only))
sp, prot, region = collections.Counter(), collections.Counter(), collections.defaultdict(collections.Counter)
for m in only:
    if a[m][0] != "novel":
        continue
    h = hits.get(m)
    if not h:
        sp["nohit"] += 1
        continue
    mn = h[1].split("|")[-1]
    sp[mn.split("_")[-1]] += 1
    prot[mn] += 1
    region[mn][f"{h[8]}-{h[9]}"] += 1
print("novel, best-hit species:", sp.most_common(10))
print("novel, best-hit proteins:", prot.most_common(15))
top = prot.most_common(1)[0][0]
print(f"hit regions on {top}:", region[top].most_common(5))
