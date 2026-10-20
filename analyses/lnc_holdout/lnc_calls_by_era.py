#!/usr/bin/env python3
"""lncRNA call rate per tool, split by whether the lncRNA's gene existed in an older GENCODE release.

CPAT, CPC2 and RNAsamba were trained on GENCODE lncRNAs from their time; lncRNA genes added
later are a held-out set for them.

usage: lnc_calls_by_era.py <lnc.gtf> <old_gene_ids.txt> <runs_dir/org> <variant> <tool> [<tool> ...]
"""
import re, sys
from pathlib import Path

gtf, old_ids, runs, variant, tools = sys.argv[1], sys.argv[2], Path(sys.argv[3]), sys.argv[4], sys.argv[5:]
old = set(open(old_ids).read().split())
era = {}
for line in open(gtf):
    m = re.search(r'gene_id "([^".]+)[^"]*"; transcript_id "([^"]+)"', line)
    if m:
        era[m.group(2)] = "pre_v28" if m.group(1) in old else "post_v28"
n = {e: sum(1 for v in era.values() if v == e) for e in ("pre_v28", "post_v28")}
print("tool\tera\tn_lnc\tn_called\tfrac_called")
for t in tools:
    called = set()
    with open(runs / variant / t / "preds.tsv") as f:
        next(f)
        for line in f:
            tx = line.split("\t", 1)[0]
            if tx.startswith("lnc_"):
                called.add(tx)
    for e in ("pre_v28", "post_v28"):
        c = sum(1 for tx in called if era.get(tx) == e)
        print(f"{t}\t{e}\t{n[e]}\t{c}\t{c / n[e]:.4f}")
