#!/usr/bin/env python3
"""Compare Perl and PyTransDecoder predictions on the 84 benchmark inputs at the release commits
(tools perl_td_v7 = Perl 41fc5ee, pytd_v7 = PyTD 3a793ff; evidence E25).

usage: 01_perl_pytd_identity.py [P]      (P = /broad/hptmp/bhaas/TRANSDECODER_PAPER)
For each input in work/bench_v3/runs/<org>/<variant>/, reads preds.tsv of both tools, keys each
prediction by (transcript, start, end, strand) and writes
data/benchmark/final_v7/perl_pytd_identity.tsv: counts per tool, number shared, number only in
one tool, whether the sorted key sets are identical, and the wall-clock run times (runtime_sec.txt).
Run with $P/env/td/bin/python.
"""
import sys
from pathlib import Path

P = Path(sys.argv[1] if len(sys.argv) > 1 else "/broad/hptmp/bhaas/TRANSDECODER_PAPER")
RUNS = P / "work/bench_v3/runs"
OUT = P / "TransDecoderPaper/data/benchmark/final_v7/perl_pytd_identity.tsv"


def keys(f):
    with open(f) as fh:
        next(fh)
        return {tuple(l.rstrip("\n").split("\t")[:4]) for l in fh if l.strip()}


def rt(d):
    f = d / "runtime_sec.txt"
    return f.read_text().strip() if f.exists() else "NA"


rows = []
for py in sorted(RUNS.glob("*/*/pytd_v7/preds.tsv")):
    d = py.parent.parent
    pl = d / "perl_td_v7/preds.tsv"
    a, b = keys(pl), keys(py)
    rows.append([d.parent.name, d.name, len(a), len(b), len(a & b), len(a - b), len(b - a), a == b,
                 rt(d / "perl_td_v7"), rt(d / "pytd_v7")])
with open(OUT, "w") as fh:
    fh.write("org\tvariant\tn_perl\tn_pytd\tn_shared\tperl_only\tpytd_only\tidentical\t"
             "perl_runtime_sec\tpytd_runtime_sec\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
print(OUT, len(rows))
