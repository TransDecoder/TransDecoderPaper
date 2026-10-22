#!/usr/bin/env python3
"""Kozak score distributions on the benchmark, in bits.

Models: start_refinement.kozak.tsv of the pytd_dev_v2 runs (complete inputs; written in nats by PyTD
48c870c, converted here to bits). Per organism:
  true_start     annotated start codons of complete-type truth CDSs (complete inputs)
  internal_ATG   in-frame ATGs inside those CDSs
  frag_firstATG  first in-frame ATG of the selected 5'-open ORFs (fragments inputs, pytd_dev_v2
                 best_candidates.gff3), with the upstream length |u| in nt
Quantiles p10/p25/p50/p75/p90 and the share of sites above -2, -1, 0 and +1 bits.

usage: kozak_score_distribution.py > results/kozak_score_distribution.tsv
"""
import math
import sys

sys.path.insert(0, "/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/trees/PY_dev")
from pytransdecoder.core import kozak  # noqa: E402

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
RUNS, INPUTS = f"{P}/work/bench_v3/runs", f"{P}/data/bench_v3/inputs"
ORGS = "athal creinhardtii dmel human mouse pfalciparum spombe tthermophila".split()
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
LN2 = math.log(2)


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def fasta(path):
    d, k, b = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if k:
                d[k] = "".join(b)
            k, b = line[1:].split()[0], []
        else:
            b.append(line.strip().upper())
    if k:
        d[k] = "".join(b)
    return d


def quantiles(x):
    x = sorted(x)
    return [x[int(p * (len(x) - 1))] for p in (.1, .25, .5, .75, .9)]


def row(org, name, x, unit="bits"):
    q = "\t".join(f"{v:.2f}" for v in quantiles(x))
    share = "\t".join(f"{sum(v > t for v in x) / len(x):.2f}" for t in (-2, -1, 0, 1)) if unit == "bits" \
        else "\t".join("NA" for _ in range(4))
    print(f"{org}\t{name}\t{unit}\t{len(x)}\t{q}\t{share}")


print("org\tset\tunit\tn\tp10\tp25\tp50\tp75\tp90\tabove_-2\tabove_-1\tabove_0\tabove_+1")
for org in ORGS:
    model = kozak.read_model(f"{RUNS}/{org}/complete/pytd_dev_v2/input.fa.transdecoder_dir/start_refinement.kozak.tsv")
    model = [{b: v / LN2 for b, v in pos.items()} for pos in model]   # nats -> bits
    tx = fasta(f"{INPUTS}/{org}/complete.fa")
    true_k, int_k = [], []
    for line in open(f"{INPUTS}/{org}/complete.truth.tsv"):
        f = line.split("\t")
        if f[0] == "transcript_id" or f[7] != "complete" or f[0] not in tx:
            continue
        s, a, b = tx[f[0]], int(f[1]), int(f[2])
        if f[3] == "-":
            s = rc(s)
            a, b = len(s) - b + 1, len(s) - a + 1
        i = a - 1
        if s[i:i + 3] != "ATG":
            continue
        true_k.append(kozak.score(model, s, i))
        int_k += [kozak.score(model, s, j) for j in range(i + 3, b - 3, 3) if s[j:j + 3] == "ATG"]
    ft = fasta(f"{INPUTS}/{org}/fragments.fa")
    ulen, uk = [], []
    for line in open(f"{RUNS}/{org}/fragments/pytd_dev_v2/input.fa.transdecoder_dir/best_candidates.gff3"):
        f = line.split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        s, a, b = ft[f[0]], int(f[3]), int(f[4])
        if f[6] == "-":
            s = rc(s)
            a, b = len(s) - b + 1, len(s) - a + 1
        orf = s[a - 1:b]
        if orf[:3] == "ATG":
            continue
        atg = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"), None)
        if atg is None:
            continue
        ulen.append(atg)
        uk.append(kozak.score(model, s, a - 1 + atg))
    row(org, "true_start", true_k)
    row(org, "internal_ATG", int_k)
    row(org, "frag_firstATG", uk)
    row(org, "frag_upstream_len", ulen, unit="nt")
