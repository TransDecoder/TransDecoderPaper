#!/usr/bin/env python3
"""E36 step 2: exact F1 of every benchmark input for every alpha-form configuration on a fixed grid.

usage: python 05_tune_grid_v1.py <manifest.tsv> <labels_dir> <out.tsv.gz> <form> [<form> ...]

Benchmark rows only (48 ab initio inputs, 24 assist inputs). Forms (abl_common_v1.alpha_form; alpha clipped to [0, 5]):
  const   alpha = a0                               (the release form; a0 = 0.4 is V0)
  logL    alpha = a0 + b log2(L / 1000)            L = transcript length
  logO    alpha = a0 + b log2(O / 1000)            O = candidate (open) ORF length
  OL      alpha = a0 + b (O/L - 0.8)
  uL      alpha = a0 + b |u|/L
  stepOL  alpha = alo if O/L < r else ahi
  stepL   alpha = alo if L < l else ahi
Writes one row per (config, dataset): form, params, dataset, arm, group, org, f1_exact, tp_exact, n_changed.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abl_common_v1 as ac  # noqa: E402

A = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0]


def configs(form):
    if form == "const":
        return [dict(a0=a) for a in A]
    if form in ("logL", "logO"):
        return [dict(a0=a, b=round(b, 2)) for a in A for b in np.arange(-0.6, 1.2001, 0.05)]
    if form == "OL":
        return [dict(a0=a, b=round(b, 2)) for a in A for b in np.arange(-4, 4.001, 0.25)]
    if form == "uL":
        return [dict(a0=a, b=round(b, 2)) for a in A for b in np.arange(-6, 6.001, 0.5)]
    if form == "stepOL":
        return [dict(r=r, alo=lo, ahi=hi) for r in (0.5, 0.6, 0.7, 0.8, 0.9, 0.95) for lo in A for hi in A]
    if form == "stepL":
        return [dict(l=l, alo=lo, ahi=hi) for l in (500, 750, 1000, 1500, 2000, 3000, 4000) for lo in A for hi in A]
    raise ValueError(form)


def pstr(p):
    return ";".join(f"{k}={v:g}" for k, v in p.items())


def main():
    manifest, ldir, out = sys.argv[1:4]
    forms = sys.argv[4:]
    man = pd.read_csv(manifest, sep="\t")
    man = man[man.group.str.startswith("bench")]
    data = []
    for _, m in man.iterrows():
        lab, const = ac.load(ldir, m.dataset)
        data.append((m, lab, ac.FastScorer(lab, const)))
    rows = []
    for form in forms:
        cf = configs(form)
        for p in cf:
            ps = pstr(p)
            for m, lab, fs in data:
                d = ac.decide(lab, ac.alpha_form(form, p, lab), m.arm == "homology")
                r = fs.score(d, full=False)
                rows.append((form, ps, m.dataset, m.arm, m.group, m.org, r["f1_exact"], r["tp_exact"], r["n_changed"]))
        print(form, len(cf), "configs done", file=sys.stderr, flush=True)
    pd.DataFrame(rows, columns=["form", "params", "dataset", "arm", "group", "org", "f1_exact", "tp_exact",
                                "n_changed"]).to_csv(out, sep="\t", index=False, float_format="%.7f")


if __name__ == "__main__":
    main()
