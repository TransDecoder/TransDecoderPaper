#!/usr/bin/env python3
"""E35: check that V0 (release) offline scores equal the official scores of each run.

usage: python 03_check_official_v1.py <scores.tsv> <out.tsv>

Official tables (paper repo):
  data/start_classifier/adoption/bench_scores.tsv   benchmark complete / complete_fullutr / fragments, pytd_v7 and
                                                    pytd_v7h_homology (f1_stop, f1_exact, per-type f1_exact_*)
  data/benchmark/final_v7/mixed_scores.tsv          benchmark fragments_mixed_t1-t3, pytd_v7 (f1_stop, f1_exact; 4 dp)
  data/start_classifier/adoption/sim_scores.tsv     human d1/d8/d32 pytd_v7 and assist_on (f1, per type)
  data/simtrinity/depth/scores_by_type.tsv          human pytd_v7 per type with negatives (f1_exact_wneg; 4 dp)
  data/start_classifier/adoption/athal_depth_scores.tsv  athal A / B: f1_stop, f1_exact, f1_exact_wneg_* (4 dp)
Each row: dataset, metric, official, offline, abs diff; 'ok' when |diff| < 0.5 unit in the official's last decimal.
"""
import sys

import pandas as pd

T = "/broad/hptmp/bhaas/TRANSDECODER_PAPER/TransDecoderPaper"
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
S = pd.read_csv(sys.argv[1], sep="\t")
S = S[S.variant == "V0"].set_index("dataset")
rows = []


def add(ds, metric, off, val, dp):
    if pd.isna(off) and pd.isna(val):
        return
    rows.append(dict(dataset=ds, metric=metric, official=off, offline=round(val, 6), diff=abs(off - val),
                     ok=int(abs(off - val) < 0.5 * 10 ** -dp + 1e-9)))


b = pd.read_csv(f"{T}/data/start_classifier/adoption/bench_scores.tsv", sep="\t")
for _, r in b.iterrows():
    if r.tool not in ("pytd_v7", "pytd_v7h_homology") or r.variant not in ("complete", "complete_fullutr", "fragments"):
        continue
    ds = f"bench.{r.org}.{r.variant}" + (".H" if r.tool == "pytd_v7h_homology" else "")
    s = S.loc[ds]
    add(ds, "f1_stop", r.f1_stop, s.f1_stop, 5)
    add(ds, "f1_exact", r.f1_exact, s.f1_exact, 5)
    for t in TYPES:
        add(ds, f"f1_{t}", r[f"f1_exact_{t}"], s[f"f1_{t}"], 5)
mx = pd.read_csv(f"{T}/data/benchmark/final_v7/mixed_scores.tsv", sep="\t")
for _, r in mx[mx.tool == "pytd_v7"].iterrows():
    ds = f"bench.{r.org}.{r.variant}"
    if ds in S.index:
        add(ds, "f1_stop", r.f1_stop, S.loc[ds].f1_stop, 4)
        add(ds, "f1_exact", r.f1_exact, S.loc[ds].f1_exact, 4)
sm = pd.read_csv(f"{T}/data/start_classifier/adoption/sim_scores.tsv", sep="\t")
for _, r in sm.iterrows():
    if not r.dataset.startswith("sim.human") or r.arm not in ("pytd_v7", "assist_on"):
        continue
    ds = r.dataset + (".H" if r.arm == "assist_on" else "")
    s = S.loc[ds]
    add(ds, "f1_stop", r.f1_stop, s.f1_stop, 5)
    add(ds, "f1_exact", r.f1_exact, s.f1_exact, 5)
    for t in TYPES:
        add(ds, f"f1_{t}", r[f"f1_exact_{t}"], s[f"f1_{t}"], 5)
hd = pd.read_csv(f"{T}/data/simtrinity/depth/scores_by_type.tsv", sep="\t")
for _, r in hd[(hd.tool == "pytd_v7") & hd.level.isin(["d1", "d8", "d32"])].iterrows():
    ds = f"sim.human.{r.level}"
    add(ds, f"f1w_{r.type}", r.f1_exact_wneg, S.loc[ds][f"f1w_{r.type}"], 4)
ad = pd.read_csv(f"{T}/data/start_classifier/adoption/athal_depth_scores.tsv", sep="\t")
for _, r in ad.iterrows():
    if r.method not in ("TD_A_abinitio", "TD_B_homology_assist"):
        continue
    ds = f"sim.athal.{r.level}" + (".H" if r.method == "TD_B_homology_assist" else "")
    s = S.loc[ds]
    add(ds, "f1_stop", r.f1_stop, s.f1_stop, 4)
    add(ds, "f1_exact", r.f1_exact, s.f1_exact, 4)
    for t in TYPES:
        add(ds, f"f1w_{t}", r[f"f1_exact_wneg_{t}"], s[f"f1w_{t}"], 4)
out = pd.DataFrame(rows)
out.to_csv(sys.argv[2], sep="\t", index=False, float_format="%.6f")
print(f"{len(out)} comparisons over {out.dataset.nunique()} datasets; not ok: {int((out.ok == 0).sum())}")
print(out[out.ok == 0].to_string())
