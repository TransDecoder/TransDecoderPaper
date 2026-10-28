#!/usr/bin/env python3
"""E33 / E32: Arabidopsis depth series with the PyTD homology start assist, next to E32's methods.

usage: python 19_athal_depth_summary.py <out.tsv> <level> [<level> ...]

Collects, with E32's collector (pipeline/simtrinity/13b_depth_collect.py, data root
$P/data/simtrinity_athal, read only):
  - E32's own runs (work root $P/work/simtrinity_athal) -> $P/work/start_classifier/athal_depth/collected/E32
  - each condition A / B / C (work root $P/work/start_classifier/athal_depth/<cond>) -> .../collected/<cond>
and writes one row per level x method: overall stop / exact F1 and per-type exact F1 with negatives
(f1_exact_wneg, as E32 reports). Methods: TD_A_abinitio, TD_B_homology_assist, TD_C_homology_noassist,
E32 pytd_v7 (release), longest_orf, and per type the best E32 method other than TransDecoder
(pytd_v7, perl_td_v600) and longest_orf ('best_other' rows give the method name per type).
Also checks that condition A reproduces E32's pytd_v7 row (overall and per type).
"""
import os
import subprocess
import sys

import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
T = f"{P}/TransDecoderPaper"
O = f"{P}/work/start_classifier/athal_depth"
out, levels = sys.argv[1], sys.argv[2:]
COL = f"{T}/pipeline/simtrinity/13b_depth_collect.py"
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]


def collect(work_root, dest):
    os.makedirs(dest, exist_ok=True)
    subprocess.run([sys.executable, COL, "--data-root", f"{P}/data/simtrinity_athal", "--work-root", work_root,
                    dest] + levels, check=True, capture_output=True)
    ov = pd.read_csv(f"{dest}/scores_overall.tsv", sep="\t")
    bt = pd.read_csv(f"{dest}/scores_by_type.tsv", sep="\t")
    return ov, bt


def rowset(ov, bt, tool, name):
    rows = []
    for lv in levels:
        o = ov[(ov.level == lv) & (ov.tool == tool)]
        if o.empty:
            continue
        r = dict(level=lv, method=name, f1_stop=o.f1_stop.iloc[0], f1_exact=o.f1_exact.iloc[0],
                 n_pred=o.n_pred.iloc[0], tp_exact=o.tp_exact.iloc[0])
        for t in TYPES:
            b = bt[(bt.level == lv) & (bt.tool == tool) & (bt.type == t)]
            r[f"f1_exact_wneg_{t}"] = b.f1_exact_wneg.iloc[0] if len(b) else float("nan")
        rows.append(r)
    return rows


e_ov, e_bt = collect(f"{P}/work/simtrinity_athal", f"{O}/collected/E32")
rows = rowset(e_ov, e_bt, "pytd_v7", "E32_pytd_v7_release") + rowset(e_ov, e_bt, "longest_orf", "longest_orf")
cond = {}
for c, name in (("A", "TD_A_abinitio"), ("B", "TD_B_homology_assist"), ("C", "TD_C_homology_noassist")):
    ov, bt = collect(f"{O}/{c}", f"{O}/collected/{c}")
    cond[c] = (ov, bt)
    rows += rowset(ov, bt, "pytd_v7", name)
# best other method (not TransDecoder, not longest ORF) per level: overall and per type
others = e_ov[~e_ov.tool.isin(["pytd_v7", "perl_td_v600", "longest_orf"])]
for lv in levels:
    o = others[others.level == lv]
    if o.empty:
        continue
    bo = o.loc[o.f1_exact.idxmax()]
    r = dict(level=lv, method="best_other", f1_stop=o.f1_stop.max(), f1_exact=bo.f1_exact,
             best_overall_exact=bo.tool)
    for t in TYPES:
        b = e_bt[(e_bt.level == lv) & (e_bt.type == t) & ~e_bt.tool.isin(["pytd_v7", "perl_td_v600", "longest_orf"])]
        if len(b):
            bb = b.loc[b.f1_exact_wneg.idxmax()]
            r[f"f1_exact_wneg_{t}"], r[f"best_{t}"] = bb.f1_exact_wneg, bb.tool
    rows.append(r)
res = pd.DataFrame(rows).sort_values(["level", "method"], key=lambda s: s.map(
    lambda v: int(v[1:]) if isinstance(v, str) and v.startswith("d") and v[1:].isdigit() else v))
res.to_csv(out, sep="\t", index=False, float_format="%.4f")

# reproduction check: A vs E32 pytd_v7
ov, bt = cond["A"]
for lv in levels:
    a = ov[(ov.level == lv) & (ov.tool == "pytd_v7")].drop(columns=["runtime_sec"], errors="ignore")
    e = e_ov[(e_ov.level == lv) & (e_ov.tool == "pytd_v7")].drop(columns=["runtime_sec"], errors="ignore")
    ab = bt[(bt.level == lv) & (bt.tool == "pytd_v7")].reset_index(drop=True)
    eb = e_bt[(e_bt.level == lv) & (e_bt.tool == "pytd_v7")].reset_index(drop=True)
    if a.empty or e.empty:
        continue
    same = a.reset_index(drop=True).equals(e.reset_index(drop=True)) and ab.equals(eb)
    print(f"{lv}: condition A reproduces E32 pytd_v7 (overall and per type): {same}")
pd.set_option("display.width", 250)
print(res.to_string(index=False))
