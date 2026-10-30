#!/usr/bin/env python3
"""E33 / E31: human (simulated GM12878) depth series with the PyTD homology start assist, next to E30/E31's methods.
Adapted from 19_athal_depth_summary.py (Arabidopsis, E32).

usage: python 22_human_depth_summary.py <out.tsv> <level> [<level> ...]

- E31's methods: the committed E31 tables data/simtrinity/depth/{scores_overall,scores_by_type,negatives}.tsv (read only).
- Conditions A / B / C (work root $P/work/start_classifier/human_depth/<cond>) are collected with E31's collector
  (pipeline/simtrinity/13b_depth_collect.py) into .../human_depth/collected/<cond>. The collector's --data-root mode
  reads every level from <root>/depth/<level>; for the human layout (d1 = E30 at $P/data/simtrinity, dN at
  $P/data/simtrinity/depth/dN) a symlink root is made at .../human_depth/dataroot (d1/labels -> E30 labels,
  d1/Trinity.fasta -> E30 assembly, dN -> E31 level dirs). Outputs are as in the default (E31) mode.
Writes one row per level x method: overall stop / exact F1, per-type exact F1 with negatives (f1_exact_wneg), and
negatives (lncRNA-only contigs) with a prediction. Methods: TD_A_abinitio, TD_B_homology_assist, TD_C_homology_noassist,
E31_pytd_v7_release, longest_orf, best_other (best E31 method other than TransDecoder (pytd_v7, perl_td_v600) and
longest_orf; overall by exact F1 and per type, the method names in best_*).
Checks that condition A reproduces E31's committed pytd_v7 rows (overall, per type, negatives) and that A's preds.tsv
and longest_orf.preds.tsv are byte-identical to the release runs ($P/work/simtrinity/{runs,depth/<level>/runs}/pytd_v7).
"""
import filecmp
import os
import subprocess
import sys

import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
T = f"{P}/TransDecoderPaper"
O = f"{P}/work/start_classifier/human_depth"
out, levels = sys.argv[1], sys.argv[2:]
COL = f"{T}/pipeline/simtrinity/13b_depth_collect.py"
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
DR = f"{O}/dataroot"


def shim():
    os.makedirs(f"{DR}/depth/d1", exist_ok=True)
    for name, target in (("labels", f"{P}/data/simtrinity/labels"), ("Trinity.fasta", f"{P}/data/simtrinity/trinity/Trinity.fasta")):
        if not os.path.lexists(f"{DR}/depth/d1/{name}"):
            os.symlink(target, f"{DR}/depth/d1/{name}")
    for lv in levels:
        if lv != "d1" and not os.path.lexists(f"{DR}/depth/{lv}"):
            os.symlink(f"{P}/data/simtrinity/depth/{lv}", f"{DR}/depth/{lv}")


def collect(work_root, dest):
    os.makedirs(dest, exist_ok=True)
    subprocess.run([sys.executable, COL, "--data-root", DR, "--work-root", work_root, dest] + levels,
                   check=True, capture_output=True)
    return tuple(pd.read_csv(f"{dest}/{n}.tsv", sep="\t") for n in ("scores_overall", "scores_by_type", "negatives"))


def rowset(ov, bt, ng, tool, name):
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
        n = ng[(ng.level == lv) & (ng.tool == tool)]
        if len(n):
            r["n_negatives"], r["negatives_with_pred"] = n.n_negatives.iloc[0], n.negatives_with_pred.iloc[0]
            r["frac_negatives_with_pred"] = n.frac_negatives_with_pred.iloc[0]
        rows.append(r)
    return rows


shim()
E = f"{T}/data/simtrinity/depth"
e_ov, e_bt, e_ng = (pd.read_csv(f"{E}/{n}.tsv", sep="\t") for n in ("scores_overall", "scores_by_type", "negatives"))
rows = rowset(e_ov, e_bt, e_ng, "pytd_v7", "E31_pytd_v7_release") + rowset(e_ov, e_bt, e_ng, "longest_orf", "longest_orf")
cond = {}
for c, name in (("A", "TD_A_abinitio"), ("B", "TD_B_homology_assist"), ("C", "TD_C_homology_noassist")):
    cond[c] = collect(f"{O}/{c}", f"{O}/collected/{c}")
    rows += rowset(*cond[c], "pytd_v7", name)
EXCL = ["pytd_v7", "perl_td_v600", "longest_orf"]
others = e_ov[~e_ov.tool.isin(EXCL)]
for lv in levels:
    o = others[others.level == lv]
    if o.empty:
        continue
    bo = o.loc[o.f1_exact.idxmax()]
    r = dict(level=lv, method="best_other", f1_stop=o.f1_stop.max(), f1_exact=bo.f1_exact, best_overall_exact=bo.tool)
    for t in TYPES:
        b = e_bt[(e_bt.level == lv) & (e_bt.type == t) & ~e_bt.tool.isin(EXCL)]
        if len(b):
            bb = b.loc[b.f1_exact_wneg.idxmax()]
            r[f"f1_exact_wneg_{t}"], r[f"best_{t}"] = bb.f1_exact_wneg, bb.tool
    rows.append(r)
res = pd.DataFrame(rows)
res["_l"] = res.level.str[1:].astype(int)
res = res.sort_values(["_l", "method"]).drop(columns="_l")
res.to_csv(out, sep="\t", index=False, float_format="%.4f")

# reproduction check: A vs E31's committed pytd_v7 rows, and A's predictions vs the release runs
ov, bt, ng = cond["A"]
for lv in levels:
    key = lambda d: d[(d.level == lv) & (d.tool == "pytd_v7")].drop(columns=["runtime_sec"], errors="ignore").reset_index(drop=True)
    a, e = key(ov), key(e_ov)
    if a.empty or e.empty:
        print(f"{lv}: missing rows")
        continue
    same = (a.equals(e.drop_duplicates().reset_index(drop=True)) and key(bt).equals(key(e_bt))
            and key(ng).equals(key(e_ng).drop_duplicates().reset_index(drop=True)))
    rel = f"{P}/work/simtrinity/runs/pytd_v7" if lv == "d1" else f"{P}/work/simtrinity/depth/{lv}/runs/pytd_v7"
    mine = f"{O}/A/depth/{lv}/runs/pytd_v7"
    ident = {f: (filecmp.cmp(f"{mine}/{f}", f"{rel}/{f}", shallow=False) if os.path.exists(f"{rel}/{f}") else "absent")
             for f in ("preds.tsv", "longest_orf.preds.tsv", "input.fa.transdecoder.gff3")}
    print(f"{lv}: condition A reproduces E31 committed pytd_v7 rows (overall, per type, negatives): {same}; "
          f"byte-identical to release run ({open(f'{rel}/pytd_version.txt').read().strip()}): {ident}")
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
print(res.to_string(index=False))
