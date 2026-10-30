#!/usr/bin/env python3
"""E36 step 4: oracle-alpha ceilings, chosen on each evaluation dataset itself.

usage: python 08_ceiling_v1.py <manifest.tsv> <labels_dir> <out_dir>

For every dataset (both arms; with the assist the homology call still overrides where it fires):
  const_oracle   best single alpha from AG for this dataset
  L_oracle       best alpha per transcript-length bin (5 bins)
  OL_oracle      best alpha per O/L bin (5 bins)
  LxOL_oracle    best alpha per (L, O/L) bin (25 bins)   <- the task's ceiling for any length-based rule
  cand_oracle    every candidate the alpha rule decides (not overridden) decided right (no rule can beat this)
AG = {0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 1.5, 2, 3}; L bins [0,500,1000,2000,4000,inf) nt; O/L bins
[0,0.5,0.7,0.9,0.95,1]. Exact tp is additive over candidates (no transcript has two candidates; labels.log), so the
best alpha per bin maximises that bin's exact count. Reported: delta tp and delta overall exact F1 (points) vs V0;
GM12878: delta of the E29 exact-given-stop shares (5' partial, complete).
Writes ceiling.tsv and ceiling_bins.tsv (chosen alpha per bin for LxOL).
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abl_common_v1 as ac  # noqa: E402

AG = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0])
LB = [0, 500, 1000, 2000, 4000, np.inf]
OB = [0, 0.5, 0.7, 0.9, 0.95, 1.0001]


def main():
    manifest, ldir, odir = sys.argv[1:4]
    man = pd.read_csv(manifest, sep="\t")
    rows, brows = [], []
    for _, m in man.iterrows():
        lab, const = ac.load(ldir, m.dataset)
        hom = m.arm == "homology"
        inc = (lab.on_excl.values == 0) & (lab.truth_type.values != "-")
        oth = lab.other_exact.values.astype(bool)
        ek, et = lab.e_keep.values.astype(bool), lab.e_trim.values.astype(bool)
        run = lab.trim_run.values.astype(bool)
        e0 = (oth | np.where(run, et, ek)) & inc
        tt = lab.truth_type.values
        lb = np.digitize(lab.tx_len.values, LB[1:-1])
        ob = np.digitize(lab.orf_len.values / lab.tx_len.values, OB[1:-1])
        E = np.array([((oth | np.where(ac.decide(lab, a, hom), et, ek)) & inc) for a in AG])  # n_alpha x n_cand
        D = E.astype(int) - e0.astype(int)
        den = const["n_truth"] + const["n_pred"]
        free = ~np.isin(lab.hcall.values, ["trim", "keep"]) if hom else np.ones(len(lab), bool)
        best_c = (oth | ek | et) & inc
        res = {}
        res["const_oracle"] = (D.sum(1).max(), None)
        for name, key in [("L_oracle", lb), ("OL_oracle", ob), ("LxOL_oracle", lb * 10 + ob)]:
            sel = np.zeros(len(lab), int)
            for k in np.unique(key):
                s = key == k
                tot = D[:, s].sum(1)
                cand = np.where(tot == tot.max())[0]
                j = cand[np.argmin(np.abs(AG[cand] - 0.4))]
                sel[s] = j
                if name == "LxOL_oracle":
                    brows.append(dict(dataset=m.dataset, arm=m.arm, Lbin=f"[{LB[k // 10]},{LB[k // 10 + 1]})",
                                      OLbin=f"[{OB[k % 10]},{min(OB[k % 10 + 1], 1)})", n_cand=int(s.sum()),
                                      n_scored=int((s & inc).sum()), best_alpha=AG[j], d_tp=int(tot[j])))
            res[name] = (int(D[sel, np.arange(len(lab))].sum()), sel)
        ce = np.where(free, best_c, e0)
        res["cand_oracle"] = (int(ce.sum() - e0.sum()), None)
        ja = int(np.argmax(D.sum(1)))
        for name, (dtp, sel) in res.items():
            r = dict(dataset=m.dataset, arm=m.arm, group=m.group, org=m.org, input=m.input, ceiling=name,
                     d_tp=int(dtp), d_f1_exact=200 * dtp / den, f1_exact_v0=2 * const["tp_exact"] / den,
                     best_const_alpha=AG[ja] if name == "const_oracle" else np.nan)
            if m.group == "gm12878":
                if name == "cand_oracle":
                    enew = ce
                elif name == "const_oracle":
                    enew = E[ja]
                else:
                    enew = E[sel, np.arange(len(lab))]
                for t in ("5prime_partial", "complete"):
                    s = tt == t
                    r[f"d_e29_{t}"] = 100 * (int(enew[s].sum()) - int(e0[s].sum())) / const[f"tpstop_{t}"]
            rows.append(r)
        print(m.dataset, "done", file=sys.stderr, flush=True)
    pd.DataFrame(rows).to_csv(f"{odir}/ceiling.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(brows).to_csv(f"{odir}/ceiling_bins.tsv", sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
