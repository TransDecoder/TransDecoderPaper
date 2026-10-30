#!/usr/bin/env python3
"""E36 step 2: leave-one-organism-out selection of each alpha form on the benchmark (E13-type objective).

usage: python 06_select_v1.py <grid_dir> <out_dir>

Per arm, per organism, deltas in exact-F1 points vs V0 (const a0 = 0.4):
  ab initio: J_org = mean of [complete, complete_fullutr, fragments, mixed (mean of t1-t3)]   (50/50 complete vs
             fragmented, as E13; fragments and mixed share the fragmented half)
  assist:    J_org = mean of [complete, complete_fullutr, fragments, fragments]                (no mixed homology runs)
  E13 (strict, reported) = 0.5 (complete + fragments).
Selection on the training organisms: 'free' = argmax mean J; 'accept' = argmax mean J subject to the E35 acceptance
criteria on the training inputs (mean complete delta >= -0.1, no training input < -0.5). Ties -> fewer changed decisions.
Folds: each organism held out (7 train); 'all' = all 8 (the configuration carried to the held-out evaluation).
Writes loo_folds.tsv, loo_summary.tsv, selected.tsv, best_by_form_insample.tsv.
"""
import glob
import sys

import numpy as np
import pandas as pd

grid_dir, odir = sys.argv[1:3]
G = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{grid_dir}/*.tsv.gz"))], ignore_index=True)
G["input"] = G.dataset.str.split(".").str[2]
v0 = G[(G.form == "const") & (G.params == "a0=0.4")].set_index("dataset").f1_exact
G["d"] = 100 * (G.f1_exact - G.dataset.map(v0))
G["cls"] = G.input.str.replace(r"_t[123]$", "", regex=True)
orgs = sorted(G.org.unique())


def per_org(g, arm):
    t = g.pivot_table(index=["form", "params", "org"], columns="cls", values="d", aggfunc="mean")
    if arm == "abinitio":
        t["J"] = t[["complete", "complete_fullutr", "fragments", "fragments_mixed"]].mean(axis=1)
    else:
        t["J"] = 0.25 * t.complete + 0.25 * t.complete_fullutr + 0.5 * t.fragments
    t["E13"] = 0.5 * (t.complete + t.fragments)
    return t


def fold_stats(g, T, tr):
    """mean over training orgs of J etc. and acceptance on the training inputs, per config."""
    s = T[T.index.get_level_values("org").isin(tr)].groupby(level=["form", "params"]).mean()
    gg = g[g.org.isin(tr)]
    s["min_input"] = gg.groupby(["form", "params"]).d.min()
    s["changed"] = gg.groupby(["form", "params"]).n_changed.sum()
    s["accept"] = (s.complete >= -0.1) & (s.min_input >= -0.5)
    return s


def pick(s, mode):
    c = s if mode == "free" else s[s.accept]
    if c.empty:
        return None
    best = c.J.max()
    c = c[c.J >= best - 1e-12].sort_values("changed")
    return c.index[0]


folds, sel, ins = [], [], []
for arm in ["abinitio", "homology"]:
    g = G[G.arm == arm]
    T = per_org(g, arm)
    for form in sorted(g.form.unique()):
        gf, Tf = g[g.form == form], T[T.index.get_level_values("form") == form]
        best_in = fold_stats(gf, Tf, orgs)
        b = best_in.sort_values("J", ascending=False).head(1).reset_index()
        b.insert(0, "arm", arm)
        ins.append(b)
        for mode in ["free", "accept"]:
            for ho in orgs + ["all"]:
                tr = [o for o in orgs if o != ho]
                s = fold_stats(gf, Tf, tr)
                k = pick(s, mode)
                row = dict(arm=arm, form=form, mode=mode, held_out=ho)
                if k is None:
                    row["params"] = "none_acceptable"
                    folds.append(row)
                    continue
                row.update(params=k[1], train_J=s.loc[k].J)
                if ho != "all":
                    h = Tf.loc[(k[0], k[1], ho)]
                    hi = gf[(gf.params == k[1]) & (gf.org == ho)]
                    row.update(ho_J=h.J, ho_E13=h.E13, ho_complete=h.complete, ho_complete_fullutr=h.complete_fullutr,
                               ho_fragments=h.fragments, ho_mixed=h.get("fragments_mixed", np.nan),
                               ho_min_input=hi.d.min())
                    folds.append(row)
                else:
                    a = fold_stats(gf, Tf, orgs).loc[k]
                    gi = gf[gf.params == k[1]]
                    worst = gi.loc[gi.d.idxmin()]
                    sel.append(dict(arm=arm, form=form, mode=mode, params=k[1], J=a.J, E13=a.E13,
                                    complete=a.complete, complete_fullutr=a.complete_fullutr, fragments=a.fragments,
                                    mixed=a.get("fragments_mixed", np.nan), n_inputs=len(gi),
                                    n_below_05=int((gi.d < -0.5).sum()), worst_input=f"{worst.dataset} {worst.d:+.2f}",
                                    n_changed=int(gi.n_changed.sum()),
                                    acceptable=bool(a.complete >= -0.1 and gi.d.min() >= -0.5),
                                    both_up=bool(a.complete > 0 and a.fragments > 0 and
                                                 (arm == "homology" or a.fragments_mixed > 0))))
F = pd.DataFrame(folds)
F.to_csv(f"{odir}/loo_folds.tsv", sep="\t", index=False, float_format="%.3f")
L = F[F.held_out != "all"].groupby(["arm", "form", "mode"]).agg(
    loo_J=("ho_J", "mean"), loo_E13=("ho_E13", "mean"), loo_complete=("ho_complete", "mean"),
    loo_fragments=("ho_fragments", "mean"), loo_mixed=("ho_mixed", "mean"), folds_J_down=("ho_J", lambda x: int((x < 0).sum())),
    worst_fold_J=("ho_J", "min"), worst_ho_input=("ho_min_input", "min"),
    distinct_params=("params", "nunique")).reset_index()
L.to_csv(f"{odir}/loo_summary.tsv", sep="\t", index=False, float_format="%.3f")
pd.DataFrame(sel).to_csv(f"{odir}/selected.tsv", sep="\t", index=False, float_format="%.3f")
pd.concat(ins).to_csv(f"{odir}/best_by_form_insample.tsv", sep="\t", index=False, float_format="%.3f")
print(L.to_string())
print(pd.DataFrame(sel).to_string())
