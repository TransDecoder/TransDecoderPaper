#!/usr/bin/env python3
"""E33 step 5: compact tables for the write-up.

usage: python step5_12_tables_v1.py <step5 data dir in the repo>

Reads main_by_dataset.tsv, selected_by_dataset.tsv, hit_stats.tsv, full_scores.tsv there; writes
isolate_compact.tsv (per dataset: exact F1 Swiss-Prot / UniRef50 a5b24c3 / best UniRef50 top-n variant, per-type
deltas, rule firings, wrong trims / keeps, hit coverage, top hit >= 95% identity) and full_mode_compact.tsv
(per full-mode dataset: UniRef50 minus Swiss-Prot in exact F1, stop F1, stop sensitivity and precision, per-type
exact F1, predictions on negatives), and prints both.
"""
import sys

import pandas as pd

S = sys.argv[1]
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
M = pd.read_csv(f"{S}/main_by_dataset.tsv", sep="\t")
H = pd.read_csv(f"{S}/hit_stats.tsv", sep="\t")
SEL = pd.read_csv(f"{S}/selected_by_dataset.tsv", sep="\t")
best = SEL[SEL.selection == "best_uniref50"]
order = list(dict.fromkeys(M.dataset))
rows = []
for ds in order:
    m = M[M.dataset == ds].set_index("variant")
    sp, ur = m.loc["sprot|k1|top1|none"], m.loc["uniref50|k1|top1|none"]
    h = H[(H.dataset == ds) & (H.search == "k1")].set_index("db")
    b = best[best.dataset == ds].iloc[0]
    r = dict(dataset=ds, f1_release=m.loc["release", "f1_exact"], f1_sprot=sp.f1_exact, f1_uniref50=ur.f1_exact,
             d_uniref50=100 * (ur.f1_exact - sp.f1_exact), d_best_uniref50=b.d_f1_exact,
             d_perfect=100 * (m.loc["perfect", "f1_exact"] - sp.f1_exact))
    for t in TYPES:
        r[f"d_uniref50_{t}"] = 100 * (ur[f"f1_{t}"] - sp[f"f1_{t}"])
    for k in ("n_R1", "n_R2", "wrong_trim", "wrong_keep"):
        r[f"{k}_sprot"], r[f"{k}_uniref50"] = sp[k], ur[k]
    r["best_wrong_trim"], r["best_wrong_keep"] = b.wrong_trim, b.wrong_keep
    r["n_r3"] = sp.n_r3
    for db in ("sprot", "uniref50"):
        r[f"hit_{db}"] = h.loc[db, "frac_with_hit"]
        r[f"pid95_{db}"] = h.loc[db, "top_pid_ge95"]
    r["uniref50_top_cluster_own_member"] = h.loc["uniref50", "top_cluster_has_own_member"]
    r["uniref50_top_own_member_pid95"] = h.loc["uniref50", "top_own_member_pid_ge95"]
    rows.append(r)
I = pd.DataFrame(rows)
I.to_csv(f"{S}/isolate_compact.tsv", sep="\t", index=False, float_format="%.4f")

F = pd.read_csv(f"{S}/full_scores.tsv", sep="\t")
rows = []
for ds in dict.fromkeys(F.dataset):
    f = F[F.dataset == ds].set_index("arm")
    sp, ur, ab = f.loc["sprot"], f.loc["uniref50"], f.loc["abinitio"]
    r = dict(dataset=ds, f1_exact_abinitio=ab.f1_exact, f1_exact_sprot=sp.f1_exact, f1_exact_uniref50=ur.f1_exact,
             d_f1_exact=100 * (ur.f1_exact - sp.f1_exact), d_f1_stop=100 * (ur.f1_stop - sp.f1_stop),
             sens_stop_sprot=sp.sens_stop, sens_stop_uniref50=ur.sens_stop,
             prec_stop_sprot=sp.prec_stop, prec_stop_uniref50=ur.prec_stop,
             n_pred_sprot=sp.n_pred, n_pred_uniref50=ur.n_pred,
             neg_pred_abinitio=ab.n_pred_on_negatives, neg_pred_sprot=sp.n_pred_on_negatives,
             neg_pred_uniref50=ur.n_pred_on_negatives, n_negatives=sp.n_negatives,
             orfs_blast_hit_sprot=sp.orfs_blast_hit, orfs_blast_hit_uniref50=ur.orfs_blast_hit,
             R1_sprot=sp.R1, R1_uniref50=ur.R1, R2_sprot=sp.R2, R2_uniref50=ur.R2)
    for t in TYPES:
        r[f"d_{t}"] = 100 * (ur[f"f1_exact_{t}"] - sp[f"f1_exact_{t}"])
    rows.append(r)
G = pd.DataFrame(rows)
G.to_csv(f"{S}/full_mode_compact.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 300)
pd.set_option("display.max_columns", 60)
print(I.to_string(float_format=lambda z: f"{z:.3f}"))
print(G.to_string(float_format=lambda z: f"{z:.4f}"))
