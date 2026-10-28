#!/usr/bin/env python3
"""E33 step 3: can BLAST / Pfam hits of the 5'-open ORF decide keep-open vs trim-to-first-ATG?

usage: python 08_homology_rules.py <features_dir> <manifest.tsv> <out_dir>

manifest columns: dataset, homology_dir (contains blastp.outfmt6 and pfam.domtblout for the ORF ids of
the pytd_v7 run; benchmark: the pytd_v7_homology run, whose LongOrfs output is identical; simulated
assemblies: 07_homology_candidates_v1.sh), run_dir, truth, excluded (as 05_evaluate.py).

Per candidate (01_features.py rows; ORF peptide position 1 = first whole codon at the open 5' end;
first in-frame ATG at peptide position a = seg/3 + 1):
  BLAST: the first (best) line of the query in blastp.outfmt6 (DIAMOND -k 1); qstart, sstart.
  Pfam: smallest envelope start (env_from) over the ORF's domains (hmmsearch --cut_ga).
Rules (each a call on the candidates where it fires):
  R1_k[_s]   keep open if qstart <= a - k (k = 1, 5, 10, 20); _s: also sstart > 1
  R2_s_d     trim if sstart <= s and the subject N-terminus maps within d aa of the ATG:
             |qstart - (sstart - 1) - a| <= d  (s = 1, 3, 5, 10; d = 2, 5, 10)
  R3_k       keep open if a Pfam envelope starts at <= a - k (k = 1, 5)
Policies (decision columns; homology call where it fires, the run's release decision elsewhere):
  each rule alone, and the a-priori combination H = R2_3_5 (trim) + R1_5_s (keep) + R3_1 (keep);
  conflicting calls fall back to the release decision. 'oracle_on_H' = truth decision where H fires.
Writes rule_metrics.tsv (per dataset x rule: fires, precision and recall of its call, decisions changed
vs release and how many of those fix / break an exact ORF) and <out_dir>/decisions/<dataset>.* in the
format of 04_selftrain.py, for 05_evaluate.py.
"""
import os
import sys

import numpy as np
import pandas as pd

fdir, manifest, odir = sys.argv[1:4]
os.makedirs(f"{odir}/decisions", exist_ok=True)
man = pd.read_csv(manifest, sep="\t")

R1 = [(k, s) for k in (1, 5, 10, 20) for s in (False, True)]
R2 = [(s, d) for s in (1, 3, 5, 10) for d in (2, 5, 10)]
R3 = [1, 5]
rows = []
for _, m in man.iterrows():
    ds = m.dataset
    F = pd.read_csv(f"{fdir}/{ds}.tsv.gz", sep="\t", low_memory=False)
    b = pd.read_csv(f"{m.homology_dir}/blastp.outfmt6", sep="\t", header=None,
                    names=["q", "s", "pid", "alen", "mm", "gap", "qstart", "qend", "sstart", "send", "evalue", "bits"])
    b = b.drop_duplicates("q", keep="first").set_index("q")
    pf = {}
    for line in open(f"{m.homology_dir}/pfam.domtblout"):
        if line.startswith("#"):
            continue
        x = line.split()
        pf[x[0]] = min(pf.get(x[0], 10 ** 9), int(x[19]))
    F["a"] = F.seg // 3 + 1
    F["qstart"] = F.orf_id.map(b.qstart)
    F["sstart"] = F.orf_id.map(b.sstart)
    F["pf_env"] = F.orf_id.map(pf)
    hasb, hasp = F.qstart.notna().values, F.pf_env.notna().values
    qs, ss, a, pe = F.qstart.values, F.sstart.values, F.a.values, F.pf_env.values
    rel_keep = F.trim_run.values == 0
    lab = F.label.values
    rank = F.true_atg_rank.values
    trunc = lab == "truncated"
    utr = (lab == "utr") & (rank >= 1)
    utr1 = (lab == "utr") & (rank == 1)
    labelled = trunc | utr

    rules = {}
    with np.errstate(invalid="ignore"):
        for k, s in R1:
            f = hasb & (qs <= a - k) & ((ss > 1) if s else True)
            rules[f"R1_{k}" + ("_s" if s else "")] = (f, True)
        for s, d in R2:
            f = hasb & (ss <= s) & (np.abs(qs - (ss - 1) - a) <= d)
            rules[f"R2_{s}_{d}"] = (f, False)
        for k in R3:
            rules[f"R3_{k}"] = (hasp & (pe <= a - k), True)
    n_tr, n_ut = (trunc).sum(), (utr).sum()
    rows.append(dict(dataset=ds, rule="(coverage)", n_cand=len(F), n_truncated=n_tr, n_utr=n_ut,
                     blast_truncated=(hasb & trunc).sum() / max(n_tr, 1), blast_utr=(hasb & utr).sum() / max(n_ut, 1),
                     pfam_truncated=(hasp & trunc).sum() / max(n_tr, 1), pfam_utr=(hasp & utr).sum() / max(n_ut, 1)))
    for name, (f, keep) in rules.items():
        fl = f & labelled
        correct = (trunc if keep else utr) & fl
        changed = f & (rel_keep != keep)
        # exact effect of the changes: keep correct for truncated; trim exact only for utr rank 1
        fix = changed & ((trunc) if keep else utr1)
        brk = changed & ((utr1) if keep else trunc)
        rows.append(dict(dataset=ds, rule=name, call="keep" if keep else "trim", fires=f.sum(), fires_labelled=fl.sum(),
                         precision=correct.sum() / fl.sum() if fl.sum() else np.nan,
                         recall=correct.sum() / ((trunc if keep else utr).sum() or np.nan),
                         changed=changed.sum(), fixes_exact=fix.sum(), breaks_exact=brk.sum()))
    # decision columns
    H_keep = rules["R1_5_s"][0] | rules["R3_1"][0]
    H_trim = rules["R2_3_5"][0]
    out = pd.DataFrame({"orf_id": F.orf_id})
    for name, (f, keep) in rules.items():
        out[name] = np.where(f, 1.0 if keep else 0.0, rel_keep.astype(float))
    conflict = H_keep & H_trim
    out["H"] = np.where(conflict, rel_keep, np.where(H_trim, 0.0, np.where(H_keep, 1.0, rel_keep))).astype(float)
    out["H_noR3"] = np.where(rules["R1_5_s"][0] & H_trim, rel_keep,
                             np.where(H_trim, 0.0, np.where(rules["R1_5_s"][0], 1.0, rel_keep))).astype(float)
    fired = H_keep | H_trim
    truth_keep = np.where(trunc, 1.0, np.where(utr1, 0.0, rel_keep.astype(float)))
    out["oracle_on_H"] = np.where(fired, truth_keep, rel_keep.astype(float))
    rows.append(dict(dataset=ds, rule="H", call="mixed", fires=fired.sum(), fires_labelled=(fired & labelled).sum(),
                     changed=(fired & ((out.H.values >= .5) != rel_keep)).sum(), conflicts=conflict.sum()))
    out.to_csv(f"{odir}/decisions/{ds}.decisions.tsv.gz", sep="\t", index=False)
    with open(f"{odir}/decisions/{ds}.log.tsv", "w") as o:
        o.write(f"dataset\t{ds}\n")
    print(ds, "rules done", file=sys.stderr)

pd.DataFrame(rows).to_csv(f"{odir}/rule_metrics.tsv", sep="\t", index=False, float_format="%.4f")
