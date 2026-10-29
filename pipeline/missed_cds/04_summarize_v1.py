#!/usr/bin/env python3
"""E34: summary tables from 01_trace_missed_v1.py outputs (and the fresh runs of 02/03).

usage: python 04_summarize_v1.py <manifest.tsv[,manifest2.tsv]> <big_dir> <out_dir>
Writes to <out_dir>:
  classification.tsv          dataset x truth type: truth, found (stop), missed, miss classes a-e with sub-reasons
  b_characterisation.tsv      class (b) misses vs found CDSs of the same type (length, contrast, GC3, expression, homology, q)
  b_vs_found_by_length.tsv    the same contrast / fit measures within CDS-length bins
  b_by_stratum.tsv            class (b) by the candidate's decoy stratum (ineligible = in-frame score or contrast <= 0)
  c_winners.tsv               class (c): what won
  counterfactuals.tsv         offline selection variants + fresh runs (min length 30, --sensitive check)
  start5open_outcomes.tsv     3' partials / complete CDSs whose candidate is 5'-open: outcome counts
  start5open_kept_open.tsv    kept-open ones: scores, |u|, Kozak, first-ATG, reason
  start5open_by_u.tsv         trim-needed candidates (true start = first ATG) by |u| bin: kept-open rate and scores
"""
import os
import sys

import numpy as np
import pandas as pd

man_path, big, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
man = pd.concat([pd.read_csv(x, sep="\t") for x in man_path.split(",")], ignore_index=True)
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
LBINS = [(-1, 50, "<50"), (50, 100, "50-99"), (100, 150, "100-149"), (150, 200, "150-199"), (200, 300, "200-299"),
         (300, 10 ** 9, ">=300")]
UBINS = [(0, 30, "3-29"), (30, 60, "30-59"), (60, 150, "60-149"), (150, 300, "150-299"), (300, 10 ** 9, ">=300")]
ALPHA = 0.4


def lbin(n):
    return next(lab for lo, hi, lab in LBINS if lo <= n < hi)


def ubin(u):
    if u != u:
        return "no_atg"
    return next(lab for lo, hi, lab in UBINS if lo <= u < hi)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def gc3(s):
    c = s[2::3]
    return (c.count("G") + c.count("C")) / len(c)


def q(x, p):
    x = x.dropna()
    return float(np.quantile(x, p)) if len(x) else np.nan


cls_rows, bchar, blen, cwin, s5o, s5k, s5u, bstrat = [], [], [], [], [], [], [], []
cf = []
for _, m in man.iterrows():
    ds = m.dataset
    D = pd.read_csv(f"{big}/{ds}/per_cds.tsv.gz", sep="\t")
    S = pd.read_csv(f"{big}/{ds}/start5open.tsv.gz", sep="\t")
    tr = read_fasta(f"{m.run_dir}/input.fa.transdecoder_dir/__checkpoints_TDpredict/top_training_orfs.cds")
    tg = np.array([gc3(s) for s in tr.values()])
    t05, t50, t95 = np.quantile(tg, [0.05, 0.5, 0.95])
    D["gc3_outside_train"] = (D.gc3 < t05) | (D.gc3 > t95)
    D["hom_any"] = D.cand_blast.fillna(False).astype(bool) | D.cand_pfam.fillna(False).astype(bool)
    D["lbin"] = D.n_codons.map(lbin)
    # ---- classification ----
    for typ in TYPES + ["all"]:
        d = D if typ == "all" else D[D.type == typ]
        mm = d[d.status == "missed"]
        r = dict(dataset=ds, type=typ, n_truth=len(d), exact=(d.status == "exact").sum(),
                 stop_not_exact=(d.status == "stop").sum(), missed=len(mm),
                 missed_share_of_nonexact=round(len(mm) / max(1, len(d) - (d.status == "exact").sum()), 3))
        for c in "abcde":
            r[c] = (mm.miss_class == c).sum()
        for sub in ["s0<=0", "contrast<=0", "q>0.05"]:
            r[f"b_{sub}"] = ((mm.miss_class == "b") & (mm.miss_sub == sub)).sum()
        b = mm[mm.miss_class == "b"]
        r["b_no_pred_on_tx"] = b.win_id.isna().sum()
        r["b_other_orf_predicted"] = b.win_id.notna().sum()
        for sub in ["nonoverlapping_orf", "other_frame_overlap", "homology_ranked", "antisense"]:
            r[f"c_{sub}"] = ((mm.miss_class == "c") & (mm.miss_sub == sub)).sum()
        r["a_subs"] = ";".join(f"{k}:{v}" for k, v in mm[mm.miss_class == "a"].miss_sub.value_counts().items())
        cls_rows.append(r)
    # ---- class (b) vs found ----
    for typ in TYPES + ["all"]:
        d = D if typ == "all" else D[D.type == typ]
        for grp, g in (("missed_b", d[(d.status == "missed") & (d.miss_class == "b")]),
                       ("found", d[d.status != "missed"])):
            r = dict(dataset=ds, type=typ, group=grp, n=len(g))
            if not len(g):
                bchar.append(r)
                continue
            r.update(codons_median=g.n_codons.median(), codons_q25=q(g.n_codons, .25), codons_q75=q(g.n_codons, .75),
                     **{f"len_{lab}": round((g.lbin == lab).mean(), 3) for _, _, lab in LBINS},
                     cand_nt_median=g.cand_len.median(),
                     contrast_pb_median=round(g.cand_contrast_pb.median(), 4),
                     contrast_pb_q25=round(q(g.cand_contrast_pb, .25), 4),
                     true_contrast_pb_median=round(g.true_contrast_pb.median(), 4),
                     true_s0_pb_median=round(g.true_s0_pb.median(), 4),
                     true_s0_le0=round((g.true_s0_pb.dropna() <= 0).mean(), 3),
                     n_true_scored=g.true_contrast_pb.notna().sum(),
                     gc3_median=round(g.gc3.median(), 3), train_gc3_median=round(t50, 3),
                     train_gc3_5_95=f"{t05:.3f}-{t95:.3f}", gc3_outside_train=round(g.gc3_outside_train.mean(), 3),
                     tx_len_median=g.tx_len.median(),
                     **{f"tpm_{lab}": round((g.expr_bin == lab).mean(), 3) for lab in ["<1", "1-5", "5-20", "20-100", ">100"]},
                     hom_hit=round(g.hom_any.mean(), 3) if m.condition == "B" else np.nan)
            if grp == "missed_b":
                qq = g.cand_q
                r.update(ineligible=round((~g.cand_eligible.astype(bool)).mean(), 3),
                         q_le_0_1=round((qq <= 0.1).mean(), 3), q_le_0_2=round((qq <= 0.2).mean(), 3),
                         q_le_0_5=round((qq <= 0.5).mean(), 3), q_median=round(qq.median(), 3),
                         accepted_at_q0_2=int(g.cand_accept_q02.sum()), accepted_at_q0_5=int(g.cand_accept_q05.sum()))
            bchar.append(r)
        # within length bins
        for _, _, lab in LBINS:
            for grp, g in (("missed_b", d[(d.status == "missed") & (d.miss_class == "b") & (d.lbin == lab)]),
                           ("found", d[(d.status != "missed") & (d.lbin == lab)])):
                blen.append(dict(dataset=ds, type=typ, codons=lab, group=grp, n=len(g),
                                 contrast_pb_median=round(g.cand_contrast_pb.median(), 4) if len(g) else np.nan,
                                 true_contrast_pb_median=round(g.true_contrast_pb.median(), 4) if len(g) else np.nan,
                                 true_s0_pb_median=round(g.true_s0_pb.median(), 4) if len(g) else np.nan,
                                 gc3_outside_train=round(g.gc3_outside_train.mean(), 3) if len(g) else np.nan,
                                 tpm_lt1=round((g.expr_bin == "<1").mean(), 3) if len(g) else np.nan,
                                 found_rate=np.nan))
            nb = ((d.status == "missed") & (d.miss_class == "b") & (d.lbin == lab)).sum()
            nall = (d.lbin == lab).sum()
            blen[-1]["found_rate"] = blen[-2]["found_rate"] = round(1 - (
                (d.status == "missed") & (d.lbin == lab)).sum() / nall, 3) if nall else np.nan
    # ---- class (b) by decoy stratum (ORF-length bin, transcript-length bin, complete ORF) ----
    bb = D[(D.status == "missed") & (D.miss_class == "b")]
    fd = D[D.status != "missed"]
    for k, g in bb.groupby(bb.cand_stratum.fillna("ineligible")):
        f = fd[fd.cand_stratum == k]
        bstrat.append(dict(dataset=ds, stratum=k, n_missed_b=len(g), n_found=len(f),
                           missed_b_q_median=round(g.cand_q.median(), 3),
                           missed_b_contrast_median=round(g.cand_contrast.median(), 2),
                           found_contrast_median=round(f.cand_contrast.median(), 2) if len(f) else np.nan,
                           missed_b_accept_q02=int(g.cand_accept_q02.sum()),
                           missed_b_accept_q05=int(g.cand_accept_q05.sum())))
    # ---- class (c) winners ----
    c = D[(D.status == "missed") & (D.miss_class == "c")]
    for typ, g in c.groupby("type"):
        cwin.append(dict(dataset=ds, type=typ, n=len(g), **g.miss_sub.value_counts().to_dict(),
                         winner_longer=int(g.win_longer.sum()), winner_higher_contrast=int(g.win_higher_contrast.sum()),
                         winner_same_frame=int(g.win_same_frame.sum()),
                         winner_more_homology=int((g.win_hom > g.cand_hom).sum()),
                         cand_codons_median=g.n_codons.median(), winner_nt_median=g.win_len.median()))
    # ---- counterfactuals ----
    C = pd.read_csv(f"{big}/{ds}/counterfactuals.tsv", sep="\t")
    C.insert(0, "dataset", ds)
    C["source"] = "offline (run's intermediates + PyTD refinement code)"
    cf.append(C)
    # ---- start decision on 5'-open candidates ----
    if m.condition == "B":
        H = pd.read_csv(f"{big}/{ds}/homology_decisions.tsv.gz", sep="\t", dtype=str).set_index("orf_id")
        S["decision_source"] = S.cand_id.map(H.source)
        S["rules_fired"] = S.cand_id.map(H.rules_fired)
    else:
        S["decision_source"] = "ab_initio"
        S["rules_fired"] = "."
    S["has_atg"] = S.atg_first.notna()
    for typ in ["3prime_partial", "complete"]:
        s = S[S.type == typ]
        r = dict(dataset=ds, type=typ, n_5open_candidates=len(s), no_inframe_atg=(~s.has_atg).sum(),
                 true_is_first_atg=int(s.true_is_first_atg.fillna(False).astype(bool).sum()),
                 true_is_later_atg=int((s.has_atg & ~s.true_is_first_atg.fillna(False).astype(bool)).sum()))
        for o in ["exact", "kept_open", "trimmed_other_atg", "missed"]:
            r[o] = (s.outcome == o).sum()
        s5o.append(r)
        k = s[s.outcome == "kept_open"].copy()
        first = k.true_is_first_atg.fillna(False).astype(bool)
        k["reason"] = np.where(k.atg_first.isna(), "no_inframe_atg", np.where(~first, "true_start_not_first_atg",
                               np.where(k.decision_source == "homology", "homology_keep_R1_R3",
                                        np.where(k.coding_u >= ALPHA * k.atg_first, "u_scores_coding(>=0.4|u|)",
                                                 "kozak_veto(coding<0.4|u|,K<0)"))))
        for reason, g in list(k.groupby("reason")) + [("all", k)]:
            s5k.append(dict(dataset=ds, type=typ, reason=reason, n=len(g),
                            u_median=g.atg_first.median(), u_q25=q(g.atg_first, .25), u_q75=q(g.atg_first, .75),
                            coding_median=g.coding_u.median(), utr_median=g.utr_u.median(),
                            coding_per_u_median=round((g.coding_u / g.atg_first).median(), 3),
                            kozak_median=g.kozak.median(), kozak_lt0=round((g.kozak < 0).mean(), 3),
                            **{f"u_{lab}": int((g.atg_first.map(ubin) == lab).sum()) for _, _, lab in UBINS}))
        # trim-needed (true start = first ATG, stop found): kept-open rate by |u|
        t = s[s.true_is_first_atg.fillna(False).astype(bool) & (s.outcome != "missed")]
        for _, _, lab in UBINS + [(0, 0, "all")]:
            g = t if lab == "all" else t[t.atg_first.map(ubin) == lab]
            if not len(g):
                continue
            s5u.append(dict(dataset=ds, type=typ, u_bin=lab, n=len(g),
                            kept_open=int((g.outcome == "kept_open").sum()),
                            kept_open_rate=round((g.outcome == "kept_open").mean(), 3),
                            coding_per_u_median=round((g.coding_u / g.atg_first).median(), 3),
                            coding_ge_alpha_u=round((g.coding_u >= ALPHA * g.atg_first).mean(), 3),
                            kozak_median=round(g.kozak.median(), 2),
                            ab_initio_keep=int((g.ab_initio == "keep").sum())))

# fresh runs: min length 30 (the --sensitive program check is 07_check_sensitive_v1.py)
for mode in ("minlen30",):
    for _, m in man.iterrows():
        f = f"{big}/{mode}/{m.dataset}/counterfactuals.tsv"
        if os.path.exists(f):
            C = pd.read_csv(f, sep="\t").iloc[[0]].copy()  # the run's own release-q predictions
            base = pd.read_csv(f"{big}/{m.dataset}/per_cds.tsv.gz", sep="\t", usecols=["transcript_id", "status"])
            new = pd.read_csv(f"{big}/{mode}/{m.dataset}/per_cds.tsv.gz", sep="\t", usecols=["transcript_id", "status"])
            j = base.merge(new, on="transcript_id", suffixes=("_rel", "_new"))
            C["recovered_stop"] = ((j.status_rel == "missed") & (j.status_new != "missed")).sum()
            C["recovered_exact"] = ((j.status_rel == "missed") & (j.status_new == "exact")).sum()
            C["lost_stop"] = ((j.status_rel != "missed") & (j.status_new == "missed")).sum()
            C["lost_exact"] = ((j.status_rel == "exact") & (j.status_new != "exact")).sum()
            C["variant"] = {"minlen30": "LongOrfs -m 30 (fresh run)", "sensitive": "--sensitive (fresh Predict, check)"}[mode]
            C.insert(0, "dataset", m.dataset)
            C["source"] = "fresh PyTD a5b24c3 run"
            cf.append(C)

CF = pd.concat(cf, ignore_index=True)
rel = CF[CF.variant == "q0.05 (release)"].set_index("dataset")
CF["d_f1_exact_pts"] = [round(100 * (r.f1_exact - rel.loc[r.dataset, "f1_exact"]), 2) for r in CF.itertuples()]
CF["d_f1_stop_pts"] = [round(100 * (r.f1_stop - rel.loc[r.dataset, "f1_stop"]), 2) for r in CF.itertuples()]
CF["neg_with_pred_pct"] = (100 * CF.n_pred_neg / CF.n_neg).round(2)
CF["extra_pred_on_neg"] = [r.n_pred_neg - rel.loc[r.dataset, "n_pred_neg"] for r in CF.itertuples()]
CF["missed_release"] = [rel.loc[r.dataset, "n_truth"] - rel.loc[r.dataset, "tp_stop"] for r in CF.itertuples()]
for c in ("f1_stop", "f1_exact"):
    CF[c] = CF[c].round(4)

pd.DataFrame(cls_rows).to_csv(f"{out}/classification.tsv", sep="\t", index=False)
pd.DataFrame(bchar).to_csv(f"{out}/b_characterisation.tsv", sep="\t", index=False)
pd.DataFrame(blen).to_csv(f"{out}/b_vs_found_by_length.tsv", sep="\t", index=False)
pd.DataFrame(bstrat).to_csv(f"{out}/b_by_stratum.tsv", sep="\t", index=False)
pd.DataFrame(cwin).fillna(0).to_csv(f"{out}/c_winners.tsv", sep="\t", index=False)
CF.to_csv(f"{out}/counterfactuals.tsv", sep="\t", index=False)
pd.DataFrame(s5o).to_csv(f"{out}/start5open_outcomes.tsv", sep="\t", index=False)
pd.DataFrame(s5k).to_csv(f"{out}/start5open_kept_open.tsv", sep="\t", index=False)
pd.DataFrame(s5u).to_csv(f"{out}/start5open_by_u.tsv", sep="\t", index=False)
print("written", out)
