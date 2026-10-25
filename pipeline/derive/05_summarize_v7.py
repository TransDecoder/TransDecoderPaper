"""Summary tables from <out>/per_cds/*.per_cds.tsv.gz written by 05_long_orf_failures_v7.py (release
rule). Copy of analyses/long_orf_failures/summarize.py with these changes:
  - the coding test is S1 > 0, C > 0 and decoy q <= 0.05; there is no short-ORF guard, so
    fail_score (used for every "why" table) is every coding-test failure, including q > 0.05;
  - columns fail_short_guard -> fail_decoy_q; pass_coding_cds -> pass_scores_cds (S1 > 0 and C > 0
    on the bare CDS, which has no q-value), compared with pass_scores of the matched ORF.
Called by `05_long_orf_failures_v7.py summarize`, or: 05_summarize_v7.py <out>."""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact

ORGS = ["athal", "creinhardtii", "dmel", "human", "mouse", "pfalciparum", "spombe", "tthermophila"]
LEN_ORDER = ["0-299", "300-599", "600-999", "1000-1999", "2000-4999", ">=5000"]
MODELS = ["all", "rand500", "long500", "gc3match500"]
FEATURES = ["cds_len", "orf_len", "ext5", "utr5_len", "utr3_len", "gc", "gc3", "gc3_dev_train",
            "codon_llr_train", "s1_per_base", "contrast_rel", "aa_entropy", "top1_frac", "top3_frac",
            "seg_lowc_frac", "tandem3_frac", "disorder_frac", "wmer_unique_pct", "altframe_overlap_frac",
            "has_pfam", "has_sprot", "sprot_pident", "organelle",
            "aa_S", "aa_G", "aa_P", "aa_E", "aa_Q", "aa_N", "aa_K", "aa_C", "aa_H", "aa_L", "aa_I", "aa_A"]
KEYWORDS = ["hypothetical", "uncharacterized", "unknown", "putative", "conserved", "ribosomal",
            "histone", "zinc finger", "transposon", "transposable", "retrotransposon", "mucin", "collagen",
            "keratin", "defensin", "erythrocyte membrane protein", "rifin", "stevor", "PfEMP1",
            "surface", "repeat", "kinase", "cysteine-rich", "glycine-rich", "proline-rich",
            "serine-rich", "low complexity", "chloroplast", "mitochondrial", "pseudo", "LCR", "SCR",
            "F-box", "PPR", "pentatricopeptide", "leucine-rich", "meiotically", "wtf"]


def load(out):
    frames, metas = [], []
    for org in ORGS:
        p = out / "per_cds" / f"{org}.per_cds.tsv.gz"
        if not p.exists():
            continue
        d = pd.read_csv(p, sep="\t", low_memory=False, keep_default_na=False,
                        na_values=["NA", ""])
        frames.append(d)
        metas.append(pd.read_csv(out / "per_cds" / f"{org}.meta.tsv", sep="\t"))
    d = pd.concat(frames, ignore_index=True)
    meta = pd.concat(metas, ignore_index=True).set_index("org")
    d["gc3_dev_train"] = d["gc3"] - d["org"].map(meta["training_gc3_mean"])
    d["matched"] = d["category"] != "no_orf"
    d["fail"] = d["category"].isin(["fail_coding_long_rescued", "fail_coding_lost"])
    # release rule: every coding-test failure (S1 <= 0, a shifted frame scores higher, or q > 0.05)
    d["fail_score"] = d["fail"]
    d["fail_frame"] = d["fail"] & d["fail_reason"].isin(["S1<=0", "shifted_frame>=S1"])
    return d, meta


def pct(a, b):
    return round(100.0 * a / b, 2) if b else np.nan


def by_org(d, meta):
    rows = []
    for org, g in list(d.groupby("org")) + [("ALL", d)]:
        m = g[g.matched]
        fail = m[m.fail]
        r = dict(org=org, n_truth=len(g), no_orf=int((~g.matched).sum()), n_matched=len(m),
                 L_star=meta.loc[org, "L_star"] if org in meta.index else np.nan,
                 pass_coding=int((m.category == "pass_coding").sum()),
                 fail_coding=len(fail), fail_coding_pct=pct(len(fail), len(m)),
                 fail_S1_le0=int((fail.fail_reason == "S1<=0").sum()),
                 fail_shifted_frame=int((fail.fail_reason == "shifted_frame>=S1").sum()),
                 fail_decoy_q=int((fail.fail_reason == "decoy_q>0.05").sum()),
                 fail_frame_pct=pct(int(m.fail_frame.sum()), len(m)),
                 fail_ge300=int((fail.orf_len >= 300).sum()),
                 fail_ge300_pct_of_ge300=pct(int((fail.orf_len >= 300).sum()), int((m.orf_len >= 300).sum())),
                 rescued_by_long_rule=int((fail.category == "fail_coding_long_rescued").sum()),
                 rescued_and_predicted=int(((fail.category == "fail_coding_long_rescued") &
                                            (fail.outcome == "predicted")).sum()),
                 fail_lost=int((fail.category == "fail_coding_lost").sum()),
                 fail_lost_pct_of_matched=pct(int((fail.category == "fail_coding_lost").sum()), len(m)),
                 pass_but_outcompeted=int(((m.category == "pass_coding") & (m.outcome == "lost_outcompeted")).sum()),
                 rescued_but_outcompeted=int(((m.category == "fail_coding_long_rescued") &
                                              (m.outcome == "lost_outcompeted")).sum()),
                 predicted_stop=int((m.outcome == "predicted").sum()),
                 predicted_stop_pct_of_truth=pct(int((m.outcome == "predicted").sum()), len(g)),
                 fail_in_training=int(fail.in_training.sum()),
                 orf_scores_fail_pct=pct(int((m.pass_scores == 0).sum()), len(m)),
                 cds_only_scores_fail_pct=pct(int((m.pass_scores_cds == 0).sum()), len(m)))
        rows.append(r)
    return pd.DataFrame(rows)


def by_length(d):
    rows = []
    for org, g in list(d.groupby("org")) + [("ALL", d)]:
        m = g[g.matched]
        for b in LEN_ORDER:
            x = m[m.len_bin == b]
            if not len(x):
                continue
            f = x[x.fail]
            rows.append(dict(org=org, cds_len_bin=b, n=len(x), fail_coding=len(f),
                             fail_coding_pct=pct(len(f), len(x)),
                             fail_S1_le0=int((f.fail_reason == "S1<=0").sum()),
                             fail_shifted_frame=int((f.fail_reason == "shifted_frame>=S1").sum()),
                             fail_decoy_q=int((f.fail_reason == "decoy_q>0.05").sum()),
                             rescued_by_long_rule=int((f.category == "fail_coding_long_rescued").sum()),
                             rescued_and_predicted=int(((f.category == "fail_coding_long_rescued") &
                                                        (f.outcome == "predicted")).sum()),
                             still_lost=int((f.category == "fail_coding_lost").sum()),
                             still_lost_pct=pct(int((f.category == "fail_coding_lost").sum()), len(x)),
                             pass_but_outcompeted=int(((x.category == "pass_coding") &
                                                       (x.outcome == "lost_outcompeted")).sum())))
    return pd.DataFrame(rows)


def oracle(d):
    rows = []
    for org, g in list(d.groupby("org")) + [("ALL", d)]:
        m = g[g.matched]
        r = dict(org=org, n=len(m), td_fail_pct=pct(int(m.fail.sum()), len(m)),
                 td_fail_score_pct=pct(int(m.fail_score.sum()), len(m)))
        fs = m[m.fail_score]
        for k in MODELS:
            col = f"pass_{k}"
            r[f"{k}_fail_pct"] = pct(int((m[col] == 0).sum()), len(m))
            r[f"{k}_rescues_td_scorefail_pct"] = pct(int((fs[col] == 1).sum()), len(fs))
        r["td_scorefail_n"] = len(fs)
        # ORFs >= 300 nt only (no short-ORF guard)
        L = m[m.orf_len >= 300]
        r["n_ge300"] = len(L)
        r["td_fail_ge300_pct"] = pct(int(L.fail.sum()), len(L))
        for k in MODELS:
            r[f"{k}_fail_ge300_pct"] = pct(int((L[f"pass_{k}"] == 0).sum()), len(L))
        rows.append(r)
    return pd.DataFrame(rows)


def matched_controls(m, fail_mask, seed=1, k=5, nbins=20):
    """For each failing CDS, up to k passing CDSs of the same organism and CDS-length quantile bin."""
    rng = np.random.default_rng(seed)
    ctrl_idx = []
    for org, g in m.groupby("org"):
        bins = pd.qcut(np.log(g.cds_len), nbins, labels=False, duplicates="drop")
        for b in np.unique(bins):
            gb = g[bins == b]
            f = gb[fail_mask.loc[gb.index]]
            p = gb[(gb.category == "pass_coding")]
            if not len(f) or not len(p):
                continue
            n = min(len(p), k * len(f))
            ctrl_idx.extend(rng.choice(p.index.values, n, replace=False))
    return m.loc[ctrl_idx]


def features_fail_vs_pass(d):
    m = d[d.matched & (d.orf_len >= 300)]
    rows = []
    for org, g in list(m.groupby("org")) + [("ALL", m)]:
        fail = g[g.fail_score]
        if len(fail) < 3:
            continue
        ctrl = matched_controls(g, g.fail_score)
        for feat in FEATURES:
            a, b = fail[feat].dropna(), ctrl[feat].dropna()
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                p = mannwhitneyu(a, b).pvalue
            except ValueError:
                p = np.nan
            rows.append(dict(org=org, feature=feat, n_fail=len(a), n_ctrl=len(b),
                             fail_mean=a.mean(), ctrl_mean=b.mean(), fail_median=a.median(),
                             ctrl_median=b.median(), mannwhitney_p=p))
    return pd.DataFrame(rows)


def rate_by_quantile(d, feats=("gc3", "gc3_dev_train", "codon_llr_train", "seg_lowc_frac",
                               "aa_entropy", "top1_frac", "disorder_frac", "tandem3_frac",
                               "s1_per_base")):
    m = d[d.matched & (d.orf_len >= 300)]
    rows = []
    for org, g in m.groupby("org"):
        for feat in feats:
            x = g[[feat, "fail_score"]].dropna()
            try:
                q = pd.qcut(x[feat], 5, labels=False, duplicates="drop")
            except ValueError:
                continue
            for qi in sorted(q.unique()):
                s = x[q == qi]
                rows.append(dict(org=org, feature=feat, quintile=int(qi) + 1, lo=s[feat].min(),
                                 hi=s[feat].max(), n=len(s), fail_score=int(s.fail_score.sum()),
                                 fail_score_pct=pct(int(s.fail_score.sum()), len(s))))
        for feat in ("has_pfam", "has_sprot", "organelle"):
            x = g[[feat, "fail_score"]].dropna()
            for v in (0, 1):
                s = x[x[feat] == v]
                if len(s):
                    rows.append(dict(org=org, feature=feat, quintile=v, lo=v, hi=v, n=len(s),
                                     fail_score=int(s.fail_score.sum()),
                                     fail_score_pct=pct(int(s.fail_score.sum()), len(s))))
    return pd.DataFrame(rows)


def name_prefix(n):
    n = str(n)
    m = re.match(r"([A-Za-z]+)", n)
    return m.group(1).upper() if m else n


def enrichment(d):
    """Gene-name prefix (families) and description keyword enrichment among coding-test failures
    (ORFs >= 300 nt)."""
    m = d[d.matched & (d.orf_len >= 300)].copy()
    m["prefix"] = m.gene_name.map(name_prefix)
    rows = []
    for org, g in m.groupby("org"):
        tot_f, tot = int(g.fail_score.sum()), len(g)
        for pre, s in g.groupby("prefix"):
            nf = int(s.fail_score.sum())
            if nf < 3:
                continue
            odds, p = fisher_exact([[nf, len(s) - nf], [tot_f - nf, tot - len(s) - tot_f + nf]],
                                   alternative="greater")
            rows.append(dict(org=org, kind="name_prefix", term=pre, n=len(s), n_fail=nf,
                             fail_pct=pct(nf, len(s)), org_fail_pct=pct(tot_f, tot), fisher_p=p,
                             examples=",".join(s[s.fail_score].gene_name.astype(str).head(8))))
        desc = g.description.fillna("").astype(str).str.lower()
        for kw in KEYWORDS:
            hit = desc.str.contains(kw.lower(), regex=False)
            n = int(hit.sum())
            if n < 5:
                continue
            nf = int((hit & g.fail_score).sum())
            odds, p = fisher_exact([[nf, n - nf], [tot_f - nf, tot - n - tot_f + nf]], alternative="greater")
            rows.append(dict(org=org, kind="description", term=kw, n=n, n_fail=nf, fail_pct=pct(nf, n),
                             org_fail_pct=pct(tot_f, tot), fisher_p=p,
                             examples=",".join(g[hit & g.fail_score].gene_name.astype(str).head(8))))
    return pd.DataFrame(rows).sort_values(["org", "fisher_p"])


def examples(d, per_org=15):
    cols = ["org", "transcript_id", "gene_name", "description", "cds_len", "orf_len", "category",
            "outcome", "fail_reason", "S1", "S2", "S3", "contrast_rel", "gc3", "gc3_dev_train",
            "codon_llr_train", "aa_entropy", "top3_aa", "top1_frac", "seg_lowc_frac", "tandem3_frac",
            "has_pfam", "pfam_domains", "has_sprot", "sprot_hit", "pass_all", "pass_gc3match500"]
    m = d[d.fail_score]
    out = []
    for org, g in m.groupby("org"):
        out.append(g.sort_values("cds_len", ascending=False).head(per_org)[cols])
    return pd.concat(out)


def gm_summary(out):
    p = out / "gm12878_mane_assembled.tsv.gz"
    if not p.exists():
        return None, None
    g = pd.read_csv(p, sep="\t", keep_default_na=False, na_values=["NA", ""])
    rows = []
    for asm, x in g.groupby("assembly"):
        ne = x[x.eval_exact == 0]
        r = dict(assembly=asm, mane_assembled=len(x), exact=int((x.eval_exact == 1).sum()),
                 not_exact=len(ne))
        for st in ["predicted_exact", "predicted_stop", "lost_outcompeted", "fail_coding_lost", "no_orf"]:
            r[f"not_exact_{st}"] = int((ne.status == st).sum())
        f = ne[ne.status == "fail_coding_lost"]
        r["not_exact_fail_coding_pct"] = pct(len(f), len(ne))
        r["fail_S1_le0"] = int((f.fail_reason == "S1<=0").sum())
        r["fail_shifted_frame"] = int((f.fail_reason == "shifted_frame>=S1").sum())
        r["fail_short_guard"] = int((f.fail_reason == "short_orf_guard").sum())
        allfail = x[x.pass_coding == 0]
        r["all_assembled_fail_coding"] = len(allfail)
        r["all_assembled_fail_coding_pct"] = pct(len(allfail), int(x.pass_coding.notna().sum()))
        r["all_assembled_fail_but_long_rescued"] = int((allfail.long_rule == 1).sum())
        if "bench_category" in f:
            r["fail_also_fail_in_bench_human"] = int(f.bench_category.isin(
                ["fail_coding_lost", "fail_coding_long_rescued"]).sum())
            r["fail_bench_pass"] = int((f.bench_category == "pass_coding").sum())
        rows.append(r)
    ex = g[(g.eval_exact == 0) & (g.status == "fail_coding_lost")].sort_values(["assembly", "cds_len"],
                                                                               ascending=[True, False])
    return pd.DataFrame(rows), ex


def flags(d):
    """Candidate explanations as flags, for coding-test failures (ORF >= 300 nt) and their
    length-matched passing controls. Percentile thresholds are per organism over matched CDSs >= 300 nt."""
    m = d[d.matched & (d.orf_len >= 300)].copy()
    for org, g in m.groupby("org"):
        q = g[["seg_lowc_frac", "top3_frac", "tandem3_frac", "gc3"]].quantile([0.05, 0.95])
        i = g.index
        m.loc[i, "f_low_complexity"] = ((g.seg_lowc_frac > q.loc[0.95, "seg_lowc_frac"]) |
                                        (g.top3_frac > q.loc[0.95, "top3_frac"]) |
                                        (g.tandem3_frac > q.loc[0.95, "tandem3_frac"]))
        m.loc[i, "f_gc3_low"] = g.gc3 < q.loc[0.05, "gc3"]
        m.loc[i, "f_gc3_high"] = g.gc3 > q.loc[0.95, "gc3"]
    m["f_gc3_extreme"] = m.f_gc3_low.astype(bool) | m.f_gc3_high.astype(bool)
    m["f_rescued_gc3_model"] = m.pass_gc3match500 == 1
    m["f_rescued_any_truth_model"] = (m[[f"pass_{k}" for k in MODELS]] == 1).any(axis=1)
    m["f_no_sprot"] = m.has_sprot == 0
    m["f_no_pfam"] = m.has_pfam == 0
    m["f_orphan"] = (m.has_sprot == 0) & (m.has_pfam == 0)
    m["f_altframe_orf_50pct"] = m.altframe_overlap_frac >= 0.5
    # an out-of-frame LongOrfs ORF covering >= 50% of the CDS that scores higher (and > 0): a possible
    # frame error in the annotation, or overlapping coding
    m["f_altframe_orf_stronger"] = (m.altframe_overlap_frac >= 0.5) & (m.altframe_S1 > 0) & (m.altframe_S1 > m.S1)
    m["f_organelle"] = m.organelle == 1

    def attribute(r):
        if r.f_low_complexity:
            return "low_complexity_or_repetitive"
        if r.f_gc3_extreme or r.f_rescued_gc3_model:
            return "atypical_gc3_composition"
        if r.f_altframe_orf_stronger:
            return "stronger_out_of_frame_orf"
        if r.f_orphan:
            return "no_homology_orphan"
        return "unexplained"
    m["attribution"] = m.apply(attribute, axis=1)
    flag_cols = [c for c in m.columns if c.startswith("f_")]
    rows = []
    for org, g in list(m.groupby("org")) + [("ALL", m)]:
        fail = g[g.fail_score]
        if not len(fail):
            continue
        ctrl = matched_controls(g, g.fail_score)
        for c in flag_cols + ["attr:" + a for a in ["low_complexity_or_repetitive", "atypical_gc3_composition",
                                                     "stronger_out_of_frame_orf", "no_homology_orphan",
                                                     "unexplained"]]:
            if c.startswith("attr:"):
                a, b = (fail.attribution == c[5:]), (ctrl.attribution == c[5:])
            else:
                a, b = fail[c].astype(bool), ctrl[c].astype(bool)
            # attribution is defined for failures only (controls pass, so 'rescued' flags are moot)
            rows.append(dict(org=org, flag=c, n_fail=len(fail), fail_pct=pct(int(a.sum()), len(fail)),
                             n_ctrl=len(ctrl),
                             ctrl_pct=np.nan if c.startswith("attr:") else pct(int(b.sum()), len(ctrl))))
    return pd.DataFrame(rows), m


RP_RE = re.compile(r"^(RP[LS]\d+[A-Z]?\d*|Rp[ls]\d+[A-Za-z]?\d*|rp[ls]\d+)$")


def ribosomal(d):
    """Cytoplasmic ribosomal proteins (highly expressed, strong codon bias)."""
    m = d[d.matched].copy()
    desc = m.description.fillna("").astype(str).str.lower()
    m["is_rp"] = m.gene_name.astype(str).str.match(RP_RE) | (
        desc.str.contains("ribosomal protein") & ~desc.str.contains("mitochondri|chloroplast|plastid|apicoplast|kinase|methyltransferase|biogenesis|modification|acetyltransferase|hydroxylase"))
    rows = []
    for org, g in m.groupby("org"):
        rp = g[g.is_rp]
        if not len(rp):
            continue
        rows.append(dict(org=org, n_rp=len(rp), rp_fail=int(rp.fail.sum()), rp_fail_pct=pct(int(rp.fail.sum()), len(rp)),
                         rp_fail_score=int(rp.fail_score.sum()),
                         rp_fail_score_pct=pct(int(rp.fail_score.sum()), len(rp)),
                         rp_rescued_gc3_model=int((rp.fail_score & (rp.pass_gc3match500 == 1)).sum()),
                         rp_median_cds_len=rp.cds_len.median(), rp_median_gc3=rp.gc3.median(),
                         all_median_gc3=g.gc3.median(),
                         org_fail_pct_same_len=pct(int(g[(g.cds_len >= rp.cds_len.min()) & (g.cds_len <= rp.cds_len.max())].fail.sum()),
                                                   len(g[(g.cds_len >= rp.cds_len.min()) & (g.cds_len <= rp.cds_len.max())])),
                         examples_failing=",".join(rp[rp.fail].gene_name.astype(str).head(10))))
    return pd.DataFrame(rows)


def utr_effect(out, d):
    """complete (UTRs <= 250 nt) vs complete_fullutr (full UTRs): failures and outcomes by UTR length."""
    rows = []
    for org in ORGS:
        p = out / "per_cds" / f"{org}.complete_fullutr.per_cds.tsv.gz"
        if not p.exists():
            continue
        f = pd.read_csv(p, sep="\t", low_memory=False)
        f = f[f.category != "no_orf"]
        c = d[(d.org == org) & d.matched]
        j = c[["transcript_id", "category", "outcome"]].merge(
            f[["transcript_id", "category", "outcome", "utr5_len", "utr3_len", "ext5"]],
            on="transcript_id", suffixes=("_trim", "_full"))
        failf = f.category.isin(["fail_coding_long_rescued", "fail_coding_lost"])
        r = dict(org=org, n_trim=len(c), n_full=len(f),
                 fail_pct_trim=pct(int(c.fail.sum()), len(c)), fail_pct_full=pct(int(failf.sum()), len(f)),
                 predicted_pct_trim=pct(int((c.outcome == "predicted").sum()), len(c)),
                 predicted_pct_full=pct(int((f.outcome == "predicted").sum()), len(f)),
                 outcompeted_trim=int((c.outcome == "lost_outcompeted").sum()),
                 outcompeted_full=int((f.outcome == "lost_outcompeted").sum()),
                 category_changed=int((j.category_trim != j.category_full).sum()),
                 lost_only_in_full=int(((j.outcome_trim == "predicted") & (j.outcome_full != "predicted")).sum()),
                 lost_only_in_full_outcompeted=int(((j.outcome_trim == "predicted") & (j.outcome_full == "lost_outcompeted")).sum()),
                 ext5_gt0_full_pct=pct(int((f.ext5 > 0).sum()), len(f)))
        for lo, hi in [(0, 250), (250, 1000), (1000, 10 ** 9)]:
            for side in ("utr5_len", "utr3_len"):
                x = f[(f[side] >= lo) & (f[side] < hi)]
                tag = f"{side[:4]}_{lo}-{hi if hi < 10 ** 9 else 'inf'}"
                r[f"{tag}_n"] = len(x)
                r[f"{tag}_fail_pct"] = pct(int(x.category.isin(["fail_coding_long_rescued", "fail_coding_lost"]).sum()), len(x))
                r[f"{tag}_outcompeted_pct"] = pct(int((x.outcome == "lost_outcompeted").sum()), len(x))
        rows.append(r)
    return pd.DataFrame(rows)


def gm_features(out):
    p = out / "gm12878_mane_assembled.tsv.gz"
    if not p.exists():
        return None
    g = pd.read_csv(p, sep="\t")
    g = g[g.pass_coding.notna() & (g.orf_len >= 300)]
    rows = []
    for asm, x in g.groupby("assembly"):
        f, ps = x[x.pass_coding == 0], x[x.pass_coding == 1]
        for feat in ["cds_len", "gc3", "aa_entropy", "top1_frac", "seg_lowc_frac", "tandem3_frac", "disorder_frac", "S1", "contrast_rel"]:
            a, b = f[feat].dropna(), ps[feat].dropna()
            rows.append(dict(assembly=asm, feature=feat, n_fail=len(a), n_pass=len(b), fail_median=a.median(),
                             pass_median=b.median(),
                             mannwhitney_p=mannwhitneyu(a, b).pvalue if len(a) > 2 and len(b) > 2 else np.nan))
    return pd.DataFrame(rows)


def by_chromosome(d, min_n=20):
    """Coding-test failure rate (ORFs >= 300 nt) by chromosome / scaffold."""
    m = d[d.matched & (d.orf_len >= 300)]
    rows = []
    for (org, seqid), g in m.groupby(["org", "seqid"]):
        if len(g) < min_n:
            continue
        rows.append(dict(org=org, seqid=seqid, n=len(g), fail_score=int(g.fail_score.sum()),
                         fail_score_pct=pct(int(g.fail_score.sum()), len(g)), median_gc3=g.gc3.median(),
                         rescued_gc3_model=int((g.fail_score & (g.pass_gc3match500 == 1)).sum()),
                         examples=",".join(g[g.fail_score].gene_name.astype(str).head(10))))
    return pd.DataFrame(rows).sort_values(["org", "fail_score_pct"], ascending=[True, False])


def main(out):
    out = Path(out)
    d, meta = load(out)
    meta.reset_index().to_csv(out / "org_meta.tsv", sep="\t", index=False)
    tabs = dict(summary_by_org=by_org(d, meta), summary_by_length=by_length(d), oracle_models=oracle(d),
                features_fail_vs_pass=features_fail_vs_pass(d), fail_rate_by_feature=rate_by_quantile(d),
                family_keyword_enrichment=enrichment(d), examples_failing=examples(d))
    ft, fm = flags(d)
    tabs["explanation_flags"] = ft
    tabs["ribosomal_proteins"] = ribosomal(d)
    tabs["fail_by_chromosome"] = by_chromosome(d)
    u = utr_effect(out, d)
    if len(u):
        tabs["utr_effect"] = u
    gf = gm_features(out)
    if gf is not None:
        tabs["gm12878_fail_vs_pass_features"] = gf
    fm[fm.fail_score][["org", "transcript_id", "gene_name", "description", "cds_len", "fail_reason", "attribution"]
                      + [c for c in fm.columns if c.startswith("f_")]].to_csv(
        out / "failing_cds_attribution.tsv.gz", sep="\t", index=False)
    gs, gex = gm_summary(out)
    if gs is not None:
        tabs["gm12878_summary"] = gs
        tabs["gm12878_fail_examples"] = gex
    for name, t in tabs.items():
        t.to_csv(out / f"{name}.tsv", sep="\t", index=False, float_format="%.4g")
        print(f"wrote {out / name}.tsv ({len(t)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/long_orf_failures_v7"))
