#!/usr/bin/env python3
"""Compact tables for notebooks/05_coding_failures.Rmd (Fig. 5 and the coding-failure section),
from the release-rule analysis in $P/work/long_orf_failures_v7/ (05_long_orf_failures_v7.sh, then
05_long_orf_failures_v7.py summarize and 05_length_vs_homology_v7.py).

Writes data/long_orf_failures/v7/:
  - the summary tables of 05_summarize_v7.py, copied unchanged (summary_by_org, summary_by_length,
    oracle_models, explanation_flags, fail_rate_by_feature, fail_by_chromosome, ribosomal_proteins,
    features_fail_vs_pass, family_keyword_enrichment, utr_effect, org_meta), plus
    length_vs_homology_fallback, length_rule_only_predictions and alt_rules_negatives;
  - dmel_gc3_rescue.tsv: D. melanogaster failures (ORFs >= 300 nt) by GC3 quintile and for
    chromosomes 4 and Y, with the number that pass under the GC3-matched model (Fig. 5c; port of
    figures/fig5_failures/derive_dmel_gc3_rescue.py, with every coding-test failure counted, as the
    release rule has no separate short-ORF test);
  - failing_ge1kb.tsv: matched true CDSs >= 1 kb that fail the coding test (a few hundred rows; gene
    names quoted in the text);
  - reasons_by_length.tsv: failures by CDS length bin and reason, all organisms.
The per-CDS tables stay in $P/work (large).

usage: 05_fig5_tables_v7.py [work_dir]
"""
import shutil
import sys
from pathlib import Path

import pandas as pd

P = Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER")
W = Path(sys.argv[1]) if len(sys.argv) > 1 else P / "work/long_orf_failures_v7"
REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "data/long_orf_failures/v7"
OUT.mkdir(parents=True, exist_ok=True)
ORGS = ["athal", "creinhardtii", "dmel", "human", "mouse", "pfalciparum", "spombe", "tthermophila"]

COPY = ["summary_by_org", "summary_by_length", "oracle_models", "explanation_flags", "fail_rate_by_feature",
        "fail_by_chromosome", "ribosomal_proteins", "features_fail_vs_pass", "family_keyword_enrichment",
        "utr_effect", "org_meta", "length_vs_homology_fallback", "length_rule_only_predictions",
        "alt_rules_negatives"]
for name in COPY:
    shutil.copyfile(W / f"{name}.tsv", OUT / f"{name}.tsv")


def load(org):
    return pd.read_csv(W / "per_cds" / f"{org}.per_cds.tsv.gz", sep="\t", low_memory=False,
                       keep_default_na=False, na_values=["NA", ""])


FAIL = ["fail_coding_long_rescued", "fail_coding_lost"]

# Fig. 5c
d = load("dmel")
m = d[(d.category != "no_orf") & (d.orf_len >= 300)].copy()
m["fail_score"] = m.category.isin(FAIL)
m["q"] = pd.qcut(m.gc3, 5, labels=[1, 2, 3, 4, 5])
rows = []
for q, g in m.groupby("q", observed=True):
    rows.append(dict(group=f"gc3_q{q}", lo=round(g.gc3.min(), 4), hi=round(g.gc3.max(), 4), n=len(g),
                     fail_score=int(g.fail_score.sum()),
                     rescued_gc3_model=int((g.fail_score & (g.pass_gc3match500 == 1)).sum())))
for name, sel in [("chr4", m.seqid == "4"), ("chrY", m.seqid == "Y"), ("other_chr", ~m.seqid.isin(["4", "Y"]))]:
    g = m[sel]
    rows.append(dict(group=name, lo="NA", hi="NA", n=len(g), fail_score=int(g.fail_score.sum()),
                     rescued_gc3_model=int((g.fail_score & (g.pass_gc3match500 == 1)).sum())))
pd.DataFrame(rows).to_csv(OUT / "dmel_gc3_rescue.tsv", sep="\t", index=False)

# long failures and reasons by length
long_rows, reason_rows = [], []
for org in ORGS:
    d = load(org)
    mm = d[d.category != "no_orf"]
    f = mm[mm.category.isin(FAIL)]
    x = f[f.cds_len >= 1000]
    long_rows.append(x[["org", "transcript_id", "gene_name", "seqid", "cds_len", "orf_len", "L_star", "category",
                        "outcome", "fail_reason", "decoy_q", "has_sprot", "has_pfam"]])
    for (b, r), g in f.groupby(["len_bin", "fail_reason"]):
        reason_rows.append(dict(org=org, cds_len_bin=b, fail_reason=r, n=len(g)))
pd.concat(long_rows).sort_values(["org", "cds_len"], ascending=[True, False]).to_csv(
    OUT / "failing_ge1kb.tsv", sep="\t", index=False)
pd.DataFrame(reason_rows).to_csv(OUT / "reasons_by_length.tsv", sep="\t", index=False)

for p in sorted(OUT.iterdir()):
    print(f"{p.name}\t{p.stat().st_size}", file=sys.stderr)
