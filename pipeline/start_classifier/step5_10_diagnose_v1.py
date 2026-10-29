#!/usr/bin/env python3
"""E33 step 5: why do UniRef50 top hits change start decisions? Candidates whose a5b24c3 decision differs between
Swiss-Prot and UniRef50 -k 1 hits (start-decision isolate), with the UniRef50 top hit's properties.

usage: python step5_10_diagnose_v1.py <step5_manifest.tsv> <names.dmp> <out_dir>

For each dataset and each labelled candidate (truncated: keep is right; UTR extension with the true start at
the first ATG: trim is right), the a5b24c3 call (R1 keep / R2 trim / none) from each database's top hit.
Writes
  diag_flips.tsv     per dataset: candidates right with Swiss-Prot and wrong with UniRef50 and vice versa,
                     split by the UniRef50 top-hit call (R1 / R2 / none) and by whether R3 fired
  diag_r1_wrong.tsv  UniRef50 R1 keeps on candidates whose right call is trim (UTR extensions): top-hit
                     properties (cluster size, subject begins with M, identity, how far upstream the subject's
                     N-terminus maps, representative species) summarised per dataset
  diag_r1_species.tsv representative species of those wrong R1 hits (top 15 per organism)
"""
import collections
import sys

import numpy as np
import pandas as pd

sys.argv[0] = sys.argv[0]
man, names_dmp, odir = sys.argv[1:4]
P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
B = f"{P}/data/start_classifier/step5"
C12 = ["q", "s", "pid", "alen", "mm", "gap", "qstart", "qend", "sstart", "send", "evalue", "bits"]
UCOLS = C12 + ["qlen", "slen", "staxids", "s_first", "s_n", "s_lca", "s_rep_taxid", "s_own"]
M = pd.read_csv(man, sep="\t")
names = {}
for line in open(names_dmp):
    f = line.split("\t|\t")
    if f[3].startswith("scientific name"):
        names[f[0]] = f[1]


def call(b, amap):
    a = b.q.map(amap).values
    r2 = (b.sstart.values <= 3) & (np.abs(b.qstart.values - (b.sstart.values - 1) - a) <= 5)
    r1 = (b.qstart.values <= a - 5) & (b.sstart.values > 1)
    return pd.Series(np.where(r2, "R2", np.where(r1, "R1", "none")), index=b.q.values)


flips, r1rows, species = [], [], collections.defaultdict(collections.Counter)
for _, m in M.iterrows():
    F = pd.read_csv(m.features, sep="\t", low_memory=False)
    F["a"] = F.seg // 3 + 1
    amap = dict(zip(F.orf_id, F.a))
    lab = np.where(F.label == "truncated", "keep", np.where((F.label == "utr") & (F.true_atg_rank == 1), "trim", "-"))
    F["right"] = lab
    sp = pd.read_csv(f"{m.homology_dir}/blastp.outfmt6", sep="\t", header=None, names=C12)
    sp = sp[sp.q.isin(amap)].drop_duplicates("q")
    sp = sp[sp.evalue <= 1e-5]
    u = pd.read_csv(f"{B}/hits/uniref50/{m.dataset}/k1.tsv.gz", sep="\t", header=None, names=UCOLS,
                    dtype={"staxids": str, "s_lca": str, "s_rep_taxid": str, "s_own": str})
    u = u[u.q.isin(amap)].drop_duplicates("q")
    u = u[u.evalue <= 1e-5]
    F["sp"] = F.orf_id.map(call(sp, amap)).fillna("nohit")
    F["ur"] = F.orf_id.map(call(u, amap)).fillna("nohit")
    good = {"keep": "R1", "trim": "R2"}
    bad = {"keep": "R2", "trim": "R1"}
    L = F[F.right != "-"]
    for side in ("sprot_right_uniref_wrong", "uniref_right_sprot_wrong"):
        a_, b_ = ("sp", "ur") if side.startswith("sprot") else ("ur", "sp")
        x = L[(L[a_] == L.right.map(good)) & (L[b_] == L.right.map(bad))]
        y = L[(L[a_] == L.right.map(good)) & (L[b_].isin(["none", "nohit"]))]
        flips.append(dict(dataset=m.dataset, comparison=side, right_vs_opposite_call=len(x),
                          right_vs_no_call=len(y)))
    w = L[(L.right == "trim") & (L.ur == "R1")].merge(u, left_on="orf_id", right_on="q")
    if len(w):
        up = w.a - (w.qstart - (w.sstart - 1))
        r1rows.append(dict(dataset=m.dataset, n_wrong_R1=len(w), frac_sprot_also_R1=(w.sp == "R1").mean(),
                           frac_sprot_R2=(w.sp == "R2").mean(), frac_sprot_nohit=(w.sp == "nohit").mean(),
                           median_pid=w.pid.median(), frac_pid_ge95=(w.pid >= 95).mean(),
                           frac_subject_M=(w.s_first == "M").mean(), median_cluster_n=w.s_n.median(),
                           frac_singleton=(w.s_n == 1).mean(), median_aa_upstream=up.median(),
                           frac_upstream_ge20=(up >= 20).mean(),
                           frac_cluster_has_own_member=w.s_own.astype(str).str.contains(
                               {"human": "9606", "mouse": "10090", "athal": "3702", "dmel": "7227",
                                "spombe": "284812", "creinhardtii": "3055", "pfalciparum": "36329",
                                "tthermophila": "312017"}[m.org]).mean()))
        for t in w.s_rep_taxid:
            species[m.org][names.get(str(t), str(t))] += 1
pd.DataFrame(flips).to_csv(f"{odir}/diag_flips.tsv", sep="\t", index=False)
pd.DataFrame(r1rows).to_csv(f"{odir}/diag_r1_wrong.tsv", sep="\t", index=False, float_format="%.3f")
rows = []
for org, c in species.items():
    tot = sum(c.values())
    for sp_, n in c.most_common(15):
        rows.append(dict(org=org, rep_species=sp_, n=n, frac=n / tot))
pd.DataFrame(rows).to_csv(f"{odir}/diag_r1_species.tsv", sep="\t", index=False, float_format="%.3f")
