#!/usr/bin/env python3
"""E33 step 4: does deeper / looser search find evidence for the candidates a5b24c3 decides wrongly?

usage: python step4_05_evidence_v1.py <manifest.tsv> <out_tsv> [dataset ...]

Labelled candidates (01_features.py): truncated (keep open is right) and UTR extensions whose true start is
the first ATG (trim is right). Decision of a5b24c3 = rules on the -k 1 / --cut_ga inputs (as
step4_03_rule_variants_v1.py 'k1|1e-05|best|top1|ga'). For each candidate and evidence source, whether any
hit supports the right call and whether any hit supports the wrong call:
  k1_top1      the a5b24c3 BLAST hit (R1 / R2) and Pfam GA (R3)
  <mode>_<E>   DIAMOND <mode> (fast, msens) top 25 subjects, every HSP with e-value <= E (1e-5, 10) is
               tested with R1 / R2; Pfam domains with P <= 1e-4 (i-Evalue 1 at domZ 10,000) tested with R3
Rows: dataset x source x (a5b24c3 right / wrong) with counts of candidates that have support for the right
call only, the wrong call only, both, neither.
"""
import gzip
import sys

import numpy as np
import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
S4 = f"{P}/data/start_classifier/step4/searches"
FD = f"{P}/data/start_classifier/features"
C12 = ["q", "s", "pid", "alen", "mm", "gap", "qstart", "qend", "sstart", "send", "evalue", "bits"]
C14 = C12 + ["qlen", "slen"]
manifest, out = sys.argv[1:3]
only = set(sys.argv[3:])
man = pd.read_csv(manifest, sep="\t")


def pfam_min_env(path, maxE=None):
    pf = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            x = line.split()
            if maxE is not None and float(x[12]) > maxE:
                continue
            pf[x[0]] = min(pf.get(x[0], 10 ** 9), int(x[19]))
    return pf


def support(b, amap, ids):
    a = b.q.map(amap).values
    qs, ss = b.qstart.values, b.sstart.values
    r2 = (ss <= 3) & (np.abs(qs - (ss - 1) - a) <= 5)
    r1 = (qs <= a - 5) & (ss > 1)
    t = ids.isin(set(b.q[r2])).values
    k = ids.isin(set(b.q[r1])).values
    return t, k


rows = []
for _, m in man.iterrows():
    ds = m.dataset
    if only and ds not in only:
        continue
    F = pd.read_csv(f"{FD}/{ds}.tsv.gz", sep="\t", low_memory=False)
    F["a"] = F.seg // 3 + 1
    amap = dict(zip(F.orf_id, F.a))
    ids = F.orf_id
    a = F.a.values
    trunc = (F.label == "truncated").values
    utr1 = ((F.label == "utr") & (F.true_atg_rank == 1)).values
    lab = trunc | utr1
    rel_keep = F.trim_run.values == 0
    # a5b24c3
    bk = pd.read_csv(f"{m.homology_dir}/blastp.outfmt6", sep="\t", header=None, names=C12)
    bk = bk[bk.q.isin(amap)].drop_duplicates("q", keep="first")
    t1, k1 = support(bk, amap, ids)
    with np.errstate(invalid="ignore"):
        pe = ids.map(pfam_min_env(f"{m.homology_dir}/pfam.domtblout")).values.astype(float)
        r3 = ~np.isnan(pe) & (pe <= a - 1)
        pe2 = ids.map(pfam_min_env(f"{S4}/{ds}/pfam.domP.domtblout.gz", 1e-4)).values.astype(float)
        r3d = ~np.isnan(pe2) & (pe2 <= a - 1)
    keep = k1 | r3
    dec = np.where(t1 & keep, rel_keep, np.where(t1, False, np.where(keep, True, rel_keep)))
    right = np.where(trunc, dec, ~dec)
    srcs = {"k1_top1": (t1, k1 | r3)}
    for mode in ("fast", "msens"):
        b = pd.read_csv(f"{S4}/{ds}/blastp.{mode}.k25.e10.tsv.gz", sep="\t", header=None, names=C14)
        b = b[b.q.isin(amap)]
        for E in (1e-5, 10.0):
            t, k = support(b[b.evalue <= E], amap, ids)
            srcs[f"{mode}_{E:g}"] = (t, k | r3d)
    for name, (t, k) in srcs.items():
        sup_right = np.where(trunc, k, t)
        sup_wrong = np.where(trunc, t, k)
        for st, msk in (("a5b24c3_right", lab & right), ("a5b24c3_wrong", lab & ~right)):
            rows.append(dict(dataset=ds, source=name, a5b24c3=st, n=int(msk.sum()),
                             right_only=int((msk & sup_right & ~sup_wrong).sum()),
                             wrong_only=int((msk & ~sup_right & sup_wrong).sum()),
                             both=int((msk & sup_right & sup_wrong).sum()),
                             neither=int((msk & ~sup_right & ~sup_wrong).sum())))
    print(ds, "done", file=sys.stderr)
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
