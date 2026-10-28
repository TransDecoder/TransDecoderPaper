#!/usr/bin/env python3
"""E33 step 4: offline variants of the a5b24c3 homology start rules (hit depth, e-value, HSPs, Pfam threshold).

usage: python step4_03_rule_variants_v1.py <manifest.tsv> <out_dir> [dataset ...]

manifest columns (as step3_manifest.tsv): dataset, homology_dir (the -k 1 1e-5 DIAMOND + --cut_ga Pfam
results used by a5b24c3: blastp.outfmt6, pfam.domtblout), run_dir (release pytd_v7 run), truth, excluded.
Deep searches (step4_01_deep_search_v1.sh) in $P/data/start_classifier/step4/searches/<dataset>/:
  blastp.{fast,msens}.k25.e10.tsv.gz (DIAMOND -k 25 --max-hsps 0 --evalue 10, default / --more-sensitive),
  pfam.domP.domtblout.gz (hmmsearch --domE 1 with -Z 1 --domZ 1, so i-Evalue = per-domain P-value).
Rules as 08_homology_rules.py / PyTD a5b24c3, per hit (a = peptide position of the first in-frame ATG):
  R2 trim  if sstart <= 3 and |qstart - (sstart-1) - a| <= 5
  R1 keep  if qstart <= a - 5 and sstart > 1
  R3 keep  if a Pfam envelope starts at <= a - 1
  BLAST call + R3: trim with R3 = conflict -> ab initio (release) rule; trim -> trim; keep or R3 -> keep;
  nothing -> ab initio rule.
Variant dimensions:
  search  k1    : the a5b24c3 inputs (DIAMOND -k 1 --evalue 1e-5, first HSP); only with evalue 1e-5, top1
          fast / msens : the deep searches
  evalue  HSP e-value cutoff 1e-10, 1e-5, 1e-3, 1e-1, 1, 10
  hsp     best : per subject the highest-scoring HSP
          all  : per subject the HSP reaching furthest to the subject N-terminus (min sstart; ties: bits)
  depth   top1 (best subject); top-n (n = 3, 5, 10, 25 subjects ranked by best HSP bitscore) with schemes
          maj      majority of the firing hits (trim vs keep; tie = no call)
          unan     all firing hits agree
          bits     bitscore-weighted: larger summed bitscore of trim vs keep hits
          anyup    any R2 hit -> trim; else any R1 hit -> keep (keep only when no hit supports a trim)
          first    the best-ranked hit that fires decides
  pfam    ga (--cut_ga, the a5b24c3 input); dE1, dE1e-1, dE1e-3, dE1e-5: domain i-Evalue <= 1 ... 1e-5
          at a reference domZ of 10,000 (P <= 1e-4 ... 1e-9; the step-3 candidate searches had domZ ~4,400,
          the number of sequences per chunk); none (no R3)
Scoring: the release run's predictions with each candidate's keep/trim decision substituted (as
05_evaluate.py); exact = stop and start within 2 nt; F1 = 2 tp / (n_truth + n_pred) overall and per truth
type (n_pred = predictions on transcripts of the type), as 16_score_adoption.py; excluded contigs removed.
Writes <out_dir>/variants_f1.tsv.gz (dataset x variant), <out_dir>/top1_vs_k1.tsv (agreement of the deep
search's top hit at 1e-5 with the -k 1 hit).
"""
import collections
import gzip
import os
import sys

import numpy as np
import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
S4 = f"{P}/data/start_classifier/step4/searches"
FD = f"{P}/data/start_classifier/features"
TOL = 2
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
EVALS = [1e-10, 1e-5, 1e-3, 1e-1, 1.0, 10.0]
NS = [3, 5, 10, 25]
SCHEMES = ["maj", "unan", "bits", "anyup", "first"]
PFAMS = ["ga", "dE1", "dE1e-1", "dE1e-3", "dE1e-5", "none"]
DOMZ_REF = 1e4
PF_E = {"dE1": 1.0 / DOMZ_REF, "dE1e-1": 1e-1 / DOMZ_REF, "dE1e-3": 1e-3 / DOMZ_REF, "dE1e-5": 1e-5 / DOMZ_REF}

manifest, odir = sys.argv[1:3]
only = set(sys.argv[3:])
os.makedirs(odir, exist_ok=True)
man = pd.read_csv(manifest, sep="\t")


def ends(s, e, strand):
    return (s, e) if strand == "+" else (e, s)


def scorer(F, m):
    """Precompute per-candidate exactness under keep / trim; return a function keep(bool array) -> row."""
    excl = set()
    if m.excluded != "-":
        excl = {l.split("\t")[0] for l in open(m.excluded)}
    truth = {}
    with open(m.truth) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            truth[r["transcript_id"]] = (int(r["cds_start"]), int(r["cds_end"]), r["strand"], r["type"])
    preds = collections.defaultdict(list)
    with open(os.path.join(m.run_dir, "preds.tsv")) as fh:
        fh.readline()
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if x[0] not in excl:
                preds[x[0]].append((int(x[1]), int(x[2]), x[3]))
    L, o5, ol, sg = F.tx_len.values, F.orf5.values, F.orf_len.values, F.seg.values
    plus = (F.strand == "+").values
    op = np.where(plus, o5 + 1, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5)
    tr = np.where(plus, o5 + 1 + sg, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5 - sg)
    tids, strands, trim_run = F.transcript_id.values, F.strand.values, F.trim_run.values
    n = len(F)
    valid = np.zeros(n, bool)
    ex_keep, ex_trim = np.zeros(n, bool), np.zeros(n, bool)
    ctype = np.full(n, -1)
    for i in range(n):  # remove the candidates' release predictions
        tid = tids[i]
        if tid in excl:
            continue
        old = (int(tr[0][i]), int(tr[1][i])) if trim_run[i] else (int(op[0][i]), int(op[1][i]))
        preds[tid].remove((old[0], old[1], strands[i]))
        valid[i] = True

    def is_ex(t, q):
        if q[2] != t[2]:
            return False
        q5, q3 = ends(*q)
        t5, t3 = ends(*t[:3])
        return abs(q3 - t3) <= TOL and abs(q5 - t5) <= TOL

    n_truth = np.zeros(4, int)
    n_pred_t = np.zeros(4, int)
    fixed = np.zeros(4, int)
    fixed_tid = set()
    for tid, t in truth.items():
        if tid in excl:  # (truth on excluded contigs: none expected, kept as 05_evaluate does)
            pass
        k = TYPES.index(t[3])
        n_truth[k] += 1
        if any(is_ex(t, q) for q in preds.get(tid, [])):
            fixed[k] += 1
            fixed_tid.add(tid)
    seen = collections.Counter()
    for i in np.where(valid)[0]:
        tid = tids[i]
        if tid not in truth or tid in fixed_tid:
            continue
        t = truth[tid]
        ctype[i] = TYPES.index(t[3])
        ex_keep[i] = is_ex(t, (int(op[0][i]), int(op[1][i]), strands[i]))
        ex_trim[i] = is_ex(t, (int(tr[0][i]), int(tr[1][i]), strands[i]))
        if ex_keep[i] or ex_trim[i]:
            seen[tid] += 1
    assert max(seen.values(), default=0) <= 1
    # predictions per type (unchanged by decisions)
    n_pred = 0
    for tid, lst in preds.items():
        n_pred += len(lst)
    for i in np.where(valid)[0]:
        n_pred += 1
    for tid, t in truth.items():
        k = TYPES.index(t[3])
        n_pred_t[k] += len(preds.get(tid, []))
    for i in np.where(valid)[0]:
        if tids[i] in truth:
            n_pred_t[TYPES.index(truth[tids[i]][3])] += 1
    N = n_truth.sum()

    def score(keep):
        ex = np.where(keep, ex_keep, ex_trim)
        tp = fixed + np.array([(ex & (ctype == k)).sum() for k in range(4)])
        row = dict(n_truth=N, n_pred=n_pred, tp_exact=tp.sum(), f1_exact=2 * tp.sum() / (N + n_pred))
        for k, t in enumerate(TYPES):
            row[f"f1_{t}"] = 2 * tp[k] / (n_truth[k] + n_pred_t[k]) if n_truth[k] + n_pred_t[k] else np.nan
            row[f"tp_{t}"] = tp[k]
        return row

    return score, valid


def read_blast(path, cols):
    b = pd.read_csv(path, sep="\t", header=None, names=cols)
    b["order"] = np.arange(len(b))
    return b


C12 = ["q", "s", "pid", "alen", "mm", "gap", "qstart", "qend", "sstart", "send", "evalue", "bits"]
C14 = C12 + ["qlen", "slen"]


def hit_table(b, amap, E, hsp):
    """Per (query, subject): representative HSP, rank by best bitscore, outcome (+1 keep, -1 trim, 0)."""
    d = b[b.evalue <= E]
    if not len(d):
        return d.assign(rank=[], out=[])
    sb = d.groupby(["q", "s"], sort=False).bits.transform("max")
    d = d.assign(sbits=sb)
    if hsp == "best":
        rep = d.sort_values(["q", "s", "bits", "order"], ascending=[True, True, False, True])
    else:
        rep = d.sort_values(["q", "s", "sstart", "bits"], ascending=[True, True, True, False])
    rep = rep.drop_duplicates(["q", "s"], keep="first")
    first = d.drop_duplicates(["q", "s"], keep="first")[["q", "s", "order"]].rename(columns={"order": "o1"})
    rep = rep.merge(first, on=["q", "s"])
    rep = rep.sort_values(["q", "sbits", "o1"], ascending=[True, False, True])
    rep["rank"] = rep.groupby("q").cumcount()
    a = rep.q.map(amap).values
    qs, ss = rep.qstart.values, rep.sstart.values
    r2 = (ss <= 3) & (np.abs(qs - (ss - 1) - a) <= 5)
    r1 = (qs <= a - 5) & (ss > 1)
    rep["out"] = np.where(r2, -1, np.where(r1, 1, 0))
    return rep[["q", "s", "rank", "out", "sbits", "qstart", "sstart"]]


def blast_calls(rep, ids, n, scheme):
    """Call per candidate id: -1 trim, +1 keep, 0 none."""
    h = rep[rep["rank"] < n]
    if n == 1 or scheme == "top1":
        c = h[h["rank"] == 0].set_index("q").out
        return ids.map(c).fillna(0).values.astype(int)
    f = h[h.out != 0]
    cnt = f.assign(t=(f.out < 0).astype(int), k=(f.out > 0).astype(int)).groupby("q")[["t", "k"]].sum()
    nt = ids.map(cnt.t).fillna(0).values
    nk = ids.map(cnt.k).fillna(0).values
    if scheme == "maj":
        return np.where(nt > nk, -1, np.where(nk > nt, 1, 0))
    if scheme == "unan":
        return np.where((nt > 0) & (nk == 0), -1, np.where((nk > 0) & (nt == 0), 1, 0))
    if scheme == "anyup":
        return np.where(nt > 0, -1, np.where(nk > 0, 1, 0))
    if scheme == "bits":
        bt = ids.map(f[f.out < 0].groupby("q").sbits.sum()).fillna(0).values
        bk = ids.map(f[f.out > 0].groupby("q").sbits.sum()).fillna(0).values
        return np.where(bt > bk, -1, np.where(bk > bt, 1, 0))
    if scheme == "first":
        c = f.drop_duplicates("q", keep="first").set_index("q").out  # rep is sorted by rank within q
        return ids.map(c).fillna(0).values.astype(int)
    raise ValueError(scheme)


def read_pfam(path, maxE=None):
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


rows, agree = [], []
for _, m in man.iterrows():
    ds = m.dataset
    if only and ds not in only:
        continue
    sd = f"{S4}/{ds}"
    need = [f"{sd}/blastp.fast.done", f"{sd}/blastp.msens.done", f"{sd}/pfam.domP.done"]
    if not all(os.path.exists(x) for x in need):
        print(ds, "deep searches missing, skipped", file=sys.stderr)
        continue
    F = pd.read_csv(f"{FD}/{ds}.tsv.gz", sep="\t", low_memory=False)
    F["a"] = F.seg // 3 + 1
    amap = dict(zip(F.orf_id, F.a))
    ids = F.orf_id
    rel_keep = F.trim_run.values == 0
    score, valid = scorer(F, m)
    # Pfam R3 vectors
    a = F.a.values
    r3 = {}
    pga = read_pfam(f"{m.homology_dir}/pfam.domtblout")
    r3["ga"] = ids.map(pga).values
    pde = f"{sd}/pfam.domP.domtblout.gz"
    for k, e in PF_E.items():
        r3[k] = ids.map(read_pfam(pde, e)).values
    with np.errstate(invalid="ignore"):
        r3 = {k: (~np.isnan(v.astype(float))) & (v.astype(float) <= a - 1) for k, v in r3.items()}
    r3["none"] = np.zeros(len(F), bool)

    def decide(bc, pk):
        trim, keep = bc < 0, (bc > 0) | pk
        return np.where(trim & keep, rel_keep, np.where(trim, False, np.where(keep, True, rel_keep)))

    base = dict(dataset=ds)
    rows.append({**base, "variant": "release", **score(rel_keep)})
    # a5b24c3 inputs (k1)
    bk1 = read_blast(f"{m.homology_dir}/blastp.outfmt6", C12)
    bk1 = bk1[bk1.q.isin(amap)]
    bk1 = bk1.drop_duplicates("q", keep="first")
    rep_k1 = hit_table(bk1, amap, 1e-5, "best")
    bc_k1 = blast_calls(rep_k1, ids, 1, "top1")
    for pf in PFAMS:
        rows.append({**base, "variant": f"k1|1e-05|best|top1|{pf}", "search": "k1", "evalue": 1e-5, "hsp": "best",
                     "depth": "top1", "pfam": pf, "n_changed_vs_release": int((decide(bc_k1, r3[pf]) != rel_keep).sum()),
                     **score(decide(bc_k1, r3[pf]))})
    for mode in ("fast", "msens"):
        b = read_blast(f"{sd}/blastp.{mode}.k25.e10.tsv.gz", C14)
        b = b[b.q.isin(amap)]
        for E in EVALS:
            for hsp in ("best", "all"):
                rep = hit_table(b, amap, E, hsp)
                if mode == "fast" and E == 1e-5 and hsp == "best":
                    t1 = rep[rep["rank"] == 0].set_index("q")
                    k = rep_k1.set_index("q")
                    common = t1.index.intersection(k.index)
                    agree.append(dict(dataset=ds, k1_hits=len(k), deep_top1_hits=len(t1), both=len(common),
                                      same_subject=(t1.loc[common, "s"] == k.loc[common, "s"]).sum(),
                                      same_outcome=(t1.loc[common, "out"] == k.loc[common, "out"]).sum(),
                                      same_call_all=(blast_calls(rep, ids, 1, "top1") == bc_k1).sum(), n_cand=len(F)))
                depths = [("top1", 1, "top1")] + [(f"top{n}_{sc}", n, sc) for n in NS for sc in SCHEMES]
                for dname, n, sc in depths:
                    bc = blast_calls(rep, ids, n, sc)
                    for pf in PFAMS:
                        keep = decide(bc, r3[pf])
                        rows.append({**base, "variant": f"{mode}|{E:g}|{hsp}|{dname}|{pf}", "search": mode, "evalue": E,
                                     "hsp": hsp, "depth": dname, "pfam": pf,
                                     "n_changed_vs_release": int((keep != rel_keep).sum()), **score(keep)})
        print(ds, mode, "done", file=sys.stderr)
    rows.append({**base, "variant": "perfect",
                 **score(np.where(F.label == "truncated", True,
                                  np.where((F.label == "utr") & (F.true_atg_rank == 1), False, rel_keep)))})
    pd.DataFrame(rows).to_csv(f"{odir}/variants_f1.partial.tsv.gz", sep="\t", index=False, float_format="%.6g")

pd.DataFrame(rows).to_csv(f"{odir}/variants_f1.tsv.gz", sep="\t", index=False, float_format="%.6g")
pd.DataFrame(agree).to_csv(f"{odir}/top1_vs_k1.tsv", sep="\t", index=False)
print("done", len(rows), "rows", file=sys.stderr)
