#!/usr/bin/env python3
"""E33 step 5: UniRef50 vs Swiss-Prot hits for the homology start assist (start-decision isolate).

usage: python step5_06_rule_variants_v1.py <step5_manifest.tsv> <out_dir> [dataset ...]

Manifest columns: dataset, org, homology_dir (Swiss-Prot -k 1 1e-5 blastp.outfmt6 + --cut_ga pfam.domtblout,
the a5b24c3 inputs), sprot_k25 (DIAMOND fast -k 25 --max-hsps 0 --evalue 10 vs sprot_minus_<taxid>),
run_dir (release pytd_v7 run), truth, excluded, features (01_features.py table).
UniRef50 hits (step5_05_expand_hits_v1.py): $B/hits/uniref50/<dataset>/{k1,k25}.tsv.gz, DIAMOND vs UniRef50
2026_03 with --taxon-exclude <organism>, -k 1 --evalue 1e-5 / -k 25 --max-hsps 0 --evalue 1e-5; columns
C14 + staxids, s_first (subject residue 1), s_n (cluster size), s_lca, s_rep_taxid, s_own (benchmark
organisms among the cluster members).
Same Pfam (--cut_ga) hits for both databases.

Rules (a = peptide position of the first in-frame ATG; per hit = best HSP of the subject, e <= 1e-5):
  a5b24c3   R2 trim if sstart <= 3 and |qstart - (sstart-1) - a| <= 5; R1 keep if qstart <= a-5 and sstart > 1;
            R3 keep if a Pfam envelope starts at <= a-1. BLAST call + R3: trim with R3 = conflict -> ab initio;
            trim -> trim; keep or R3 -> keep; nothing -> ab initio. (step4_03_rule_variants_v1.py)
  top-n     the top n = 3, 5, 10, 25 subjects (by best-HSP bitscore), R1/R2 per hit, combined by maj, unan,
            bits, anyup, first (as step 4)
  cons      consensus N-terminus: each hit projects the subject's residue 1 onto the query,
            p = qstart - (sstart - 1); class trim if |p - a| <= 5, keep if p < a - 5 (i.e. >= 6 aa upstream;
            p = a - 5 counts as trim), else neutral. Bitscore-weighted: trim if the trim hits carry > 1/2 of
            the weight, keep if the keep hits do; weight denominator = all hits considered (`all`) or the
            trim + keep hits only (`fire`). n = 1, 3, 5, 10, 25; `cap30` uses only hits with sstart <= 30.
  M filter  none; `r2`: a trim class (R2 or cons trim) needs a subject beginning with M; `all`: hits whose
            subject does not begin with M are dropped before ranking.
Sources: sprot|k1 (the a5b24c3 inputs), sprot|k25 (step-4 fast search, filtered to e <= 1e-5),
uniref50|k1, uniref50|k25.
Scoring as step4_03 (decision substitution in the release run's predictions; exact F1 overall and per
type with negatives). Decision metrics on labelled candidates: wrong trims (trimmed, truth truncated) and
wrong keeps (kept open, truth UTR extension with the true start at the first ATG).
Writes <out_dir>/variants_f1.tsv.gz and <out_dir>/hit_stats.tsv (coverage, identity, cluster stats).
"""
import collections
import gzip
import os
import sys

import numpy as np
import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
B = f"{P}/data/start_classifier/step5"
HITS = os.environ.get("STEP5_HITS", f"{B}/hits/uniref50")  # override only for tests
TOL = 2
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
NS = [3, 5, 10, 25]
SCHEMES = ["maj", "unan", "bits", "anyup", "first"]
CONS_NS = [1, 3, 5, 10, 25]
MFILTS = ["none", "r2", "all"]
ORG_TAXID = {"human": "9606", "mouse": "10090", "athal": "3702", "dmel": "7227", "spombe": "284812",
             "creinhardtii": "3055", "pfalciparum": "36329", "tthermophila": "312017"}

manifest, odir = sys.argv[1:3]
only = set(sys.argv[3:])
os.makedirs(odir, exist_ok=True)
man = pd.read_csv(manifest, sep="\t")


def ends(s, e, strand):
    return (s, e) if strand == "+" else (e, s)


def scorer(F, m):
    """As step4_03_rule_variants_v1.py: per-candidate exactness under keep / trim; keep(bool array) -> row."""
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
    for i in range(n):
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
    n_pred = sum(len(v) for v in preds.values()) + int(valid.sum())
    for tid, t in truth.items():
        n_pred_t[TYPES.index(t[3])] += len(preds.get(tid, []))
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


C12 = ["q", "s", "pid", "alen", "mm", "gap", "qstart", "qend", "sstart", "send", "evalue", "bits"]
C14 = C12 + ["qlen", "slen"]
UCOLS = C14 + ["staxids", "s_first", "s_n", "s_lca", "s_rep_taxid", "s_own"]

SP_FIRST = {}
with gzip.open(f"{B}/db/sprot_seqinfo.tsv.gz", "rt") as fh:
    for line in fh:
        x = line.rstrip("\n").split("\t")
        SP_FIRST[x[0]] = x[1]


def read_pfam(path):
    pf = {}
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            x = line.split()
            pf[x[0]] = min(pf.get(x[0], 10 ** 9), int(x[19]))
    return pf


def hit_table(b, amap):
    """Per (query, subject): the best-bitscore HSP at e <= 1e-5, ranked by bitscore within the query (file order
    breaks ties, as DIAMOND's own order)."""
    d = b[b.evalue <= 1e-5].copy()
    d["order"] = np.arange(len(d))
    rep = d.sort_values(["q", "s", "bits", "order"], ascending=[True, True, False, True])
    rep = rep.drop_duplicates(["q", "s"], keep="first")
    first = d.drop_duplicates(["q", "s"], keep="first")[["q", "s", "order"]].rename(columns={"order": "o1"})
    rep = rep.drop(columns="order").merge(first, on=["q", "s"])
    rep = rep.sort_values(["q", "bits", "o1"], ascending=[True, False, True]).reset_index(drop=True)
    rep["a"] = rep.q.map(amap).values
    return rep


def classify(rep, mfilt):
    """Add rank, R1/R2 outcome `out` (+1 keep, -1 trim, 0) and consensus class `cls` under an M filter."""
    r = rep
    if mfilt == "all":
        r = r[r.s_first == "M"]
    r = r.copy()
    r["rank"] = r.groupby("q").cumcount()
    a, qs, ss = r.a.values, r.qstart.values, r.sstart.values
    m_ok = (r.s_first.values == "M") if mfilt == "r2" else np.ones(len(r), bool)
    r2 = (ss <= 3) & (np.abs(qs - (ss - 1) - a) <= 5) & m_ok
    r1 = (qs <= a - 5) & (ss > 1)
    r["out"] = np.where(r2, -1, np.where(r1, 1, 0))
    p = qs - (ss - 1)
    ct = (np.abs(p - a) <= 5) & m_ok
    ck = p < a - 5
    r["cls"] = np.where(ct, -1, np.where(ck, 1, 0))
    return r


def blast_calls(rep, ids, n, scheme):
    h = rep[rep["rank"] < n]
    if n == 1:
        c = h.set_index("q").out
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
        bt = ids.map(f[f.out < 0].groupby("q").bits.sum()).fillna(0).values
        bk = ids.map(f[f.out > 0].groupby("q").bits.sum()).fillna(0).values
        return np.where(bt > bk, -1, np.where(bk > bt, 1, 0))
    if scheme == "first":
        c = f.drop_duplicates("q", keep="first").set_index("q").out
        return ids.map(c).fillna(0).values.astype(int)
    raise ValueError(scheme)


def cons_calls(rep, ids, n, denom, cap):
    h = rep[rep["rank"] < n]
    if cap:
        h = h[h.sstart <= cap]
    if not len(h):
        return np.zeros(len(ids), int)
    w = h.bits.values
    g = pd.DataFrame(dict(q=h.q.values, wt=np.where(h.cls.values < 0, w, 0.0), wk=np.where(h.cls.values > 0, w, 0.0),
                          w=w)).groupby("q").sum()
    wt = ids.map(g.wt).fillna(0).values
    wk = ids.map(g.wk).fillna(0).values
    W = ids.map(g.w).fillna(0).values if denom == "all" else wt + wk
    return np.where((W > 0) & (wt > W / 2), -1, np.where((W > 0) & (wk > W / 2), 1, 0))


def load_source(m, db, search, amap):
    ds = m.dataset
    if db == "sprot":
        if search == "k1":
            b = pd.read_csv(f"{m.homology_dir}/blastp.outfmt6", sep="\t", header=None, names=C12)
            b = b[b.q.isin(amap)].drop_duplicates("q", keep="first")
        else:
            b = pd.read_csv(m.sprot_k25, sep="\t", header=None, names=C14)
            b = b[b.q.isin(amap)]
        b["s_first"] = b.s.map(SP_FIRST).fillna("?")
        b["s_own"] = "-"
        b["s_n"] = np.nan
    else:
        b = pd.read_csv(f"{HITS}/{ds}/{search}.tsv.gz", sep="\t", header=None, names=UCOLS,
                        dtype={"staxids": str, "s_lca": str, "s_rep_taxid": str, "s_own": str})
        b = b[b.q.isin(amap)]
        if search == "k1":
            b = b.drop_duplicates("q", keep="first")
    return b


rows, stats = [], []
for _, m in man.iterrows():
    ds = m.dataset
    if only and ds not in only:
        continue
    F = pd.read_csv(m.features, sep="\t", low_memory=False)
    F["a"] = F.seg // 3 + 1
    amap = dict(zip(F.orf_id, F.a))
    ids = F.orf_id
    rel_keep = F.trim_run.values == 0
    trunc = (F.label == "truncated").values
    utr1 = ((F.label == "utr") & (F.true_atg_rank == 1)).values
    score, valid = scorer(F, m)
    pga = read_pfam(f"{m.homology_dir}/pfam.domtblout")
    pv = ids.map(pga).values.astype(float)
    with np.errstate(invalid="ignore"):
        r3 = (~np.isnan(pv)) & (pv <= F.a.values - 1)

    def decide(bc):
        trim, keep = bc < 0, (bc > 0) | r3
        return np.where(trim & keep, rel_keep, np.where(trim, False, np.where(keep, True, rel_keep)))

    def emit(variant, bc, extra):
        keep = decide(bc)
        rows.append(dict(dataset=ds, variant=variant, **extra,
                         n_changed_vs_release=int((keep != rel_keep).sum()),
                         n_blast_trim=int((bc < 0).sum()), n_blast_keep=int((bc > 0).sum()), n_r3=int(r3.sum()),
                         n_conflict=int(((bc < 0) & r3).sum()),
                         wrong_trim=int((~keep & trunc).sum()), wrong_keep=int((keep & utr1).sum()),
                         right_trim=int((~keep & utr1).sum()), right_keep=int((keep & trunc).sum()),
                         **score(keep)))

    rows.append(dict(dataset=ds, variant="release", n_changed_vs_release=0,
                     wrong_trim=int((~rel_keep & trunc).sum()), wrong_keep=int((rel_keep & utr1).sum()),
                     right_trim=int((~rel_keep & utr1).sum()), right_keep=int((rel_keep & trunc).sum()),
                     **score(rel_keep)))
    perf = np.where(trunc, True, np.where(utr1, False, rel_keep))
    rows.append(dict(dataset=ds, variant="perfect", n_changed_vs_release=int((perf != rel_keep).sum()), **score(perf)))
    nlab_t, nlab_u = int(trunc.sum()), int(utr1.sum())
    for db in ("sprot", "uniref50"):
        for search in ("k1", "k25"):
            b = load_source(m, db, search, amap)
            rep = hit_table(b, amap)
            top = rep.drop_duplicates("q", keep="first")
            hit = ids.isin(top.q).values
            st = dict(dataset=ds, db=db, search=search, n_cand=len(F), n_with_hit=int(hit.sum()),
                      frac_with_hit=hit.mean(), n_trunc=nlab_t, frac_trunc_with_hit=hit[trunc].mean() if nlab_t else np.nan,
                      n_utr1=nlab_u, frac_utr1_with_hit=hit[utr1].mean() if nlab_u else np.nan,
                      top_pid_ge95=(top.pid >= 95).mean(), top_pid_ge99=(top.pid >= 99).mean(),
                      top_pid_median=top.pid.median(), top_sstart_le3=(top.sstart <= 3).mean(),
                      top_subject_M=(top.s_first == "M").mean(),
                      hits_per_query=rep.groupby("q").size().mean() if len(rep) else 0)
            if db == "uniref50":
                own = top.s_own.astype(str).str.split(",").apply(lambda z: ORG_TAXID[m.org] in z)
                st.update(top_cluster_has_own_member=own.mean(), top_cluster_singleton=(top.s_n == 1).mean(),
                          top_cluster_n_median=top.s_n.median(),
                          top_own_member_pid_ge95=(own & (top.pid >= 95)).mean())
            stats.append(st)
            for mf in MFILTS:
                r = classify(rep, mf)
                base = dict(db=db, search=search, mfilt=mf)
                bc1 = blast_calls(r, ids, 1, "top1")
                t1 = r[r["rank"] == 0].set_index("q")
                emit(f"{db}|{search}|top1|{mf}", bc1,
                     dict(base, depth="top1", n_R1=int((ids.map(t1.out) == 1).sum()), n_R2=int((ids.map(t1.out) == -1).sum())))
                if search == "k1":
                    continue
                for n in NS:
                    for sc in SCHEMES:
                        emit(f"{db}|{search}|top{n}_{sc}|{mf}", blast_calls(r, ids, n, sc), dict(base, depth=f"top{n}_{sc}"))
                for n in CONS_NS:
                    for den in ("all", "fire"):
                        for cap in (None, 30):
                            dn = f"cons{n}_{den}" + ("_cap30" if cap else "")
                            emit(f"{db}|{search}|{dn}|{mf}", cons_calls(r, ids, n, den, cap), dict(base, depth=dn))
            print(ds, db, search, "done", file=sys.stderr, flush=True)
    pd.DataFrame(rows).to_csv(f"{odir}/variants_f1.partial.tsv.gz", sep="\t", index=False, float_format="%.6g")

pd.DataFrame(rows).to_csv(f"{odir}/variants_f1.tsv.gz", sep="\t", index=False, float_format="%.6g")
pd.DataFrame(stats).to_csv(f"{odir}/hit_stats.tsv", sep="\t", index=False, float_format="%.4f")
print("done", len(rows), "rows", file=sys.stderr)
