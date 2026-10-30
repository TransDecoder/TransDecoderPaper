#!/usr/bin/env python3
"""E35: one-sided Kozak term variants of the start rule, scored offline on finished PyTD runs.

usage: python 02_evaluate_v1.py <manifest.tsv> <candidates_dir> <out_dir> [dataset ...]

Rule (release, V0): trim the 5'-open ORF to its first in-frame ATG when
    UTR(u) = 0.4|u| + 17 K  >  coding(u) = C(u)/mu.
Variants change only the Kozak term (all other quantities as the run computed them, 01_candidates_v1.py):
    V0      release: 17 K
    V1      17 max(K, 0)
    V2_05   17 max(K, -0.5)
    V2_1    17 max(K, -1)
    V3_8.5  17 K if K >= 0 else 8.5 K
    V3_4    17 K if K >= 0 else 4 K
    V4      17 max(K, 0) if coding(u) <= 0 else 17 K
Every variant can only raise the UTR score, so it can only turn keep-open into trim.
Homology arm: the run's homology call (trim / keep from R1-R3, via PyTD's homology_start_call) overrides
the ab initio rule where it fires; conflicts and no call fall back to the (variant) ab initio rule.

Scoring (as score_predictions.py / 16_score_adoption.py / 13b_depth_collect.py): the run's preds.tsv with
each candidate's decision substituted (stops never change); predictions on excluded contigs removed; a
prediction matches a truth CDS when its stop is within 2 nt (same strand), exact when its start is too;
one credit per truth CDS. Overall F1 = 2 tp / (n_truth + all scored predictions). Per type:
f1_<type> = 2 tp_exact / (n_truth of type + predictions on transcripts of that type); f1w_<type> ("with
negatives", 13b_depth_collect.py) also counts in the denominator the predictions on negative contigs whose
own ORF is of that type (ATG start / stop end). Benchmark inputs have no negatives (f1w = f1).
Check: V0 decisions equal the run's own decisions for every candidate, and V0 predictions equal preds.tsv.

Writes <out_dir>/scores.tsv (dataset x variant), changes.tsv (dataset x variant x truth type x class of the
changed decisions), changed_candidates.tsv.gz (every changed decision with K, |u|, coding, class), checks.tsv.
Change classes (truth CDS on the transcript, judged on the transcript's best match before/after):
    w2r   wrong -> right (not exact before, exact after)
    r2w   right -> wrong
    w2w   stop matched, exact neither before nor after; split by where the true start is relative to the
          first ATG: 'later' (true start downstream: a later ATG), 'upstream' (true start inside u or beyond
          the open end, e.g. a truncated CDS)
    nostop_truth  the candidate's stop does not match the transcript's truth CDS
    no_truth      no truth CDS on the transcript (negatives, unlabelled)
"""
import collections
import gzip
import os
import sys

import numpy as np
import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
sys.path.insert(0, f"{P}/TransDecoderPaper/benchmark_v2/scripts")
from score_predictions import TOL, ends, read_preds, read_truth  # noqa: E402

TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTN", "TGCAN")
VARIANTS = ["V0", "V1", "V2_05", "V2_1", "V3_8.5", "V3_4", "V4"]


def kozak_term(v, K, coding):
    if v == "V0":
        return 17 * K
    if v == "V1":
        return 17 * np.maximum(K, 0.0)
    if v == "V2_05":
        return 17 * np.maximum(K, -0.5)
    if v == "V2_1":
        return 17 * np.maximum(K, -1.0)
    if v == "V3_8.5":
        return np.where(K >= 0, 17 * K, 8.5 * K)
    if v == "V3_4":
        return np.where(K >= 0, 17 * K, 4 * K)
    if v == "V4":
        return np.where(coding <= 0, 17 * np.maximum(K, 0.0), 17 * K)
    raise ValueError(v)


def decide(v, C):
    seg, coding, K = C.seg.values.astype(float), C.coding.values, C.kozak.values
    trim_ab = coding < 0.4 * seg + kozak_term(v, K, coding)
    h = C.hcall.values
    return np.where(h == "trim", True, np.where(h == "keep", False, trim_ab))


def read_fasta(path, keep):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in keep:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in keep:
                buf.append(line.strip())
    if name in keep:
        seqs[name] = "".join(buf).upper()
    return seqs


def orf_type(seq, s, e, strand):
    c = seq[s - 1:e]
    if strand == "-":
        c = c.translate(COMP)[::-1]
    return {(True, True): "complete", (False, True): "5prime_partial",
            (True, False): "3prime_partial", (False, False): "internal"}[(c[:3] == "ATG", c[-3:] in STOPS)]


def best_level(t, plist):
    best = None
    for p in plist:
        if p[2] != t[2]:
            continue
        p5, p3 = ends(*p)
        t5, t3 = ends(*t[:3])
        if abs(p3 - t3) <= TOL:
            lvl = "exact" if abs(p5 - t5) <= TOL else "stop"
            if best != "exact":
                best = lvl
    return best


def score(truth, preds, negs, negseq):
    n_pred = sum(len(v) for v in preds.values())
    row = dict(n_truth=len(truth), n_pred=n_pred, tp_stop=0, tp_exact=0)
    per = {t: [0, 0, 0] for t in TYPES}
    for tid, t in truth.items():
        b = best_level(t, preds.get(tid, []))
        row["tp_stop"] += b is not None
        row["tp_exact"] += b == "exact"
        a = per[t[3]]
        a[0] += 1
        a[1] += len(preds.get(tid, []))
        a[2] += b == "exact"
    negt = collections.Counter()
    for tid in negs:
        for q in preds.get(tid, []):
            negt[orf_type(negseq[tid], *q)] += 1
    row["n_pred_neg"] = sum(negt.values())
    row["f1_stop"] = 2 * row["tp_stop"] / (row["n_truth"] + n_pred)
    row["f1_exact"] = 2 * row["tp_exact"] / (row["n_truth"] + n_pred)
    for t, (n, npr, tp) in per.items():
        row[f"n_{t}"], row[f"tp_{t}"], row[f"npred_{t}"], row[f"negpred_{t}"] = n, tp, npr, negt[t]
        row[f"f1_{t}"] = 2 * tp / (n + npr) if n + npr else np.nan
        # 13b_depth_collect.py: F1 from sensitivity tp/n and precision tp/(npr + negt)
        sens = tp / n if n else 0.0
        den = npr + negt[t]
        prec = tp / den if den else 0.0
        row[f"f1w_{t}"] = 2 * sens * prec / (sens + prec) if sens + prec else (np.nan if not n else 0.0)
    return row


def main():
    manifest, cdir, odir = sys.argv[1:4]
    only = set(sys.argv[4:])
    os.makedirs(odir, exist_ok=True)
    man = pd.read_csv(manifest, sep="\t")
    srows, crows, checks, chg = [], [], [], []
    for _, m in man.iterrows():
        ds = m.dataset
        if only and ds not in only:
            continue
        C = pd.read_csv(f"{cdir}/{ds}.tsv.gz", sep="\t", keep_default_na=False,
                        dtype={"coding": float, "kozak": float})
        excl = set() if m.excluded == "-" else {l.split("\t")[0].rstrip("\n") for l in open(m.excluded)}
        negs = set() if m.negatives == "-" else {l.strip() for l in open(m.negatives) if l.strip()}
        negs -= excl
        truth = read_truth(m.truth)
        base = {k: v for k, v in read_preds(os.path.join(m.run_dir, "preds.tsv")).items() if k not in excl}
        negseq = read_fasta(os.path.join(m.run_dir, "input.fa"), negs)
        negs = {t for t in negs if t in negseq}
        L, o5, ol, sg = C.tx_len.values, C.orf5.values, C.orf_len.values, C.seg.values
        plus = (C.strand == "+").values
        op = np.where(plus, o5 + 1, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5)
        tr = np.where(plus, o5 + 1 + sg, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5 - sg)
        tids, strands, trim_run = C.transcript_id.values, C.strand.values, C.trim_run.values.astype(bool)
        inc = np.array([t not in excl for t in tids])
        # true start relative to the first ATG (for changed-decision classes)
        d0 = decide("V0", C)
        n_mis = int((d0 != trim_run).sum())
        n_mis_rel = int((d0 != C.trim_rel.values.astype(bool)).sum())
        # Kozak-veto candidates under V0 (kept open by the ab initio rule only because K < 0)
        seg_f, cod, K = sg.astype(float), C.coding.values, C.kozak.values
        ab0 = cod < 0.4 * seg_f + 17 * K
        veto = (~ab0) & (K < 0) & (cod < 0.4 * seg_f)
        checks.append(dict(dataset=ds, candidates=len(C), on_excluded=int((~inc).sum()), v0_vs_run_mismatch=n_mis,
                           v0_vs_extract_mismatch=n_mis_rel, kept_open_run=int((~trim_run & inc).sum()),
                           kept_ab_initio=int((~ab0 & inc).sum()), kozak_veto=int((veto & inc).sum()),
                           hcall_trim=int((C.hcall == "trim").sum()), hcall_keep=int((C.hcall == "keep").sum()),
                           hcall_conflict=int((C.hcall == "conflict").sum())))
        for v in VARIANTS:
            dv = decide(v, C)
            preds = {k: list(x) for k, x in base.items()}
            changed = np.where((dv != trim_run) & inc)[0]
            for i in changed:
                tid = tids[i]
                old = (int(tr[0][i]), int(tr[1][i])) if trim_run[i] else (int(op[0][i]), int(op[1][i]))
                new = (int(tr[0][i]), int(tr[1][i])) if dv[i] else (int(op[0][i]), int(op[1][i]))
                lst = preds[tid]
                j = lst.index((old[0], old[1], strands[i]))
                lst[j] = (new[0], new[1], strands[i])
            if v == "V0":
                same = all(sorted(base.get(k, [])) == sorted(preds.get(k, [])) for k in set(base) | set(preds))
                checks[-1]["v0_preds_equal_run"] = int(same and len(changed) == 0)
            row = dict(dataset=ds, arm=m.arm, group=m.group, org=m.org, input=m.input, variant=v,
                       n_changed=len(changed), n_changed_to_trim=int(dv[changed].sum()),
                       **score(truth, preds, negs, negseq))
            srows.append(row)
            cnt = collections.Counter()
            for i in changed:
                tid = tids[i]
                t = truth.get(tid)
                if t is None:
                    cls, typ, where = ("no_truth_neg" if tid in negs else "no_truth"), "-", "-"
                else:
                    typ = t[3]
                    t5, t3 = ends(*t[:3])
                    a3 = op[1][i] if strands[i] == "+" else op[0][i]
                    if t[2] != strands[i] or abs(a3 - t3) > TOL:
                        cls, where = "nostop_truth", "-"
                    else:
                        b0 = best_level(t, base.get(tid, []))
                        b1 = best_level(t, preds.get(tid, []))
                        atg5 = tr[0][i] if strands[i] == "+" else tr[1][i]
                        off = (t5 - atg5) if strands[i] == "+" else (atg5 - t5)  # >0: true start downstream of ATG
                        where = "at_atg" if abs(off) <= TOL else ("later" if off > 0 else "upstream")
                        if b0 != "exact" and b1 == "exact":
                            cls = "w2r"
                        elif b0 == "exact" and b1 != "exact":
                            cls = "r2w"
                        else:
                            cls = "w2w_" + where if b0 != "exact" else "r2r"
                cnt[(typ, cls)] += 1
                if v != "V0":
                    chg.append(dict(dataset=ds, variant=v, transcript_id=tid, orf_id=C.orf_id.values[i],
                                    truth_type=typ, cls=cls, where=where, seg=int(sg[i]), coding=cod[i], kozak=K[i],
                                    hcall=C.hcall.values[i], new_trim=int(dv[i])))
            for (typ, cls), n in sorted(cnt.items()):
                crows.append(dict(dataset=ds, arm=m.arm, group=m.group, variant=v, truth_type=typ, cls=cls, n=n))
        print(ds, "done", checks[-1], file=sys.stderr, flush=True)
    sfx = "" if not only else ".partial"
    pd.DataFrame(srows).to_csv(f"{odir}/scores{sfx}.tsv", sep="\t", index=False, float_format="%.6f")
    pd.DataFrame(crows).to_csv(f"{odir}/changes{sfx}.tsv", sep="\t", index=False)
    pd.DataFrame(checks).to_csv(f"{odir}/checks{sfx}.tsv", sep="\t", index=False)
    pd.DataFrame(chg).to_csv(f"{odir}/changed_candidates{sfx}.tsv.gz", sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
