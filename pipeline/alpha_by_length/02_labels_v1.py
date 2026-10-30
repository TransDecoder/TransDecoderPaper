#!/usr/bin/env python3
"""E36: per-candidate oracle labels and per-dataset constants for fast rescoring.

usage: python 02_labels_v1.py <manifest.tsv> <out_dir> [dataset ...]

For every manifest row: the candidate table (E35 01_candidates_v1.py output), the run's preds.tsv (excluded contigs
removed), the truth, negatives. Per candidate (5'-open ORF with a first in-frame ATG):
  tx_len (L), orf_len (O, open ORF incl. u), seg (|u|), coding, kozak, hcall, trim_run (run's decision)
  truth_type   type of the transcript's truth CDS ('-' none)
  stop_match   candidate stop = truth stop (2 nt, same strand)
  e_keep       keeping the 5' end open gives the exact CDS
  e_trim       trimming to the first ATG gives the exact CDS
  oracle       keep / trim / neither (stop-matched, neither exact) / nostop / no_truth
  other_exact  another prediction on the transcript already matches the truth exactly
  is_neg, on_excl, orf_has_stop (needed to type predictions on negative contigs)
  where        true start relative to the first ATG (at_atg / later / upstream), as E35
Constants (V0 = the run): n_truth, n_pred, tp_exact, per type n, npred, tp, negpred (E35 score()), tpstop per type.
"""
import gzip
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
spec = importlib.util.spec_from_file_location("e35", f"{P}/TransDecoderPaper/pipeline/kozak_floor/02_evaluate_v1.py")
e35 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e35)
TOL, ends, read_preds, read_truth = e35.TOL, e35.ends, e35.read_preds, e35.read_truth
STOPS = e35.STOPS
COMP = e35.COMP


def main():
    manifest, odir = sys.argv[1:3]
    only = set(sys.argv[3:])
    os.makedirs(odir, exist_ok=True)
    man = pd.read_csv(manifest, sep="\t")
    for _, m in man.iterrows():
        ds = m.dataset
        if only and ds not in only:
            continue
        C = pd.read_csv(m.cand, sep="\t", keep_default_na=False, dtype={"coding": float, "kozak": float})
        excl = set() if m.excluded == "-" else {l.split("\t")[0].rstrip("\n") for l in open(m.excluded)}
        negs = set() if m.negatives == "-" else {l.strip() for l in open(m.negatives) if l.strip()}
        negs -= excl
        truth = read_truth(m.truth)
        base = {k: v for k, v in read_preds(os.path.join(m.run_dir, "preds.tsv")).items() if k not in excl}
        negseq = e35.read_fasta(os.path.join(m.run_dir, "input.fa"), negs)
        negs = {t for t in negs if t in negseq}
        # constants (V0)
        row = e35.score(truth, base, negs, negseq)
        for t in e35.TYPES:
            row[f"tpstop_{t}"] = sum(1 for tid, tr_ in truth.items() if tr_[3] == t and
                                     e35.best_level(tr_, base.get(tid, [])) is not None)
        row["dataset"] = ds
        L, o5, ol, sg = C.tx_len.values, C.orf5.values, C.orf_len.values, C.seg.values
        plus = (C.strand == "+").values
        op = np.where(plus, o5 + 1, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5)
        tr = np.where(plus, o5 + 1 + sg, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5 - sg)
        trim_run = C.trim_run.values.astype(bool)
        n = len(C)
        tt = np.full(n, "-", dtype=object)
        stop_m = np.zeros(n, int); ek = np.zeros(n, int); et = np.zeros(n, int); oth = np.zeros(n, int)
        where = np.full(n, "-", dtype=object)
        isneg = np.zeros(n, int); onex = np.zeros(n, int); hstop = np.full(n, -1)
        n_multi = 0
        seen = {}
        for i in range(n):
            tid, s = C.transcript_id.values[i], C.strand.values[i]
            seen[tid] = seen.get(tid, 0) + 1
            if tid in excl:
                onex[i] = 1
                continue
            cur = (int(tr[0][i]), int(tr[1][i]), s) if trim_run[i] else (int(op[0][i]), int(op[1][i]), s)
            lst = list(base.get(tid, []))
            lst.remove(cur)  # raises if the run's prediction is not there
            if tid in negs:
                isneg[i] = 1
                q = negseq[tid][op[0][i] - 1:op[1][i]]
                if s == "-":
                    q = q.translate(COMP)[::-1]
                hstop[i] = int(q[-3:] in STOPS)
            t = truth.get(tid)
            if t is None:
                continue
            tt[i] = t[3]
            oth[i] = int(e35.best_level(t, lst) == "exact")
            t5, t3 = ends(*t[:3])
            a3 = op[1][i] if s == "+" else op[0][i]
            if t[2] != s or abs(a3 - t3) > TOL:
                continue
            stop_m[i] = 1
            open5 = op[0][i] if s == "+" else op[1][i]
            atg5 = tr[0][i] if s == "+" else tr[1][i]
            ek[i] = int(abs(open5 - t5) <= TOL)
            et[i] = int(abs(atg5 - t5) <= TOL)
            off = (t5 - atg5) if s == "+" else (atg5 - t5)
            where[i] = "at_atg" if abs(off) <= TOL else ("later" if off > 0 else "upstream")
        n_multi = sum(1 for v in seen.values() if v > 1)
        oracle = np.where(tt == "-", "no_truth", np.where(stop_m == 0, "nostop",
                          np.where(ek == 1, "keep", np.where(et == 1, "trim", "neither"))))
        out = pd.DataFrame(dict(dataset=ds, transcript_id=C.transcript_id, orf_id=C.orf_id, strand=C.strand,
                                tx_len=L, orf5=o5, orf_len=ol, seg=sg, coding=C.coding, kozak=C.kozak,
                                hcall=C.hcall, trim_run=C.trim_run, truth_type=tt, stop_match=stop_m, e_keep=ek,
                                e_trim=et, oracle=oracle, other_exact=oth, where=where, is_neg=isneg, on_excl=onex,
                                orf_has_stop=hstop))
        out.to_csv(f"{odir}/{ds}.labels.tsv.gz", sep="\t", index=False, float_format="%.17g")
        row["n_cand"] = n
        row["n_tx_multi_cand"] = n_multi
        pd.DataFrame([row]).to_csv(f"{odir}/{ds}.const.tsv", sep="\t", index=False)
        print(ds, "candidates", n, "multi-candidate transcripts", n_multi, "both exact",
              int((ek & et).sum()), file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
