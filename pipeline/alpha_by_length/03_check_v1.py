#!/usr/bin/env python3
"""E36 checks of the fast rescoring (abl_common_v1.FastScorer).

usage: python 03_check_v1.py <manifest.tsv> <labels_dir> <e35_scores.tsv> <out.tsv>

1. V0 (the run's own decisions) on every dataset = E35 scores.tsv V0 (90 datasets) / the constants.
2. E35 variants V1 and V3_8.5 recomputed with the fast scorer = E35's full-substitution scores (90 datasets x 2).
3. An alpha-form variant (logL a0 0.4, b 0.3) scored by fast scorer and by full substitution (E35 score() on the
   run's preds with every changed candidate's coordinates substituted) on 10 datasets incl. negatives.
Columns compared: f1_exact and per type f1 and f1w; tolerance 1e-9.
"""
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abl_common_v1 as ac  # noqa: E402

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
spec = importlib.util.spec_from_file_location("e35", f"{P}/TransDecoderPaper/pipeline/kozak_floor/02_evaluate_v1.py")
e35 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e35)

COLS = ["f1_exact"] + [f"f1_{t}" for t in ac.TYPES] + [f"f1w_{t}" for t in ac.TYPES]


def cmp(a, b):
    worst = 0.0
    for c in COLS:
        x, y = a.get(c, np.nan), b.get(c, np.nan)
        if (pd.isna(x) and pd.isna(y)) or (pd.isna(x) and y == 0) or (x == 0 and pd.isna(y)):
            continue
        worst = max(worst, abs(float(x) - float(y)))
    return worst


def full_score(m, lab, d):
    excl = set() if m.excluded == "-" else {l.split("\t")[0].rstrip("\n") for l in open(m.excluded)}
    negs = set() if m.negatives == "-" else {l.strip() for l in open(m.negatives) if l.strip()}
    negs -= excl
    truth = e35.read_truth(m.truth)
    base = {k: v for k, v in e35.read_preds(os.path.join(m.run_dir, "preds.tsv")).items() if k not in excl}
    negseq = e35.read_fasta(os.path.join(m.run_dir, "input.fa"), negs)
    negs = {t for t in negs if t in negseq}
    L, o5, ol, sg = lab.tx_len.values, lab.orf5.values, lab.orf_len.values, lab.seg.values
    plus = (lab.strand == "+").values
    op = np.where(plus, o5 + 1, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5)
    tr = np.where(plus, o5 + 1 + sg, L - (o5 + ol) + 1), np.where(plus, o5 + ol, L - o5 - sg)
    run = lab.trim_run.values.astype(bool)
    preds = {k: list(x) for k, x in base.items()}
    for i in np.where((d != run) & (lab.on_excl.values == 0))[0]:
        tid, s = lab.transcript_id.values[i], lab.strand.values[i]
        old = (int(tr[0][i]), int(tr[1][i]), s) if run[i] else (int(op[0][i]), int(op[1][i]), s)
        new = (int(tr[0][i]), int(tr[1][i]), s) if d[i] else (int(op[0][i]), int(op[1][i]), s)
        lst = preds[tid]
        lst[lst.index(old)] = new
    return e35.score(truth, preds, negs, negseq)


def main():
    manifest, ldir, e35s, out = sys.argv[1:5]
    man = pd.read_csv(manifest, sep="\t")
    S = pd.read_csv(e35s, sep="\t").set_index(["dataset", "variant"])
    rows = []
    for _, m in man.iterrows():
        lab, const = ac.load(ldir, m.dataset)
        fs = ac.FastScorer(lab, const)
        hom = m.arm == "homology"
        d0 = lab.trim_run.values.astype(bool)
        r0 = fs.score(d0)
        rows.append(dict(dataset=m.dataset, check="V0_fast_vs_const", max_abs_diff=cmp(r0, const),
                         n_changed=r0["n_changed"]))
        if (m.dataset, "V0") in S.index:
            rows.append(dict(dataset=m.dataset, check="V0_fast_vs_E35", max_abs_diff=cmp(r0, S.loc[(m.dataset, "V0")].to_dict()),
                             n_changed=r0["n_changed"]))
            # V0 decision rule = run
            dd = ac.decide(lab, 0.4, hom)
            rows.append(dict(dataset=m.dataset, check="V0_rule_vs_run_decisions", max_abs_diff=float((dd != d0).sum()),
                             n_changed=0))
            for v in ["V1", "V3_8.5"]:
                K = lab.kozak.values
                ab = lab.coding.values < 0.4 * lab.seg.values + e35.kozak_term(v, K, lab.coding.values)
                h = lab.hcall.values
                dv = np.where(h == "trim", True, np.where(h == "keep", False, ab)) if hom else ab
                r = fs.score(dv)
                rows.append(dict(dataset=m.dataset, check=f"{v}_fast_vs_E35_full", n_changed=r["n_changed"],
                                 max_abs_diff=cmp(r, S.loc[(m.dataset, v)].to_dict())))
        else:
            dd = ac.decide(lab, 0.4, hom)
            rows.append(dict(dataset=m.dataset, check="V0_rule_vs_run_decisions", max_abs_diff=float((dd != d0).sum()),
                             n_changed=0))
    full_sets = ["bench.human.complete", "bench.athal.fragments", "bench.mouse.fragments_mixed_t2",
                 "bench.pfalciparum.fragments.H", "sim.human.d1", "sim.human.d16", "sim.athal.d8.H", "sim.human.d32.H",
                 "gm12878.trinity", "gm12878.trinity.H"]
    for ds in full_sets:
        m = man.set_index("dataset").loc[ds]
        m = m.copy(); m["dataset"] = ds
        lab, const = ac.load(ldir, ds)
        fs = ac.FastScorer(lab, const)
        a = ac.alpha_form("logL", dict(a0=0.4, b=0.3), lab)
        d = ac.decide(lab, a, m.arm == "homology")
        r = fs.score(d)
        rf = full_score(m, lab, d)
        rows.append(dict(dataset=ds, check="logL_fast_vs_full_substitution", max_abs_diff=cmp(r, rf),
                         n_changed=r["n_changed"]))
        print(ds, "done", file=sys.stderr, flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(out, sep="\t", index=False, float_format="%.3g")
    print(R.groupby("check").agg(n=("dataset", "size"), max_diff=("max_abs_diff", "max"),
                                 changed=("n_changed", "sum")).to_string())


if __name__ == "__main__":
    main()
