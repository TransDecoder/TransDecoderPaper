#!/usr/bin/env python3
"""E36 step 3: score the selected alpha-form configurations on every dataset (benchmark in-sample; simulated
assemblies and GM12878 held out: never used for tuning).

usage: python 07_evaluate_v1.py <manifest.tsv> <labels_dir> <selected.tsv> <out_dir>

Configurations: V0 (const a0 = 0.4) and every row of selected.tsv (form, mode = free / accept, tuned on arm X with
all 8 organisms). Each configuration is applied in both arms (ab initio rule; with the assist, the homology call
overrides where it fires). Scores: overall exact F1, per type f1 and f1w (with negatives), tp; GM12878: the E29
statistic (exact given stop for MANE 5' partials and complete references).
Changed decisions (vs the run) by oracle class: w2r (now exact), r2w (exact lost), w2w (stop matched, neither exact),
masked (another prediction already exact), nostop (stop is not the truth stop), no_truth / no_truth_neg.
Writes scores.tsv, changes.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abl_common_v1 as ac  # noqa: E402


def parse(ps):
    return {k: float(v) for k, v in (x.split("=") for x in ps.split(";"))}


def main():
    manifest, ldir, selp, odir = sys.argv[1:5]
    os.makedirs(odir, exist_ok=True)
    man = pd.read_csv(manifest, sep="\t")
    sel = pd.read_csv(selp, sep="\t")
    sel = sel[sel.params != "none_acceptable"]
    cfgs = [("V0", "const", "-", "-", "a0=0.4")]
    for _, s in sel.iterrows():
        if s.form == "const" and s.params == "a0=0.4":
            continue
        cfgs.append((f"{s.form}.{s['mode']}.{'ab' if s.arm == 'abinitio' else 'H'}", s.form, s["mode"], s.arm, s.params))
    srows, crows = [], []
    for _, m in man.iterrows():
        lab, const = ac.load(ldir, m.dataset)
        fs = ac.FastScorer(lab, const)
        hom = m.arm == "homology"
        run = lab.trim_run.values.astype(bool)
        inc = lab.on_excl.values == 0
        orc, oth = lab.oracle.values, lab.other_exact.values.astype(bool)
        tt = lab.truth_type.values
        for cid, form, mode, tarm, ps in cfgs:
            d = ac.decide(lab, ac.alpha_form(form, parse(ps), lab), hom)
            r = fs.score(d)
            r.update(dataset=m.dataset, arm=m.arm, group=m.group, org=m.org, input=m.input, config=cid, form=form,
                     mode=mode, tuned_arm=tarm, params=ps)
            if m.group == "gm12878":
                for t in ("5prime_partial", "complete"):
                    r[f"e29_exact_given_stop_{t}"] = r[f"tp_{t}"] / const[f"tpstop_{t}"]
            srows.append(r)
            ch = (d != run) & inc
            to_right = np.where(orc == "keep", ~d, d)  # new decision is the exact one
            cls = np.select([orc == "no_truth", orc == "nostop", oth, orc == "neither", to_right],
                            ["no_truth", "nostop", "masked", "w2w", "w2r"], "r2w")
            cls = np.where((orc == "no_truth") & (lab.is_neg.values == 1), "no_truth_neg", cls)
            if ch.any():
                c = pd.DataFrame(dict(truth_type=tt[ch], cls=cls[ch], to_trim=d[ch])).value_counts().reset_index()
                c.columns = ["truth_type", "cls", "to_trim", "n"]
                for _, x in c.iterrows():
                    crows.append(dict(dataset=m.dataset, arm=m.arm, group=m.group, config=cid, truth_type=x.truth_type,
                                      cls=x.cls, to_trim=int(x.to_trim), n=int(x.n)))
        print(m.dataset, "done", file=sys.stderr, flush=True)
    S = pd.DataFrame(srows)
    v0 = S[S.config == "V0"].set_index("dataset")
    for c in ["f1_exact"] + [f"f1w_{t}" for t in ac.TYPES] + [f"f1_{t}" for t in ac.TYPES]:
        S["d_" + c] = 100 * (S[c] - S.dataset.map(v0[c]))
    S.to_csv(f"{odir}/scores.tsv", sep="\t", index=False, float_format="%.6f")
    pd.DataFrame(crows).to_csv(f"{odir}/changes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
