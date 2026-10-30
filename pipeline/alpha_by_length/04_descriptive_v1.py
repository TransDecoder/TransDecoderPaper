#!/usr/bin/env python3
"""E36 step 1: per-candidate oracle answer (keep open = truncated vs trim = UTR) by L, O, O/L and |u|/L bins.

usage: python 04_descriptive_v1.py <manifest.tsv> <labels_dir> <out_dir>

Ab initio runs only (the candidate sets of the assist runs differ only by homology-driven selection). Groups:
bench_complete / bench_complete_fullutr / bench_fragments / bench_fragments_mixed (8 organisms pooled; mixed t1-t3
pooled; also per organism), each simulated assembly (human d1-d32, athal d1-d32), GM12878 Trinity (E29 MANE truth:
5' partials = keep, complete = trim; MANE-stop contigs only).
Per bin: n_keep (keeping open is exact: truncated CDS), n_trim (trimming to the first ATG is exact: UTR), n_neither
(stop-matched, neither exact, e.g. later true ATG), trunc_share = n_keep / (n_keep + n_trim), release_acc (share of
keep + trim decided right by the release rule), best_alpha (grid 0-3, maximises right decisions in the bin; ties to the
value closest to 0.4) and the right decisions at it.
end3 (oracle_by_bin.tsv): the candidate's 3' end from its stop-matched truth CDS, 'stop' (truth complete vs
5' partial) or 'open3' (no stop: truth 3' partial vs internal); 'all' = both plus unlabelled-type rows.
Writes oracle_by_bin.tsv (feature x bin x group), oracle_by_bin_org.tsv (benchmark per organism),
oracle_2d_L_OL.tsv (L x O/L), feature_medians.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abl_common_v1 as ac  # noqa: E402

GRID = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0])
BINS = {
    "L": ([0, 500, 1000, 1500, 2000, 3000, 5000, np.inf], "tx_len"),
    "O": ([0, 300, 600, 900, 1200, 1800, 3000, np.inf], "orf_len"),
    "OL": ([0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0001], "OL"),
    "uL": ([0, 0.02, 0.05, 0.1, 0.2, 0.4, 1.0001], "uL"),
    "u": ([0, 20, 50, 100, 200, 500, np.inf], "seg"),
}


def group_of(m):
    if m.group.startswith("bench"):
        return m.group
    if m.group == "gm12878":
        return "gm12878_trinity"
    return f"{m.group}_{m.input}"


def summarize(D, keys):
    rows = []
    for k, g in D.groupby(keys, observed=True, sort=True):
        kk = g[g.oracle.isin(["keep", "trim"])]
        nk, nt = int((kk.oracle == "keep").sum()), int((kk.oracle == "trim").sum())
        a = kk.astar.values
        keep = (kk.oracle == "keep").values
        right = np.array([((a >= x) & keep).sum() + ((a < x) & ~keep).sum() for x in GRID])
        best = right.max() if len(right) else 0
        cand = GRID[right == best]
        ba = cand[np.argmin(np.abs(cand - 0.4))] if len(kk) else np.nan
        rel = ((a >= 0.4) & keep).sum() + ((a < 0.4) & ~keep).sum()
        rows.append(dict(zip(keys, k if isinstance(k, tuple) else (k,)), n_keep=nk, n_trim=nt,
                         n_neither=int((g.oracle == "neither").sum()),
                         trunc_share=nk / (nk + nt) if nk + nt else np.nan,
                         release_acc=rel / (nk + nt) if nk + nt else np.nan, best_alpha=ba,
                         right_best=int(best), right_release=int(rel)))
    return pd.DataFrame(rows)


def main():
    manifest, ldir, odir = sys.argv[1:4]
    man = pd.read_csv(manifest, sep="\t")
    man = man[man.arm == "abinitio"]
    parts = []
    for _, m in man.iterrows():
        lab, _ = ac.load(ldir, m.dataset)
        lab = lab[lab.on_excl == 0].copy()
        lab["astar"] = ac.astar(lab)
        lab["OL"] = lab.orf_len / lab.tx_len
        lab["uL"] = lab.seg / lab.tx_len
        lab["grp"] = group_of(m)
        lab["org"] = m.org
        tt = lab.truth_type
        lab["end3"] = np.where(tt.isin(["complete", "5prime_partial"]), "stop",
                               np.where(tt.isin(["3prime_partial", "internal"]), "open3", "-"))
        parts.append(lab[["grp", "org", "dataset", "tx_len", "orf_len", "seg", "OL", "uL", "astar", "oracle",
                          "truth_type", "end3"]])
    D = pd.concat(parts, ignore_index=True)
    out1, out2 = [], []
    for f, (edges, col) in BINS.items():
        D["bin"] = pd.cut(D[col], edges, right=False)
        s = summarize(D, ["grp", "bin"])
        s.insert(1, "feature", f)
        s.insert(2, "end3", "all")
        out1.append(s)
        s = summarize(D[D.end3 != "-"], ["grp", "end3", "bin"])
        s.insert(1, "feature", f)
        out1.append(s)
        b = D[D.grp.str.startswith("bench")]
        s = summarize(b, ["grp", "org", "bin"])
        s.insert(2, "feature", f)
        out2.append(s)
        tot = summarize(D, ["grp"])
        tot.insert(1, "feature", f)
        tot.insert(2, "end3", "all")
        tot["bin"] = "all"
        out1.append(tot)
    pd.concat(out1).to_csv(f"{odir}/oracle_by_bin.tsv", sep="\t", index=False, float_format="%.4f")
    pd.concat(out2).to_csv(f"{odir}/oracle_by_bin_org.tsv", sep="\t", index=False, float_format="%.4f")
    D["Lbin"] = pd.cut(D.tx_len, [0, 1000, 2000, 4000, np.inf], right=False)
    D["OLbin"] = pd.cut(D.OL, [0, 0.5, 0.7, 0.9, 1.0001], right=False)
    summarize(D, ["grp", "Lbin", "OLbin"]).to_csv(f"{odir}/oracle_2d_L_OL.tsv", sep="\t", index=False,
                                                  float_format="%.4f")
    kk = D[D.oracle.isin(["keep", "trim"])]
    med = kk.groupby(["grp", "end3", "oracle"]).agg(n=("tx_len", "size"), L_med=("tx_len", "median"),
                                            O_med=("orf_len", "median"), OL_med=("OL", "median"),
                                            u_med=("seg", "median"), uL_med=("uL", "median"),
                                            OL_ge095=("OL", lambda x: float((x >= 0.95).mean())))
    med.reset_index().to_csv(f"{odir}/feature_medians.tsv", sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
