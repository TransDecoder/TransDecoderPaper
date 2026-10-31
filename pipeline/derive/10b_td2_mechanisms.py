#!/usr/bin/env python3
"""E39: two descriptive tallies behind the TD2 sentences of the Results.

td2_fragments_by_cds_length.tsv  stranded fragments: stop-level sensitivity of TD2 and TransDecoder for
    reference CDSs of at least / fewer than 90 codons (270 nt; TD2's LongOrfs minimum). Same end rule
    as score_predictions.py (stranded input: predicted end within 2 nt of cds_end).
td2_random_calls_by_rule.tsv     complete.wRand: TD2 predictions on random sequences, and those on
    sequences with an ORF above the PSAURON cutoff (0.5; input/psauron_score.csv); the others can only
    have passed by the long-ORF length rule.

usage (from the paper repo root): 10b_td2_mechanisms.py
"""
import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
R = f"{P}/work/bench_v3/runs"
ORGS = "pfalciparum tthermophila spombe athal dmel mouse human creinhardtii".split()
F = "data/benchmark/final_v7/"
rows, rnd = [], []
for o in ORGS:
    t = pd.read_csv(f"{P}/data/bench_v3/inputs/{o}/fragments.truth.tsv", sep="\t")
    t["L"] = t.cds_end - t.cds_start + 1
    for tool in ("td2", "pytd_v7"):
        p = pd.read_csv(f"{R}/{o}/fragments/{tool}/preds.tsv", sep="\t")
        m = t.merge(p, on="transcript_id", how="left", suffixes=("", "_p"))
        hit = (m.end - m.cds_end).abs() <= 2
        hits = hit.groupby(m.transcript_id).any()
        ge = t.set_index("transcript_id").L >= 270
        rows.append(dict(org=o, tool=tool, n_ge90=int(ge.sum()), n_lt90=int((~ge).sum()), frac_lt90=round((~ge).mean(), 4),
                         sens_stop_ge90=round(hits[ge.index][ge].mean(), 4), sens_stop_lt90=round(hits[ge.index][~ge].mean(), 4)))
    d = f"{R}/{o}/complete.wRand/td2"
    pr = pd.read_csv(f"{d}/preds.tsv", sep="\t")
    pr = pr[pr.transcript_id.str.startswith("rand_")]
    ps = pd.read_csv(f"{d}/input/psauron_score.csv", skiprows=2)
    ps["tx"] = ps.description.str.split(" ").str[0].str.rsplit(".", n=1).str[0]
    best = ps.groupby("tx")["in-frame_score"].max()
    above = int((pr.transcript_id.map(best) > 0.5).sum())
    n_seq = sum(1 for l in open(f"{P}/data/bench_v3/inputs/{o}/complete.wRand.fa") if l.startswith(">rand_"))
    rnd.append(dict(org=o, n_random=n_seq, n_called=len(pr), frac_called=round(len(pr) / n_seq, 4),
                    called_above_psauron_cutoff=above, called_length_only=len(pr) - above))
pd.DataFrame(rows).to_csv(F + "td2_fragments_by_cds_length.tsv", sep="\t", index=False)
pd.DataFrame(rnd).to_csv(F + "td2_random_calls_by_rule.tsv", sep="\t", index=False)
