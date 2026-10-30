#!/usr/bin/env python3
"""Supplementary tables assembled from committed data (E37).

Builds, from tables already in this repository (no new runs):
  perl_python_sample_sets.tsv   Perl vs Python on the five TransDecoder sample data sets, at the
                                release commits and at v6.0.0 / first PyTD port (E08)
  benchmark_per_input.tsv       every benchmark input x method that the manuscript reports
                                (E07, E25, E28); supplied as an additional data file
  start_choice_legacy.tsv       exact matches on complete inputs: first in-frame ATG vs a later ATG
                                chosen by the PWM, v6.0.0 vs the scale-corrected PWM, release (E09)
  long_orf_length_by_gc.tsv     L* by G+C, with the 2017 table it replaces (E10)
  benchmark_truth_sets.tsv      annotation file, transcripts filtered and UTRs added per organism (E07)
and prints the LaTeX table bodies used in manuscript/supplement.tex.

usage (from the paper repo root):
  $P/env/td/bin/python pipeline/derive/09_supplementary_tables.py > data/supplementary_tables/latex_bodies.txt
"""
import os

import pandas as pd

OUT = "data/supplementary_tables"
ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
ORG_LABEL = {"pfalciparum": "P. falciparum", "tthermophila": "T. thermophila", "spombe": "S. pombe",
             "athal": "A. thaliana", "dmel": "D. melanogaster", "mouse": "M. musculus",
             "human": "H. sapiens", "creinhardtii": "C. reinhardtii"}
NO_CILIATE_CODE = {"gmst", "codan", "cpc2", "rnasamba", "borf"}
METHOD = {"pytd_v7": "TransDecoder release (PyTransDecoder v0.2.0)",
          "perl_td_v7": "TransDecoder release (Perl v7.0.0)",
          "pytd_v7_homology": "TransDecoder release + homology (PyTransDecoder, DIAMOND/Swiss-Prot and Pfam)",
          "perl_td_v600": "TransDecoder v6.0.0", "prodigal": "Prodigal 2.6.3",
          "gmst": "GeneMarkS-T 5.1", "codan": "CodAn 1.1", "cpat": "CPAT 3.0.5", "cpc2": "CPC2 1.0.1",
          "rnasamba": "RNAsamba 0.2.5", "borf": "Borf 1.2", "longest_orf": "longest ORF baseline",
          "angel": "ANGEL 3.0 (not benchmarked further)"}
COLS = ["n_truth", "n_pred", "n_pred_random", "tp_stop", "tp_exact", "fp", "fp_random", "fn",
        "sensitivity_stop", "precision_stop", "f1_stop", "sensitivity_exact", "precision_exact",
        "f1_exact", "runtime_sec"]


def tsv(path):
    return pd.read_csv(path, sep="\t")


def fmt_int(x):
    return f"{int(round(x)):,}".replace(",", "{,}")


def perl_python():
    rows = []
    for run, f in [("release (Perl 41fc5ee, PyTD 3a793ff)", "release_3a793ff_41fc5ee"),
                   ("v6.0.0 (Perl 66d4712) vs first PyTD port (43f4293)", "baseline_v6.0.0")]:
        d = tsv(f"data/regression/{f}.summary.tsv").pivot(index="dataset", columns="metric", values="value")
        for ds in ["trinity", "pasa", "stringtie", "supertx", "cufflinks"]:
            r = d.loc[ds]
            rows.append({"comparison": run, "sample_set": ds,
                         **{k: int(float(r[k])) for k in ["perl_orfs", "py_orfs", "exact_match",
                                                          "same_stop_diff_start", "perl_only", "py_only"]},
                         "perl_runtime_sec": round(float(r["perl_runtime_sec"]), 1),
                         "py_runtime_sec": round(float(r["py_runtime_sec"]), 1)})
    t = pd.DataFrame(rows)
    t.to_csv(f"{OUT}/perl_python_sample_sets.tsv", sep="\t", index=False)
    name = {"trinity": "Trinity", "pasa": "PASA", "stringtie": "StringTie", "supertx": "SuperTranscripts",
            "cufflinks": "Cufflinks"}
    print("% perl_python_sample_sets.tsv")
    for comp, g in t.groupby("comparison", sort=False):
        print(f"% {comp}")
        for _, r in g.iterrows():
            print(f"{name[r.sample_set]} & {fmt_int(r.perl_orfs)} & {fmt_int(r.py_orfs)} & "
                  f"{fmt_int(r.exact_match)} & {r.same_stop_diff_start} & {r.perl_only} & {r.py_only} & "
                  f"{r.perl_runtime_sec:.0f} & {r.py_runtime_sec:.0f} \\\\")
    return t


def per_input():
    summ = tsv("data/benchmark/summary.interim_20261025.tsv")
    rel = tsv("data/benchmark/final_v7/all_scores.tsv")
    mixed = tsv("data/benchmark/final_v7/mixed_scores.tsv")
    td = ["pytd_v7", "perl_td_v7", "pytd_v7_homology"]
    comp = [t for t in METHOD if t not in td]
    d = pd.concat([rel[rel.tool.isin(td)], summ[summ.tool.isin(comp)], mixed], ignore_index=True)
    d = d.drop_duplicates(["org", "variant", "tool"], keep="first")
    d["method"] = d.tool.map(METHOD)
    d["genetic_code_supported"] = ~((d.org == "tthermophila") & d.tool.isin(NO_CILIATE_CODE))
    d["genetic_code_supported"] = d.genetic_code_supported.map({True: "yes", False: "no"})
    d["organism"] = d.org.map(ORG_LABEL)
    order = {t: i for i, t in enumerate(METHOD)}
    d["_o"] = d.org.map({o: i for i, o in enumerate(ORGS)})
    d["_t"] = d.tool.map(order)
    d = d.sort_values(["_o", "variant", "_t"])
    d["runtime_sec"] = d.runtime_sec.round(1)
    out = d[["organism", "variant", "tool", "method", "genetic_code_supported"] + COLS]
    out = out.rename(columns={"variant": "input", "tool": "run_label"})
    out.to_csv(f"{OUT}/benchmark_per_input.tsv", sep="\t", index=False)
    print(f"% benchmark_per_input.tsv: {len(out)} rows, {out.input.nunique()} input types, "
          f"{out.organism.nunique()} organisms, {out.run_label.nunique()} runs")
    print("% rows per run: " + ", ".join(f"{k} {v}" for k, v in out.run_label.value_counts().sort_index().items()))
    return out


def start_choice():
    summ = tsv("data/benchmark/summary.tsv")
    rel = tsv("data/benchmark/final_v7/all_scores.tsv")
    s = summ[summ.variant == "complete"]
    r = rel[(rel.variant == "complete") & (rel.tool == "pytd_v7")]
    tools = ["pytd_dev", "pytd_dev_pwm", "perl_td_v600", "pytd"]
    rows = []
    for o in ORGS:
        row = {"org": o, "n_truth": int(s[s.org == o].n_truth.iloc[0])}
        for t in tools:
            x = s[(s.org == o) & (s.tool == t)].iloc[0]
            row[f"{t}_tp_exact"] = int(x.tp_exact)
            row[f"{t}_f1_exact"] = x.f1_exact
        x = r[r.org == o].iloc[0]
        row["pytd_v7_tp_exact"] = int(x.tp_exact)
        row["pytd_v7_f1_exact"] = x.f1_exact
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(f"{OUT}/start_choice_legacy.tsv", sep="\t", index=False)
    print("% start_choice_legacy.tsv (exact matches; first ATG, PWM later ATG, v6.0.0, PWM scale corrected, release)")
    for _, x in t.iterrows():
        cells = [fmt_int(x[f"{k}_tp_exact"]) for k in ["pytd_dev", "pytd_dev_pwm", "perl_td_v600", "pytd", "pytd_v7"]]
        print(f"\\textit{{{ORG_LABEL[x.org]}}} & {fmt_int(x.n_truth)} & " + " & ".join(cells) + " \\\\")
    med = {k: t[f"{k}_f1_exact"].median() for k in ["pytd_dev", "pytd_dev_pwm", "perl_td_v600", "pytd", "pytd_v7"]}
    print("median exact F1 & & " + " & ".join(f"{v:.3f}" for v in med.values()) + " \\\\")
    seeds = tsv("data/start_refinement/seed_start_changes.tsv")
    print("% seed_start_changes.tsv")
    for ds, g in seeds.groupby("dataset", sort=False):
        print(f"{ds} & {fmt_int(g.n_orfs.iloc[0])} & " + " & ".join(str(v) for v in g.same_stop_diff_start)
              + " & " + " & ".join(str(v) for v in g.stop_not_shared) + " \\\\")
    return t


def long_orf():
    thr = tsv("data/long_orf_threshold/final_thresholds.tsv")
    old = tsv("data/long_orf_threshold/table_2017.tsv")[["gc", "q0.999"]].rename(columns={"q0.999": "q0.999_2017"})
    t = thr.merge(old, on="gc")
    # lengths as used in the code of v5.0.0 to v6.0.0 (@GC_TO_MIN_LONG_ORF_LENGTH in TransDecoder.Predict,
    # v5.0.0 and util/TransDecoder.Predict at 66d4712); hand-rounded from the 2017 q0.999 column
    t["v5_v6_code"] = t.gc.map({25: 465, 30: 510, 35: 555, 40: 590, 45: 645, 50: 749, 55: 797, 60: 927,
                                65: 1086, 70: 1358, 75: 1743, 80: 2422})
    t.to_csv(f"{OUT}/long_orf_length_by_gc.tsv", sep="\t", index=False)
    print("% long_orf_length_by_gc.tsv")
    for _, x in t.iterrows():
        print(f"{int(x.gc)} & {fmt_int(x.n_orfs)} & {fmt_int(x['q0.95'])} & {fmt_int(x['q0.99'])} & {fmt_int(x['q0.999'])} & "
              f"{fmt_int(x['q0.9999'])} & {fmt_int(x['v5_v6_code'])} \\\\")
    return t


def truth_sets():
    src = tsv("benchmark_v2/annotation_downloads.tsv").set_index("org")
    rows = []
    for o in ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii"]:
        log = dict(line.split("\t") for line in open(f"benchmark_v2/truth_logs_v3/{o}.log").read().split("\n") if line)
        log = {k: int(v) for k, v in log.items()}
        excl = {k: v for k, v in log.items() if k not in
                {"transcripts_with_cds", "kept_genes", "kept_with_annotated_utr5", "kept_with_annotated_utr3",
                 "utr5_added", "utr3_added"}}
        rows.append({"org": o, "release": src.loc[o, "release"], "annotation_url": src.loc[o, "annotation_url"],
                     "genome_url": src.loc[o, "genome_url"],
                     "transcripts_with_cds": log["transcripts_with_cds"],
                     "excluded": ";".join(f"{k}={v}" for k, v in excl.items()),
                     "kept_genes": log["kept_genes"],
                     "annotated_utr5": log.get("kept_with_annotated_utr5", 0),
                     "utr5_added": log.get("utr5_added", 0),
                     "annotated_utr3": log.get("kept_with_annotated_utr3", 0),
                     "utr3_added": log.get("utr3_added", 0)})
    t = pd.DataFrame(rows)
    for e in ["utr5", "utr3"]:  # 'none': no annotated UTR and none added (utr*_source in the truth table)
        t[f"{e}_none"] = t.kept_genes - t[f"annotated_{e}"] - t[f"{e}_added"]
        assert (t[f"{e}_none"] >= 0).all(), e
    t.to_csv(f"{OUT}/benchmark_truth_sets.tsv", sep="\t", index=False)
    print("% benchmark_truth_sets.tsv")
    for _, x in t.iterrows():
        print(f"\\textit{{{ORG_LABEL[x.org]}}} & {fmt_int(x.transcripts_with_cds)} & {fmt_int(x.kept_genes)} & "
              f"{fmt_int(x.annotated_utr5)} & {fmt_int(x.utr5_added)} & {fmt_int(x.annotated_utr3)} & "
              f"{fmt_int(x.utr3_added)} \\\\")
    for _, x in t.iterrows():
        print(f"% {x.org}: {x.excluded}; no 5' UTR {x.utr5_none}, no 3' UTR {x.utr3_none}")
    return t


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    perl_python()
    per_input()
    start_choice()
    long_orf()
    truth_sets()
