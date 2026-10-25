#!/usr/bin/env python3
"""Score the balanced mixed inputs fragments_mixed_t1..t3 (E28; make_mixed_inputs.py).

Writes, under <outdir> (default data/benchmark/final_v7/):
  mixed_scores.tsv
      one row per org x trial (variant) x tool; columns as summary.interim_20261025.tsv
      (score_predictions.py metrics + runtime_sec).
  five_prime_by_distance.tsv
      the 5' partials of the mixed inputs, per org x variant x tool x stratum x bin: references,
      stop- and exact-level true positives and sensitivities. On a 5' partial an exact match means
      the prediction starts at the fragment start (ORF left open, no start codon chosen); a
      stop-level match means the stop is right. Strata:
        dist_first_inframe_atg  nt from the fragment start to the first in-frame ATG:
                                0-29, 30-59, 60-119, 120-299, 300+, none;
        cut_frac                fraction of the reference CDS removed (cut_codon_5p / codons),
                                0.1-wide bins (cuts lie in 0-0.8).
Matching is that of benchmark_v2/scripts/score_predictions.py (ends within 2 nt, same strand, one
credit per reference). The longest-ORF baseline is the pytd_v7 run's longest_orf.preds.tsv. Tools
that were not run (T. thermophila for the standard-code-only tools) have no rows.

usage: score_mixed.py <inputs_dir> <runs_dir> [outdir]
"""
import os
import subprocess
import sys

T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, os.path.join(T, "benchmark_v2", "scripts"))
from score_predictions import TOL, ends, read_preds  # noqa: E402

VARIANTS = ["fragments_mixed_t1", "fragments_mixed_t2", "fragments_mixed_t3"]
TOOLS = ["pytd_v7", "perl_td_v7", "perl_td_v600", "prodigal", "gmst", "codan", "cpat", "cpc2", "rnasamba",
         "borf", "longest_orf"]
ORGS = ["human", "mouse", "dmel", "athal", "spombe", "tthermophila", "pfalciparum", "creinhardtii"]
DIST_BINS = [(0, 29, "0-29"), (30, 59, "30-59"), (60, 119, "60-119"), (120, 299, "120-299"), (300, 10**9, "300+")]
FRAC_BINS = [f"{i / 10:.1f}-{(i + 1) / 10:.1f}" for i in range(8)]


def pred_file(runs, org, variant, tool):
    if tool == "longest_orf":
        return os.path.join(runs, org, variant, "pytd_v7", "longest_orf.preds.tsv"), None
    d = os.path.join(runs, org, variant, tool)
    return os.path.join(d, "preds.tsv"), os.path.join(d, "runtime_sec.txt")


def summary_rows(inputs, runs):
    hdr, rows = None, []
    scorer = os.path.join(T, "benchmark_v2", "scripts", "score_predictions.py")
    for org in ORGS:
        for variant in VARIANTS:
            truth = os.path.join(inputs, org, variant + ".truth.tsv")
            for tool in TOOLS:
                pf, rtf = pred_file(runs, org, variant, tool)
                if not os.path.exists(pf):
                    continue
                out = subprocess.run([sys.executable, scorer, "--truth", truth, "--preds", pf],
                                     capture_output=True, text=True, check=True).stdout.strip().split("\n")
                hdr = "org\tvariant\ttool\t" + out[0].split("\t", 1)[1] + "\truntime_sec"
                rt = open(rtf).read().strip() if rtf and os.path.exists(rtf) else "NA"
                rows.append(f"{org}\t{variant}\t{tool}\t" + out[1].split("\t", 1)[1] + f"\t{rt}")
    return hdr, rows


def read_truth_rows(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]


def best_match(t, plist):
    best = None
    t5, t3 = ends(int(t["cds_start"]), int(t["cds_end"]), t["strand"])
    for p in plist:
        if p[2] != t["strand"]:
            continue
        p5, p3 = ends(*p)
        if abs(p3 - t3) <= TOL:
            lvl = "exact" if abs(p5 - t5) <= TOL else "stop"
            if best != "exact":
                best = lvl
    return best


def dist_bin(d):
    if d == "NA":
        return "none"
    return next(lab for lo, hi, lab in DIST_BINS if lo <= int(d) <= hi)


def frac_bin(t):
    cut = int(t["cut_codon_5p"])
    n_codons = cut + (int(t["cds_end"]) - int(t["cds_start"]) + 1) // 3
    return FRAC_BINS[min(7, int(10 * cut / n_codons))]


def by_distance(inputs, runs):
    cols = ["org", "variant", "tool", "stratum", "bin", "n_truth", "tp_stop", "tp_exact",
            "sensitivity_stop", "sensitivity_exact"]
    out = ["\t".join(cols)]
    order = {"dist_first_inframe_atg": [b[2] for b in DIST_BINS] + ["none"], "cut_frac": FRAC_BINS}
    for org in ORGS:
        for variant in VARIANTS:
            truth = [t for t in read_truth_rows(os.path.join(inputs, org, variant + ".truth.tsv"))
                     if t["type"] == "5prime_partial"]
            for tool in TOOLS:
                pf, _ = pred_file(runs, org, variant, tool)
                if not os.path.exists(pf):
                    continue
                preds = read_preds(pf)
                acc = {}
                for t in truth:
                    m = best_match(t, preds.get(t["transcript_id"], []))
                    for key in (("dist_first_inframe_atg", dist_bin(t["dist_first_inframe_atg"])),
                                ("cut_frac", frac_bin(t))):
                        c = acc.setdefault(key, [0, 0, 0])
                        c[0] += 1
                        c[1] += m is not None
                        c[2] += m == "exact"
                for stratum, bins in order.items():
                    for b in bins:
                        if (stratum, b) in acc:
                            n, ts, te = acc[(stratum, b)]
                            out.append("\t".join(map(str, (org, variant, tool, stratum, b, n, ts, te,
                                                           f"{ts / n:.4f}", f"{te / n:.4f}"))))
    return out


def main():
    inputs, runs = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(T, "data", "benchmark", "final_v7")
    hdr, rows = summary_rows(inputs, runs)
    with open(os.path.join(outdir, "mixed_scores.tsv"), "w") as fh:
        fh.write("\n".join([hdr] + rows) + "\n")
    with open(os.path.join(outdir, "five_prime_by_distance.tsv"), "w") as fh:
        fh.write("\n".join(by_distance(inputs, runs)) + "\n")


if __name__ == "__main__":
    sys.exit(main())
