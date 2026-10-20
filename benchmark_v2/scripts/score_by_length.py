#!/usr/bin/env python3
"""Length-stratified scoring for the minimum-ORF-length sweep.

usage: score_by_length.py <inputs_dir> <runs_dir> <tool> [<tool> ...] > by_length.tsv

For each <org>/<variant>/<tool>/preds.tsv (variants with negatives: *.wRand, i.i.d. random
sequences 'rand_*'; *.wLnc, GENCODE lncRNAs 'lnc_*'):
  truth rows      binned by true protein length (CDS nt / 3, less the stop when present);
                  sensitivity at the stop level and exact level (ends within 2 nt, as in
                  score_predictions.py)
  predictions     binned by predicted protein length; counted as TP (stop match on a real
                  transcript), FP on a real transcript (any other), or a call on a negative
Output: long format org, variant, tool, kind (truth|pred), bin, n, tp_stop, tp_exact,
fp_real, fp_negative. Bins: <50, 50-59, 60-69, ..., 90-99, 100-149, >=150 aa.
"""
import os
import sys

TOL = 2
BINS = [(0, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 100), (100, 150), (150, 10**9)]


def bin_label(aa):
    for lo, hi in BINS:
        if lo <= aa < hi:
            return f"<{hi}" if lo == 0 else (f">={lo}" if hi == 10**9 else f"{lo}-{hi - 1}")


def ends(s, e, strand):
    return (s, e) if strand == "+" else (e, s)


def main():
    inputs, runs, tools = sys.argv[1], sys.argv[2], sys.argv[3:]
    print("org\tvariant\ttool\tkind\tbin\tn\ttp_stop\ttp_exact\tfp_real\tfp_negative")
    for org in sorted(os.listdir(runs)):
        if not os.path.isdir(os.path.join(runs, org)) or org == "logs":
            continue
        for variant in sorted(os.listdir(os.path.join(runs, org))):
            if not variant.endswith(("wRand", "wLnc")):
                continue
            truth = {}
            with open(os.path.join(inputs, org, variant + ".truth.tsv")) as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                    s, e = int(r["cds_start"]), int(r["cds_end"])
                    has_stop = r["type"] in ("complete", "5prime_partial")
                    aa = (e - s + 1) // 3 - (1 if has_stop else 0)
                    truth[r["transcript_id"]] = (s, e, r["strand"], aa)
            for tool in tools:
                pf = os.path.join(runs, org, variant, tool, "preds.tsv")
                if not os.path.exists(pf):
                    continue
                tstat = {b: [0, 0, 0] for b in map(lambda x: bin_label(x[0]), BINS)}
                pstat = {b: [0, 0, 0, 0] for b in tstat}
                credited = set()
                with open(pf) as fh:
                    for line in fh:
                        if not line.strip() or line.startswith(("#", "transcript_id\t")):
                            continue
                        x = line.rstrip("\n").split("\t")
                        tid, s, e, strand = x[0], int(x[1]), int(x[2]), x[3]
                        b = bin_label((e - s + 1) // 3)
                        pstat[b][0] += 1
                        t = truth.get(tid)
                        if tid.startswith(("rand_", "lnc_")):
                            pstat[b][3] += 1
                            continue
                        if t and t[2] == strand and tid not in credited:
                            p5, p3 = ends(s, e, strand)
                            t5, t3 = ends(t[0], t[1], t[2])
                            if abs(p3 - t3) <= TOL:
                                credited.add(tid)
                                pstat[b][1] += 1
                                tb = bin_label(t[3])
                                tstat[tb][1] += 1
                                if abs(p5 - t5) <= TOL:
                                    tstat[tb][2] += 1
                                continue
                        pstat[b][2] += 1
                for t in truth.values():
                    tstat[bin_label(t[3])][0] += 1
                for b, (n, ts, te) in tstat.items():
                    print(f"{org}\t{variant}\t{tool}\ttruth\t{b}\t{n}\t{ts}\t{te}\t0\t0")
                for b, (n, tp, fpr, fpx) in pstat.items():
                    print(f"{org}\t{variant}\t{tool}\tpred\t{b}\t{n}\t{tp}\t0\t{fpr}\t{fpx}")


if __name__ == "__main__":
    main()
