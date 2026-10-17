#!/usr/bin/env python3
"""Score every completed benchmark run and write one long-format table.

usage: score_all.py <inputs_dir> <runs_dir> > summary.tsv

Walks <runs_dir>/<org>/<variant>/<tool>/preds.tsv (plus the LongOrfs baselines that
the pytd runs leave as longest_orf.preds.tsv and all_orfs.preds.tsv) and scores each
against <inputs_dir>/<org>/<variant>.truth.tsv with score_predictions.py. Adds the
run time when available.
"""
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))


def main():
    inputs, runs = sys.argv[1], sys.argv[2]
    header_done = False
    for org in sorted(os.listdir(runs)):
        odir = os.path.join(runs, org)
        if org == "logs" or not os.path.isdir(odir):
            continue
        for variant in sorted(os.listdir(odir)):
            truth = os.path.join(inputs, org, variant + ".truth.tsv")
            for tool in sorted(os.listdir(os.path.join(odir, variant))):
                rdir = os.path.join(odir, variant, tool)
                cands = [(tool, "preds.tsv")]
                if tool == "pytd":
                    cands += [("longest_orf", "longest_orf.preds.tsv"), ("all_orfs", "all_orfs.preds.tsv")]
                for label, fn in cands:
                    pf = os.path.join(rdir, fn)
                    if not os.path.exists(pf):
                        continue
                    out = subprocess.run([sys.executable, os.path.join(here, "score_predictions.py"),
                                          "--truth", truth, "--preds", pf], capture_output=True, text=True, check=True)
                    hdr, vals = out.stdout.strip().split("\n")
                    rt = "NA"
                    rtf = os.path.join(rdir, "runtime_sec.txt")
                    if label == tool and os.path.exists(rtf):
                        rt = open(rtf).read().strip()
                    if not header_done:
                        print("org\tvariant\ttool\t" + hdr.split("\t", 1)[1] + "\truntime_sec")
                        header_done = True
                    print(f"{org}\t{variant}\t{label}\t" + vals.split("\t", 1)[1] + f"\t{rt}")


if __name__ == "__main__":
    sys.exit(main())
