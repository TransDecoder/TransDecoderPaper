#!/usr/bin/env python3
"""Parity and timing of the parallel-scoring runs (pipeline/parallel_scoring/run_one_v1.sh).

Parity: every Predict output (final gff3/bed/pep/cds and the .transdecoder_dir files) of each run is
compared byte for byte with (a) the --cpu 1 run of the same input and (b) the release run
(work/release_pin, PyTD d3892b5, label/cfg rel).
Timing: Predict steps from the predict.log timestamps; LongOrfs from step_sec.tsv.

usage: summarise_v1.py <work/parallel_scoring> <work/release_pin> <out_prefix>
writes <out_prefix>.parity.tsv and <out_prefix>.timing.tsv
"""
import datetime
import filecmp
import glob
import os
import re
import sys

FILES = ["{t}.transdecoder.gff3", "{t}.transdecoder.bed", "{t}.transdecoder.pep", "{t}.transdecoder.cds",
         "{t}.transdecoder_dir/longest_orfs.gff3", "{t}.transdecoder_dir/longest_orfs.cds",
         "{t}.transdecoder_dir/longest_orfs.pep", "{t}.transdecoder_dir/hexamer.scores",
         "{t}.transdecoder_dir/longest_orfs.cds.scores", "{t}.transdecoder_dir/decoy_qvalues.tsv",
         "{t}.transdecoder_dir/best_candidates.gff3",
         "{t}.transdecoder_dir/best_candidates.gff3.revised_starts.gff3",
         "{t}.transdecoder_dir/start_refinement.kozak.tsv"]


def ts(line):
    return datetime.datetime.strptime(line[:23], "%Y-%m-%d %H:%M:%S,%f")


def steps(d):
    lines = [l for l in open(os.path.join(d, "predict.log")) if re.match(r"\d{4}-", l)]

    def at(pat):
        for l in lines:
            if pat in l:
                return ts(l)
    t0, end = ts(lines[0]), ts(lines[-1])
    sc, sel = at("Scoring all candidate ORFs"), at("Selecting best ORFs per transcript")
    gc, dec = at("Transcriptome GC") or at("Parsing CDS scores") or sel, at("Decoys:")
    ref, gen = at("Refining start"), at("Generating final")
    sec = lambda a, b: (b - a).total_seconds()
    st = dict(l.split("\t") for l in open(os.path.join(d, "step_sec.tsv")).read().split("\n") if l)
    return {"total": float(open(os.path.join(d, "runtime_sec.txt")).read()),
            "longorfs": float(st["longorfs"]), "predict": float(st["predict"]),
            "score": sec(sc, sel), "decoy": sec(gc, dec), "train": sec(t0, sc),
            "start": sec(ref, gen), "output": sec(gen, end)}


def compare(a, ta, b, tb):
    """(n identical, list of differing/missing files)"""
    same, diff = 0, []
    for f in FILES:
        fa, fb = os.path.join(a, f.format(t=ta)), os.path.join(b, f.format(t=tb))
        if not os.path.exists(fa) and not os.path.exists(fb):
            continue
        if os.path.exists(fa) and os.path.exists(fb) and filecmp.cmp(fa, fb, shallow=False):
            same += 1
        else:
            diff.append(f.split("/")[-1].format(t=""))
    return same, diff


def main(par, rel, prefix):
    runs = sorted(glob.glob(os.path.join(par, "*", "*", "cpu*", "DONE")))
    par_rows, time_rows = [], []
    for done in runs:
        d = os.path.dirname(done)
        org, inp, cpu = d.split(os.sep)[-3:]
        t = "transcripts.fa" if org == "gm" else "input.fa"
        base = os.path.join(par, org, inp, "cpu1")
        r = os.path.join(rel, "gm", inp, "rel") if org == "gm" else os.path.join(rel, "bench", org, inp, "rel")
        for ref_name, ref in (("cpu1", base), ("release", r)):
            if ref == d or not os.path.exists(os.path.join(ref, "DONE")):
                continue
            same, diff = compare(d, t, ref, t)
            par_rows.append((org, inp, cpu, ref_name, same, ",".join(diff) or "-"))
        s = steps(d)
        time_rows.append((org, inp, int(cpu[3:]), s))
    with open(prefix + ".parity.tsv", "w") as fh:
        fh.write("org\tinput\tcpu\tcompared_to\tfiles_identical\tfiles_differing\n")
        for row in par_rows:
            fh.write("\t".join(map(str, row)) + "\n")
    keys = ["total", "longorfs", "predict", "train", "score", "decoy", "start", "output"]
    one = {(o, i): s for o, i, c, s in time_rows if c == 1}
    with open(prefix + ".timing.tsv", "w") as fh:
        fh.write("org\tinput\tcpu\t" + "\t".join(k + "_s" for k in keys) +
                 "\tspeedup_total\tspeedup_score\tspeedup_decoy\n")
        for o, i, c, s in sorted(time_rows, key=lambda x: (x[0], x[1], x[2])):
            b = one.get((o, i))
            sp = [f"{b[k] / s[k]:.2f}" if b and s[k] > 0 else "NA" for k in ("total", "score", "decoy")]
            fh.write(f"{o}\t{i}\t{c}\t" + "\t".join(f"{s[k]:.1f}" for k in keys) + "\t" + "\t".join(sp) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
