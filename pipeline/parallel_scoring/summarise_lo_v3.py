#!/usr/bin/env python3
"""Parity and timing of the LongOrfs --cpu runs (pipeline/parallel_scoring/run_lo_v3.sh).

Parity: LongOrfs files (longest_orfs.gff3/cds/pep, base_freqs.dat, longorfs.params.json) of each run
compared byte for byte with the pinned release run (work/release_pin, PyTD d3892b5); pyTransdecoder runs
compared on all LongOrfs and Predict files.
Timing: LongOrfs wall time against the v0.3.0 (single-process LongOrfs) time of the same input in the
v2 series (work/fast_scoring, step_sec.tsv).

usage: summarise_lo_v3.py <out_prefix>     writes <out_prefix>.parity.tsv, .timing.tsv
"""
import filecmp
import glob
import os
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
W, REL, V2 = f"{P}/work/longorfs_cpu", f"{P}/work/release_pin", f"{P}/work/fast_scoring"
LO = ["longest_orfs.gff3", "longest_orfs.cds", "longest_orfs.pep", "base_freqs.dat", "longorfs.params.json"]
PRED = ["{t}.transdecoder.gff3", "{t}.transdecoder.bed", "{t}.transdecoder.pep", "{t}.transdecoder.cds",
        "{t}.transdecoder_dir/hexamer.scores", "{t}.transdecoder_dir/longest_orfs.cds.scores",
        "{t}.transdecoder_dir/decoy_qvalues.tsv", "{t}.transdecoder_dir/best_candidates.gff3",
        "{t}.transdecoder_dir/best_candidates.gff3.revised_starts.gff3",
        "{t}.transdecoder_dir/start_refinement.kozak.tsv"]


def compare(a, b, t, files):
    same, diff = 0, []
    for f in files:
        f = f.format(t=t)
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if os.path.exists(fa) and os.path.exists(fb) and filecmp.cmp(fa, fb, shallow=False):
            same += 1
        else:
            diff.append(os.path.basename(f))
    return same, diff


def v2_longorfs(org, inp):
    pat = (f"{V2}/gm/{inp}/rel.cpu*/step_sec.tsv" if org == "gm" else f"{V2}/bench/{org}/{inp}/rel.cpu*/step_sec.tsv")
    vals = []
    for st in glob.glob(pat):
        for line in open(st):
            k, v = line.split()
            if k == "longorfs":
                vals.append(float(v))
    return min(vals) if vals else None


def main(prefix):
    par, tim = [], {}
    for done in sorted(glob.glob(f"{W}/*/*/*/cpu*/DONE")):
        d = os.path.dirname(done)
        mode, org, inp, cpu = d.split(os.sep)[-4:]
        t = "transcripts.fa" if mode == "gm" else "input.fa"
        ref = f"{REL}/gm/{inp}/rel" if mode == "gm" else f"{REL}/bench/{org}/{inp}/rel"
        files = [f"{{t}}.transdecoder_dir/{f}" for f in LO] + (PRED if mode == "pipe" else [])
        same, diff = compare(d, ref, t, files)
        par.append((mode, f"{org}/{inp}" if mode != "gm" else f"gm/{inp}", cpu[3:], same, ",".join(diff) or "-"))
        if mode != "pipe":
            tim.setdefault((org, inp), {})[int(cpu[3:])] = float(open(f"{d}/runtime_sec.txt").read())
    with open(prefix + ".parity.tsv", "w") as fh:
        fh.write("mode\tinput\tcpu\tfiles_identical\tfiles_differing\n")
        for r in par:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(prefix + ".timing.tsv", "w") as fh:
        fh.write("input\tv030_longorfs_s\tcpu1_s\tcpu8_s\tspeedup_cpu8_vs_v030\n")
        for (org, inp), t in sorted(tim.items()):
            b = v2_longorfs(org, inp)
            f = lambda x: f"{x:.1f}" if x is not None else "NA"
            sp = f"{b / t[8]:.2f}" if b and 8 in t else "NA"
            fh.write(f"{org}/{inp}\t{f(b)}\t{f(t.get(1))}\t{f(t.get(8))}\t{sp}\n")
    print(f"{len(par)} comparisons, {sum(1 for r in par if r[4] != '-')} with differences")


if __name__ == "__main__":
    main(sys.argv[1])
