#!/usr/bin/env python3
"""Parity and timing of the dev/faster-serial runs (pipeline/parallel_scoring/run_one_v2.sh).

Parity: each run's outputs are compared byte for byte with the release-code run of the same input and
setting (PyTD d3892b5 in work/release_pin; for homology start assist, the same-code-tree a5b24c3 runs
pytd_v7h_homology in work/bench_v3/runs), and the --cpu 1 and --cpu 8 runs of an input with each other.
Preset builds: preset files compared after gunzip; preset.json compared as JSON.
Timing: total, LongOrfs and Predict times, against the v1 runs (release code at --cpu 1 and the
dev/parallel-scoring branch at --cpu 8; work/parallel_scoring).

usage: summarise_v2.py <out_prefix>     writes <out_prefix>.parity.tsv, .timing.tsv
"""
import filecmp
import glob
import gzip
import json
import os
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
W, REL, V1 = f"{P}/work/fast_scoring", f"{P}/work/release_pin", f"{P}/work/parallel_scoring"
FILES = ["{t}.transdecoder.gff3", "{t}.transdecoder.bed", "{t}.transdecoder.pep", "{t}.transdecoder.cds",
         "{t}.transdecoder_dir/longest_orfs.gff3", "{t}.transdecoder_dir/longest_orfs.cds",
         "{t}.transdecoder_dir/longest_orfs.pep", "{t}.transdecoder_dir/base_freqs.dat",
         "{t}.transdecoder_dir/hexamer.scores", "{t}.transdecoder_dir/longest_orfs.cds.scores",
         "{t}.transdecoder_dir/decoy_qvalues.tsv", "{t}.transdecoder_dir/best_candidates.gff3",
         "{t}.transdecoder_dir/best_candidates.gff3.revised_starts.gff3",
         "{t}.transdecoder_dir/start_refinement.kozak.tsv", "{t}.transdecoder_dir/longorfs.params.json"]


def compare(a, b, t):
    same, diff, n_missing_both = 0, [], 0
    for f in FILES:
        fa, fb = os.path.join(a, f.format(t=t)), os.path.join(b, f.format(t=t))
        ea, eb = os.path.exists(fa), os.path.exists(fb)
        if not ea and not eb:
            n_missing_both += 1
            continue
        if ea and eb and filecmp.cmp(fa, fb, shallow=False):
            same += 1
        else:
            diff.append(os.path.basename(f.format(t="")) + ("" if ea and eb else "(missing)"))
    return same, diff


def compare_preset(a, b):
    same, diff = 0, []
    for f in sorted(os.listdir(b)):
        fa, fb = os.path.join(a, f), os.path.join(b, f)
        if not os.path.exists(fa):
            diff.append(f + "(missing)")
            continue
        if f.endswith(".gz"):
            ok = gzip.open(fa).read() == gzip.open(fb).read()
        elif f.endswith(".json"):
            ja, jb = json.load(open(fa)), json.load(open(fb))
            ok = ja == jb
            if not ok:
                keys = sorted(k for k in set(ja) | set(jb) if ja.get(k) != jb.get(k))
                f = f + "[" + ",".join(keys) + "]"
        else:
            ok = filecmp.cmp(fa, fb, shallow=False)
        if ok:
            same += 1
        else:
            diff.append(f)
    return same, diff


def times(d):
    if not os.path.exists(os.path.join(d, "runtime_sec.txt")):
        return None
    out = {"total": float(open(os.path.join(d, "runtime_sec.txt")).read())}
    st = os.path.join(d, "step_sec.tsv")
    if os.path.exists(st):
        for line in open(st):
            k, v = line.split()
            out[k] = float(v)
    return out


def main(prefix):
    par, tim = [], []
    for done in sorted(glob.glob(f"{W}/bench/*/*/*/DONE")):
        d = os.path.dirname(done)
        org, inp, run = d.split(os.sep)[-3:]
        label, cpu = run.rsplit(".cpu", 1)
        ref = (f"{P}/work/bench_v3/runs/{org}/{inp}/pytd_v7h_homology" if label == "hassist"
               else f"{REL}/bench/{org}/{inp}/{label}")
        par.append(("bench", f"{org}/{inp}", label, cpu, "release", *compare(d, ref, "input.fa")))
        other = os.path.join(os.path.dirname(d), f"{label}.cpu{'8' if cpu == '1' else '1'}")
        if cpu == "1" and os.path.exists(os.path.join(other, "DONE")):
            par.append(("bench", f"{org}/{inp}", label, cpu, "cpu8", *compare(d, other, "input.fa")))
        if label == "rel":
            tim.append((f"{org}/{inp}", int(cpu), times(d),
                        times(f"{V1}/{org}/{inp}/cpu1"), times(f"{V1}/{org}/{inp}/cpu8")))
    for done in sorted(glob.glob(f"{W}/gm/*/*/DONE")):
        d = os.path.dirname(done)
        asm, run = d.split(os.sep)[-2:]
        cpu = run.split(".cpu")[1]
        par.append(("gm", asm, "rel", cpu, "release", *compare(d, f"{REL}/gm/{asm}/rel", "transcripts.fa")))
        if cpu == "1" and os.path.exists(f"{W}/gm/{asm}/rel.cpu8/DONE"):
            par.append(("gm", asm, "rel", cpu, "cpu8", *compare(d, f"{W}/gm/{asm}/rel.cpu8", "transcripts.fa")))
        tim.append((f"gm/{asm}", int(cpu), times(d), times(f"{V1}/gm/{asm}/cpu1"), times(f"{V1}/gm/{asm}/cpu8")))
    for done in sorted(glob.glob(f"{W}/presets/refruns/*/DONE")):
        d = os.path.dirname(done)
        name = os.path.basename(d)
        par.append(("presetbuild", name, "ref_run", "8", "release",
                    *compare(d, f"{REL}/presets/refruns/{name}", "ref.fa")))
        par.append(("presetbuild", name, "preset", "-", "release",
                    *compare_preset(f"{d}/preset", f"{REL}/presets/presets/{name}")))
    for done in sorted(glob.glob(f"{W}/presets/runs/*/*/*/DONE")):
        d = os.path.dirname(done)
        org, set_, run = d.split(os.sep)[-3:]
        setting, cpu = run.rsplit(".cpu", 1)
        par.append(("presettest", f"{org}/{set_}", setting, cpu, "release",
                    *compare(d, f"{REL}/presets/runs/{org}/{set_}/{setting}", "input.fa")))
    with open(prefix + ".parity.tsv", "w") as fh:
        fh.write("kind\tinput\tsetting\tcpu\tcompared_to\tfiles_identical\tfiles_differing\n")
        for r in par:
            fh.write("\t".join(map(str, r[:6])) + "\t" + (",".join(r[6]) or "-") + "\n")
    with open(prefix + ".timing.tsv", "w") as fh:
        fh.write("input\tcpu\ttotal_s\tlongorfs_s\tpredict_s\trelease_cpu1_total_s\tparallel_cpu8_total_s\t"
                 "speedup_vs_release_cpu1\tspeedup_vs_parallel_same_cpu\n")
        for inp, cpu, t, v1, v8 in tim:
            if t is None:
                continue
            same_cpu = v1 if cpu == 1 else v8
            f = lambda x: f"{x:.1f}" if x is not None else "NA"
            fh.write(f"{inp}\t{cpu}\t{t['total']:.1f}\t{f(t.get('longorfs'))}\t{f(t.get('predict'))}\t"
                     f"{f(v1 and v1['total'])}\t{f(v8 and v8['total'])}\t"
                     f"{(v1['total'] / t['total']) if v1 else float('nan'):.2f}\t"
                     f"{(same_cpu['total'] / t['total']) if same_cpu else float('nan'):.2f}\n")
    n_bad = sum(1 for r in par if r[6])
    print(f"{len(par)} comparisons, {n_bad} with differences")


if __name__ == "__main__":
    main(sys.argv[1])
