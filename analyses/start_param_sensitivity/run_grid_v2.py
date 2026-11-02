#!/usr/bin/env python3
"""alpha x w grid at the pinned release code (E40), on the stranded complete and fragment inputs of
the 8 benchmark organisms. Successor of run_grid.py (kept), which ran development code with w in nats.

usage: run_grid_v2.py [--jobs 4] [--dry-run]

w is in bases per bit, on the released scale. Old grid (per nat) -> this grid (per bit):
    0 -> 0, 5 -> 3.5, 10 -> 7, 25 -> 17, 50 -> 35, 100 -> 69
The default cell is exactly (0.4, 17). alpha: 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.25, 1.5, 2, 3
(78 cells), plus --no-refine-starts (79 configs x 16 inputs = 1,264 full Predict runs).

Each run: a new directory $P/work/release_pin/start_grid/runs/<org>/<input>/<config>/; LongOrfs output
of the pinned-commit release run (work/release_pin/bench/<org>/<input>/rel) symlinked: longest_orfs.{cds,gff3,pep},
base_freqs.dat, longorfs.params.json. Predict then runs in full (training, scoring, decoys, selection,
refinement). A run without score.tsv is moved to $P/_trash first (no resume).
"""
import argparse
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
import time

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
T = f"{P}/TransDecoderPaper"
PY = f"{P}/env/td/bin/python"
PYTD = f"{P}/work/trees/PY_rel_d3892b5"
INPUTS = f"{P}/data/bench_v3/inputs"
SRC = f"{P}/work/release_pin/bench"
GRID = f"{P}/work/release_pin/start_grid/runs"
SCRIPTS = f"{T}/benchmark_v2/scripts"
ORGS = "athal creinhardtii dmel human mouse pfalciparum spombe tthermophila".split()
FILES = ["longest_orfs.cds", "longest_orfs.gff3", "longest_orfs.pep", "base_freqs.dat", "longorfs.params.json"]
ALPHAS = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2, 3]
WEIGHTS = [0, 3.5, 7, 17, 35, 69]


def cfgname(a, w):
    return "norefine" if a is None else f"a{a:g}_w{w:g}"


def codes():
    out = {}
    with open(f"{T}/benchmark_v2/organisms.tsv") as fh:
        next(fh)
        for line in fh:
            x = line.rstrip("\n").split("\t")
            out[x[0]] = x[3]
    return out


def run_one(org, inp, a, w, gcode):
    cfg = cfgname(a, w)
    rdir = f"{GRID}/{org}/{inp}/{cfg}"
    if os.path.exists(f"{rdir}/score.tsv"):
        return cfg, org, inp, "skip", 0.0
    if os.path.exists(rdir):
        dst = f"{P}/_trash/release_pin_partial/grid.{org}.{inp}.{cfg}.{int(time.time())}"
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.move(rdir, dst)
    wdir = f"{rdir}/input.fa.transdecoder_dir"
    os.makedirs(wdir)
    os.symlink(f"{INPUTS}/{org}/{inp}.fa", f"{rdir}/input.fa")
    src = f"{SRC}/{org}/{inp}/rel/input.fa.transdecoder_dir"
    for fn in FILES:
        os.symlink(f"{src}/{fn}", f"{wdir}/{fn}")
    args = [PY, "-m", "pytransdecoder.cli", "predict", "-t", "input.fa", "-G", gcode]
    args += ["--no-refine-starts"] if a is None else \
        ["--start-contrast-alpha", f"{a:g}", "--start-kozak-weight", f"{w:g}"]
    env = dict(os.environ, PYTHONPATH=PYTD, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    t0 = time.time()
    with open(f"{rdir}/predict.log", "w") as log:
        rc = subprocess.run(args, cwd=rdir, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t0
    if rc != 0:
        return cfg, org, inp, f"FAILED rc={rc}", dt
    with open(f"{rdir}/preds.tsv", "w") as out:
        subprocess.run([PY, f"{SCRIPTS}/convert_predictions.py", "transdecoder_gff3",
                        "input.fa.transdecoder.gff3"], cwd=rdir, stdout=out, check=True)
    sc = subprocess.run([PY, f"{SCRIPTS}/score_predictions.py", "--truth", f"{INPUTS}/{org}/{inp}.truth.tsv",
                         "--preds", f"{rdir}/preds.tsv", "--label", cfg],
                        capture_output=True, text=True, check=True).stdout
    with open(f"{rdir}/runtime_sec.txt", "w") as fh:
        fh.write(f"{dt:.1f}\n")
    with open(f"{rdir}/score.tsv.tmp", "w") as fh:
        fh.write(sc)
    os.rename(f"{rdir}/score.tsv.tmp", f"{rdir}/score.tsv")
    return cfg, org, inp, "ok", dt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    cs = codes()
    configs = [(0.4, 17)] + [(al, w) for al in ALPHAS for w in WEIGHTS if (al, w) != (0.4, 17)] + [(None, None)]
    jobs = [(o, i) + c for c in configs for o in ORGS for i in ("complete", "fragments")]
    print(f"{len(jobs)} runs, {a.jobs} in parallel", flush=True)
    if a.dry_run:
        return
    n = 0
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(run_one, o, i, al, w, cs[o]) for o, i, al, w in jobs]
        for f in cf.as_completed(futs):
            n += 1
            cfg, o, i, st, dt = f.result()
            print(f"{time.strftime('%m-%d %H:%M:%S')} [{n}/{len(jobs)}] {o}/{i}/{cfg} {st} {dt:.0f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
