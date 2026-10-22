#!/usr/bin/env python3
"""Grid of PyTransDecoder start-refinement parameters on the stranded bench_v3 inputs.

usage: run_grid.py [--jobs 6] [--alphas 0,0.1,...] [--weights 0,5,...] [--orgs ...]
                   [--inputs complete,fragments] [--no-baseline] [--dry-run]

For every (org, input, config) a fresh run directory is made under
$P/work/start_param_grid/runs/<org>/<input>/<config>/ holding
  input.fa                                   -> symlink to the bench_v3 input
  input.fa.transdecoder_dir/longest_orfs.{cds,gff3,pep}, base_freqs.dat
                                             -> symlinks to the LongOrfs output of the
                                                final-default run (pytd_dev_v2)
Nothing else is linked: in particular not longest_orfs.cds.scores (a Predict output that
Predict would rewrite through the link) and not __checkpoints_TDpredict (which would make
Predict skip its steps). Predict then runs in full, single-threaded, with
  -G <code> --start-contrast-alpha A --start-kozak-weight W     (config a<A>_w<W>)
  -G <code> --no-refine-starts                                  (config norefine)
and its GFF3 is converted (convert_predictions.py transdecoder_gff3) and scored
(score_predictions.py) against <input>.truth.tsv. Writes preds.tsv, score.tsv, runtime_sec.txt.
Runs that already have score.tsv are skipped, so the script can be restarted.
"""
import argparse
import concurrent.futures as cf
import os
import subprocess
import sys
import time

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
T = f"{P}/TransDecoderPaper"
PY = f"{P}/env/td/bin/python"
PYTD = f"{P}/work/trees/PY_dev"
INPUTS = f"{P}/data/bench_v3/inputs"
SRC_RUNS = f"{P}/work/bench_v3/runs"
GRID = f"{P}/work/start_param_grid/runs"
SCRIPTS = f"{T}/benchmark_v2/scripts"
ORGS = "athal creinhardtii dmel human mouse pfalciparum spombe tthermophila".split()
LONGORFS_FILES = ["longest_orfs.cds", "longest_orfs.gff3", "longest_orfs.pep", "base_freqs.dat"]


def fmt(x):
    return ("%g" % x)


def config_name(alpha, weight):
    return "norefine" if alpha is None else f"a{fmt(alpha)}_w{fmt(weight)}"


def genetic_codes():
    codes = {}
    with open(f"{T}/benchmark_v2/organisms.tsv") as fh:
        next(fh)
        for line in fh:
            x = line.rstrip("\n").split("\t")
            codes[x[0]] = x[3]
    return codes


def run_one(org, inp, alpha, weight, gcode):
    cfg = config_name(alpha, weight)
    rdir = f"{GRID}/{org}/{inp}/{cfg}"
    if os.path.exists(f"{rdir}/score.tsv"):
        return cfg, org, inp, "skip", 0.0
    src = f"{SRC_RUNS}/{org}/{inp}/pytd_dev_v2/input.fa.transdecoder_dir"
    wdir = f"{rdir}/input.fa.transdecoder_dir"
    os.makedirs(wdir, exist_ok=True)
    if os.path.exists(f"{wdir}/__checkpoints_TDpredict"):
        # an interrupted earlier attempt: move it aside (no rm) so Predict starts clean
        os.rename(f"{wdir}/__checkpoints_TDpredict", f"{wdir}/__checkpoints_TDpredict.stale.{int(time.time())}")
    if not os.path.lexists(f"{rdir}/input.fa"):
        os.symlink(f"{INPUTS}/{org}/{inp}.fa", f"{rdir}/input.fa")
    for fn in LONGORFS_FILES:
        if not os.path.lexists(f"{wdir}/{fn}"):
            os.symlink(f"{src}/{fn}", f"{wdir}/{fn}")
    args = [PY, "-m", "pytransdecoder.cli", "predict", "-t", "input.fa", "-G", gcode]
    args += ["--no-refine-starts"] if alpha is None else \
        ["--start-contrast-alpha", fmt(alpha), "--start-kozak-weight", fmt(weight)]
    env = dict(os.environ, PYTHONPATH=PYTD, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    t0 = time.time()
    with open(f"{rdir}/predict.log", "w") as log:
        rc = subprocess.run(args, cwd=rdir, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t0
    if rc != 0:
        return cfg, org, inp, f"predict failed rc={rc}", dt
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


def floats(s):
    return [float(x) for x in s.split(",") if x != ""]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--alphas", default="0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0")
    ap.add_argument("--weights", default="0,5,10,25,50,100")
    ap.add_argument("--orgs", default=",".join(ORGS))
    ap.add_argument("--inputs", default="complete,fragments")
    ap.add_argument("--no-baseline", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    codes = genetic_codes()
    configs = [(al, w) for al in floats(a.alphas) for w in floats(a.weights)]
    if not a.no_baseline:
        configs.append((None, None))
    # biggest inputs first, default config first within an input (sanity check early)
    size = {(o, i): os.path.getsize(f"{SRC_RUNS}/{o}/{i}/pytd_dev_v2/input.fa.transdecoder_dir/longest_orfs.cds")
            for o in a.orgs.split(",") for i in a.inputs.split(",")}
    jobs = []
    for (o, i) in sorted(size, key=size.get, reverse=True):
        for c in sorted(configs, key=lambda c: c != (0.4, 25.0)):
            jobs.append((o, i) + c)
    print(f"{len(jobs)} runs, {a.jobs} in parallel", flush=True)
    if a.dry_run:
        for j in jobs:
            print(j)
        return
    done = 0
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(run_one, o, i, al, w, codes[o]) for o, i, al, w in jobs]
        for f in cf.as_completed(futs):
            done += 1
            cfg, o, i, status, dt = f.result()
            print(f"{time.strftime('%H:%M:%S')} [{done}/{len(jobs)}] {o}/{i}/{cfg} {status} {dt:.0f}s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
