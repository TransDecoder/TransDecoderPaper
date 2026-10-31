# Reproducing the TD2 comparison

For the TD2 authors and anyone else checking the comparison. TD2 was run ab initio with its defaults,
without its homology options. Results and caveats: `evidence/E39_td2_comparator.md`.

## Software
- `td2` 1.1.0 and `psauron` 1.1.3 from PyPI: `pip install td2==1.1.0`. Full environment: `td2_env.pip-freeze.txt`.
- Scorer: `scripts/score_predictions.py` (the same scorer as for every other method).

## Run one input
`scripts/run_td2.sh <input.fa> <outdir> <org>` does the following:

```bash
TD2.LongOrfs -t input.fa [-S] -G <NCBI table> -@ 2     # -S unless the input name contains "unstranded"
TD2.Predict  -t input.fa -G <NCBI table>               # table 6 for T. thermophila, otherwise 1
```

PSAURON is limited to two threads with `OMP_NUM_THREADS=2 MKL_NUM_THREADS=2`. `<org>` looks up the
table in `organisms.tsv`. The script converts `input.fa.TD2.gff3` to `preds.tsv` with
`scripts/convert_predictions.py transdecoder_gff3`, and stops if a transcript has more than one
prediction (none did). `BENCH_ROOT`, `PYTHON` and `TD2_ENV` override its default paths. The launcher and table scripts
below have the original project directory (`P=...`) written in, so edit it there.

## Score one input
```bash
python scripts/score_predictions.py --truth <input>.truth.tsv --preds preds.tsv | tail -1
```

A stop-level match needs the same strand and a 3' end within 2 nt of the annotated one. An exact
match also needs a 5' end within 2 nt. Scoring is described in `README.md`.

## All 120 runs
- Job list: `launch_v3/td2_jobs_v1.tsv` (input, output directory, organism): 84 benchmark inputs, 24
  mixed inputs, and 12 simulated assemblies (human and *A. thaliana*, d1-d32).
- Launcher: `launch_v3/launch_td2_v1.sh`, 4 jobs at a time.
- Inputs: built with `scripts/make_all_inputs.sh` (seed 1) from the annotations in `truth_sources.tsv`.
  The input FASTA files are not in this repository.
- Scoring and tables: `pipeline/derive/10_td2_tables.sh`, then `pipeline/derive/10b_td2_mechanisms.py`.

## Expected results
Per-input scores are in `data/supplementary_tables/benchmark_per_input.tsv` and
`data/benchmark/final_v7/all_scores.tsv` (tool `td2`).

We reran two *S. pombe* inputs from scratch. The stranded input (`fragments`) gave byte-identical
`preds.tsv`. The unstranded input (`fragments.unstranded`) differed in 3 of 3,957 predictions: the same
coordinates on the opposite strand, moving stop F1 from 0.8647 to 0.8645. Candidates and PSAURON
scores were identical. Two reruns with `PYTHONHASHSEED=1` matched each other exactly, so these
differences follow Python's hash ordering, and small differences are expected on unstranded inputs
unless `PYTHONHASHSEED` is set. It was not set for the published runs.
