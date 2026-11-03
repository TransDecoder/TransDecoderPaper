# E42. Parallel and faster PyTD scoring with byte-identical outputs (2026-11-02)

Engineering record, not a paper claim yet. Asked by Brian Haas (2026-11-02): parallelise candidate scoring and
the decoy step, then find other speedups, keeping outputs identical and deterministic.

## Release
PyTD **v0.3.0** (2026-11-02, BH approved): `devel` branch created from dev/faster-serial plus a release commit
(version 0.3.0, NumPy declared as a dependency, README `--cpu`, changelog). `main` fast-forwarded d64906d -> b3c5f15,
tag v0.3.0, GitHub release (Latest), PyPI `pytransdecoder` 0.3.0. Pre-release checks:
- wheel/sdist built from a clean `git archive` of devel; `twine check` passed;
- 84 of 84 tests pass from the installed wheel in fresh venvs on Python 3.9 (NumPy 2.0.2, Biopython 1.85),
  3.11 (2.4.6, 1.88) and 3.12 (2.5.4, 1.88);
- end to end, PyPI v0.2.1 (`--cpu 1`) vs the v0.3.0 wheel (`--cpu 1` and `--cpu 4`) on Python 3.9 and 3.12:
  human fragments.wRand, P. falciparum fragments, T. thermophila fragments (ciliate code), human fragments with
  homology start assist; 16 of 16 comparisons identical (17-18 files each);
- `spawn` start method forced, `--cpu 3`: identical to `--cpu 1` (S. pombe fragments);
- after upload, `pip install pytransdecoder==0.3.0` in a clean venv, `pyTransdecoder --cpu 4`: identical.

## Code (PyTD branches, merged into devel/main for v0.3.0)
- `dev/parallel-scoring` 98d77fb: `--cpu N` (predict and full pipeline). New `core/parallel.py`; six-frame
  candidate scoring and decoy shuffling/scoring run in N fork worker processes. Work is chunked in input order
  and reassembled in that order. Each decoy shuffle is seeded per transcript (`seed:tid:k`), so results do not
  depend on N.
- `dev/faster-serial` 5859887 (on top of the branch above, plus tests 645ad38):
  - NumPy Markov scorer (`core/markov_fast.py`). It adds the values with `np.add.accumulate` over
    `[0.0, v0, v1, ...]`, which does the reference loop's IEEE additions in the same order. Score tables with any
    other key fall back to the loop.
  - Per-codon memoised translation filled from Biopython's own codon results.
  - ORF finder: a regex scan of all frames plus `bisect` for the first stop. Both fall back to the original scans.
  - Exact rewrites of Counter k-mer and redundancy tallies (encounter order kept for ties), base counts, the
    reverse-complement table (built from Biopython) and the inlined xoshiro128** shuffle.
  - Transcripts and final models are loaded once per run.
- Tests: 84 pass, including the Perl-shared PRNG/decoy golden values and randomised bit-identity tests. A
  differential check of the new ORF finder against v0.2.1 found 122,479 ORFs identical.

## Method
- Hosts: methods01, 28 cores, shared load about 20-48 during the runs. Runs were sequential, at normal priority.
  Times are indicative.
- v1: `pipeline/parallel_scoring/{run_one_v1.sh,launch_v1.sh,summarise_v1.py}`, branch dev/parallel-scoring,
  `--cpu` 1/2/4/8 on nine benchmark inputs and GM12878 Trinity. Output: `data/parallel_scoring/v1.{parity,timing}.tsv`.
- v2: `pipeline/parallel_scoring/{run_one_v2.sh,launch_v2.sh,summarise_v2.py}`, branch dev/faster-serial, 156 runs:
  - all 108 benchmark inputs at release defaults, `--cpu 8`, and the nine timing inputs again at `--cpu 1`;
  - homology retention-only and start assist, fragments and complete inputs of each organism;
  - GM12878 on all three assemblies, plus Trinity at `--cpu 1`;
  - a human preset build and two preset-mode runs.
  Output: `data/parallel_scoring/v2.{parity,timing}.tsv`.
- Parity: 15 files per run, compared byte for byte (final gff3/bed/pep/cds and the `.transdecoder_dir` files,
  including decoy q-values and the trained models). References are the pinned release runs (`work/release_pin`,
  PyTD d3892b5) and, for the start assist, `work/bench_v3/runs/*/pytd_v7h_homology` (a5b24c3, same code tree).
  Preset files were compared after gunzip, and preset.json as parsed JSON.

## Results
- Parity: v1, 39 of 39 comparisons identical. v2, 167 of 167 identical:
  - 149 benchmark runs against release and 9 cpu1-vs-cpu8 pairs;
  - 4 GM12878 runs against release and 1 cpu1-vs-cpu8 pair;
  - preset build (reference run and 5 preset files);
  - 2 preset-mode runs.
- Wall time, total of LongOrfs and Predict, on the ten timing inputs:
  - dev/faster-serial at `--cpu 1` against release: 1.75-3.45x faster (human complete 277 -> 101 s;
    Trinity 1806 -> 767 s).
  - dev/faster-serial at `--cpu 8` against release: 5.3-10.5x faster (human complete 277 -> 40 s; Trinity
    1806 -> 222 s; T. thermophila unstranded wRand 1737 -> 166 s).
  - dev/faster-serial at `--cpu 8` against dev/parallel-scoring at `--cpu 8`: 1.6-2.5x faster.
- Remaining time at `--cpu 8` is mostly serial:
  - LongOrfs: 7.6 s on human complete, 62 s on Trinity. It has no `--cpu` yet.
  - Redundancy filter, k-mer training and output writing.
  - Decoy PRNG draws are sequential within a shuffle but split across workers.
- Not changed in behaviour: Biopython's partial-codon warning is no longer logged by the memoised path. Logs
  only, no output change.

## LongOrfs --cpu (2026-11-03; released as PyTD v0.4.0)
Asked by BH ("look into adding --cpu to LongOrfs"). In v0.3.0 at `--cpu 8`, LongOrfs was the largest remaining
single-process step (62 s of 222 s on GM12878 Trinity).
- Change:
  - transcripts are split into consecutive chunks; workers extract ORFs and format the GFF3/CDS/PEP text, and the
    chunks are written in input order;
  - the input is read once (it was read three times) and base counts are summed per chunk (integers, so exact);
  - model IDs `<acc>.pN` of different accessions never coincide, because the last `.p` of an ID starts its numeric
    suffix, so chunks number them independently. Inputs with repeated accessions fall back to one process;
  - `--cpu 1` streams in batches with one model-ID set per run;
  - `pyTransdecoder --cpu` also applies to LongOrfs. With N > 1 the transcripts are held in memory (Trinity: 570 MB
    peak against 72 MB).
- Tests: 91 pass, 7 new:
  - 1 worker against 2, 3 and 8; both strand settings; universal and ciliate codes;
  - accessions shaped like model IDs; repeated accessions, also across small serial batches; gene-trans map;
    complete-only.
  - The repeated-accession output is identical to released v0.3.0.
- Series: `pipeline/parallel_scoring/{run_lo_v3.sh,launch_lo_v3.sh,summarise_lo_v3.py}` ->
  `data/parallel_scoring/lo_v3.{parity,timing}.tsv`. 226 of 226 comparisons byte-identical with the release runs:
  - LongOrfs at `--cpu 1` and `--cpu 8` on all 108 benchmark inputs and GM12878 LRAA/StringTie/Trinity
    (5 files each);
  - `pyTransdecoder --cpu 8` on four inputs (15 LongOrfs + Predict files each).
- Timing, LongOrfs at `--cpu 8` against v0.3.0 (single process), methods01, load about 21-25:
  - median 3.4x faster (range 0.9-5.4x);
  - Trinity 62 -> 14.5 s; T. thermophila unstranded wRand 46 -> 12 s;
  - the 0.9x case is a 2.3 s input, where starting the pool costs as much as it saves.
- End to end, `pyTransdecoder --cpu 8`:
  - human complete 34 s (v0.3.0 40 s);
  - T. thermophila complete.unstranded 78 s (v0.3.0 95 s).

### Release v0.4.0 (2026-11-03, BH: "ship this as 0.4.0")
- main fast-forwarded b3c5f15 -> a285464 (= devel); tag v0.4.0; GitHub release (Latest); PyPI `pytransdecoder` 0.4.0.
- Checks:
  - wheel/sdist built from `git archive`; `twine check` passed;
  - 91 of 91 tests pass from the installed wheel on Python 3.9, 3.11 and 3.12;
  - end to end, `pyTransdecoder` v0.3.0 against the 0.4.0 wheel at `--cpu 4`: 18 of 18 files identical on human
    fragments.wRand and on T. thermophila complete.unstranded;
  - after upload, `pip install pytransdecoder==0.4.0` in a clean venv: output identical to v0.3.0.
- Python 3.12 shows a "multi-threaded ... fork()" DeprecationWarning, but only inside the pytest runner (also in
  v0.3.0). There are no warnings in CLI runs.

### Disk cleanup (2026-11-03, BH approved)
- The output files of the verification runs in `work/{parallel_scoring,fast_scoring,longorfs_cpu}` were deleted.
  Kept: runtime_sec.txt, step_sec.tsv, pytd_version.txt, logs and DONE. The parity and timing results are the
  committed `data/parallel_scoring/*.tsv`. The summarise_* scripts can no longer recompute parity without
  rerunning.
- Also deleted:
  - the Perl `*.transdecoder_dir/decoys/` subdirectories (123: perl_td_v7 benchmark runs and release_pin
    regression). No script reads them; Perl regenerates them deterministically;
  - the contents of `$P/_trash`;
  - this session's scratchpad;
  - worktrees PY_par, PY_fast and PY_lo, and their merged branches.
