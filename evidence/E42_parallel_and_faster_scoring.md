# E42. Parallel and faster PyTD scoring with byte-identical outputs (2026-11-02)

Engineering record, not a paper claim yet. Asked by Brian Haas (2026-11-02): parallelise candidate scoring and
the decoy step, then find other speedups, keeping outputs identical and deterministic.

## Code (PyTD, local branches, not pushed or merged; BH to decide)
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
