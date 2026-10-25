# E23. Version history: key TransDecoder releases on the current benchmark (supplement; 2026-10-24)

## Claims supported
- TransDecoder dates from 2011–2012. It began as a Trinity utility (`transcripts_to_best_scoring_ORFs.pl`), was released on SourceForge as transdecoder-r2012-08-07, and moved to GitHub with v2.0 (tag date 2015-01-26).
- On the current benchmark, every key release keeps calls on random sequences low (median ≤ 0.5%). Accuracy improved in steps that line up with method changes:
  - **exact start sites, complete + random:**
    - 0.666 (r2012) → 0.744 (v4.1: variable-order Markov model, PWM start refinement) → 0.758 (v5.7.1);
    - → 0.897 (v7: start refinement by frame contrast plus Kozak context).
  - **fragments, stop sensitivity:**
    - 0.76–0.77 from r2012 through v5.7.1;
    - → 0.815 (v7 length rule: 50 aa minimum with a frame-contrast test below 100 codons);
    - → 0.872 (v7 decoy FDR rule, q 0.05; E20, E22).
  - **calls on lncRNAs:**
    - 13–19% in every release up to and including the v7 length rule;
    - → 5.0 / 5.8% with the decoy FDR rule (E20). This is the first change to substantially reduce them.
- The Perl devel snapshot in the main benchmark (`perl_td_dev_v2`) and PyTD with the v7 length rule give the same median scores (to 3 decimals) on all four input types. They are not identical per input: 6 of 20 inputs match exactly, and the rest differ by at most 6 counted predictions (TP stop, TP exact, predictions, negatives with a call).
- Per the key-version analysis (`VERSIONS.md` §3), v6.0.0 gives predictions identical to v5.7.1, and is represented by it.

## Feature timeline (from `VERSIONS.md`; a code reading of each release)
| release | date | method changes that matter for these benchmarks |
|---|---|---|
| Trinity utility | 2011–2012 | ORF predictor inside Trinity (Haas et al. 2013) |
| r2012-08-07 (SourceForge) | 2012-08 | 5th-order framed Markov model trained on the 500 longest ORFs; keep an ORF if frame-1 score > 0 and > the other 5 frames, or ≥ 900 nt, or Pfam hit; partial ORFs; `-m 100` |
| v2.0 (GitHub; r2013–v3.0.1 line) | 2013–2017 | two steps (LongOrfs / Predict); training set made non-redundant with cd-hit-est; BLAST/Pfam retention |
| v4.1.0 | 2017-08 | variable-order Markov model with pseudocounts, terminal stop excluded; in-house redundancy filter; single best ORF per transcript by default; PWM start refinement on by default |
| v5.0.0 (v5.0–v5.3) | 2017-08 | same-strand 3-frame rule; GC-dependent long-ORF retention |
| v5.4–v6.0.0 (v5.7.1 ≡ v6.0.0) | 2018–2026 | multiple non-overlapping ORFs per transcript; homology-supported ORFs exempt from the overlap filter; NCBI genetic-code names (v5.1) |
| v7 (this paper; Perl + PyTD) | 2026 | start refinement by frame contrast + Kozak context; 50 aa minimum with short-ORF frame-contrast test; regenerated retention table; code-aware stops; deterministic output; then target-decoy FDR acceptance with shuffled transcripts (default q 0.05; `--precise` / `--sensitive`) and presets for small or partial inputs |

## Runs
- 80 fresh runs: 4 releases × 20 benchmark inputs (`data/benchmark/history/jobs.txt`, `scripts/hist_run.sh`). Each release used its own defaults, with `-S` (stranded inputs). `-G` was given only for non-standard codes.
  - Releases up to v5.2 need DB_File, supplied by `work/history/compat`.
  - None failed.
- v5.7.1, Perl devel and the v7 rules come from the main benchmark (`runs/*/perl_td_v571`, `runs/*/perl_td_dev_v2`, `rules/*/cliff_c200`, `rules/*/decoy_q0.05_longgc`).
- **Coordinate normalisation:** r2012 reports 3′-partial ORFs extended 3 nt beyond the transcript end. Its predictions were trimmed back inside the transcript in whole codons before scoring (`clip_preds_to_transcript.py`; `preds.clipped.tsv`). The later releases need no trimming (0 clipped in a check of v2.0 and v4.1.0). Without the trim, r2012 fragment sensitivity scores 0.37 instead of 0.76.
- Runtime, complete + random, single thread, seconds:

  | release | P. falciparum | human | A. thaliana |
  |---|---|---|---|
  | r2012 | 100 | 296 | 299 |
  | v2.0 | 6,005 | 3,849 | 2,064 |
  | v4.1.0 | 225 | 466 | 446 |
  | v5.0.0 | 233 | 518 | 478 |

  v2.0's time is dominated by cd-hit-est on the training ORFs.

## Results: medians over organisms (stop sensitivity / exact sensitivity / % negatives with a prediction)
| release | complete.wRand | fragments.wRand | complete.wLnc | fragments.wLnc |
|---|---|---|---|---|
| r2012 | 0.935 / 0.666 / 0.1 | 0.759 / 0.620 / 0.0 | 0.965 / 0.618 / 16.2 | 0.781 / 0.633 / 15.7 |
| v2.0 | 0.936 / 0.667 / 0.1 | 0.764 / 0.623 / 0.0 | 0.966 / 0.618 / 16.4 | 0.784 / 0.636 / 15.8 |
| v4.1.0 | 0.936 / 0.744 / 0.1 | 0.757 / 0.621 / 0.0 | 0.956 / 0.747 / 13.4 | 0.772 / 0.643 / 13.3 |
| v5.0.0 | 0.946 / 0.748 / 0.3 | 0.772 / 0.625 / 0.1 | 0.969 / 0.757 / 17.9 | 0.792 / 0.658 / 17.6 |
| v5.7.1 (≡ v6.0.0) | 0.946 / 0.758 / 0.3 | 0.772 / 0.623 / 0.1 | 0.969 / 0.764 / 17.9 | 0.792 / 0.653 / 17.6 |
| v7, length rule (PyTD; Perl devel equal at the median) | 0.957 / 0.897 / 0.5 | 0.815 / 0.708 / 0.1 | 0.976 / 0.905 / 18.9 | 0.858 / 0.753 / 18.6 |
| v7, decoy FDR q 0.05 | 0.965 / 0.899 / 0.4 | 0.872 / 0.773 / 0.2 | 0.964 / 0.894 / 5.0 | 0.857 / 0.755 / 5.8 |

Files: `data/benchmark/history/history_scores.tsv`, `history.medians.tsv`, `history.{pdf,png}` (`plot_history.py`), `VERSIONS.md`.

## Text guidance (Brian)
- Describe what each release changed in method. Do not discuss defects of earlier releases.
- Describe AI assistance in the v7 development in a measured, factual Methods statement, as journal policy requires.
