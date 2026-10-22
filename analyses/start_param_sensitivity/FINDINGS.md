# Start-refinement parameter sensitivity (alpha, w)

**Question.** TransDecoder trims an ORF with an open 5' end to its first in-frame ATG when
C(u) < (alpha*|u| + w*K)*mu (predict.py _revise_start_sites). The defaults are alpha = 0.4
(--start-contrast-alpha) and w = 25 (--start-kozak-weight). They were set during development
without a saved rationale (E09 Gaps). Are they well chosen, how sensitive is accuracy to them,
and are they overfit?

**Answer.** Keep alpha = 0.4, w = 25. When complete and fragment inputs count equally, the default
is the best of 78 grid cells. Leave-one-organism-out picks the default itself in 7 of 8 folds. On
GM12878 the default sits at the balance point between more exact starts and fewer wrongly trimmed
truncated ORFs. alpha carries a real trade-off: higher alpha favours full-length transcripts and
lower alpha favours fragmentary ones. The documentation should say so.

Full record: `evidence/E13_start_param_sensitivity.md`.

## What was run
- PyTD PY_dev 48c870c, unmodified. predict ran in full, single-threaded, on symlinked
  longest_orfs.{cds,gff3,pep} and base_freqs.dat from the pytd_dev_v2 runs. No checkpoints and no
  .cds.scores were copied.
- bench_v3 stranded complete and fragments, 8 organisms, -G from organisms.tsv, scored with
  score_predictions.py.
- Grid: alpha {0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0} x w {0,5,10,25,50,100}, plus
  --no-refine-starts (880 runs).
- Extension: alpha {1.25,1.5,2,3} x the same w (384 runs), added because the complete-input
  optimum sat at the alpha = 1.0 edge. Total: 1,264 runs, 41 CPU-h, 6.9 h wall at 6 jobs.
- GM12878: StringTie and LRAA final_abinitio LongOrfs output; alpha {0.2,0.4,0.6,1.0} x
  w {0,25,50}, plus no-refine (26 runs); predict -> genome-map -> eval_genome_cds.py
  (MANE v1.5, GENCODE v50).

## Checks
- (0.4, 25) preds.tsv is byte-identical to pytd_dev_v2 on all 16 inputs.
- GM12878 StringTie (0.4, 25) gff3 is byte-identical to final_abinitio's. Both assemblies' eval
  summaries match eval/<asm>_final_abinitio.summary.tsv.
- n_pred and tp_stop are the same across all 79 configs of each input, including no-refine
  (checks.tsv). Exact sensitivity and precision therefore rank cells as exact F1 does.
- The size and mtime of the source LongOrfs files and preds.tsv were unchanged after the runs.

## Results
### 1. Complete and fragment inputs want opposite alpha
Exact F1 x100 (optima.tsv):

| org | complete default | complete best (cell) | gap | fragments default | fragments best (cell) | gap | fragments no-refine |
|---|---|---|---|---|---|---|---|
| athal | 89.7 | 96.3 (3,0) | 6.6 | 72.3 | 75.3 (0,25) | 3.0 | 72.4 |
| creinhardtii | 74.2 | 78.0 (3,0) | 3.9 | 76.1 | 77.0 (0.1,10) | 1.0 | 74.6 |
| dmel | 92.2 | 94.4 (3,0) | 2.2 | 79.4 | 81.3 (0,10) | 1.8 | 78.5 |
| human | 92.2 | 96.8 (3,0) | 4.6 | 81.2 | 82.2 (0.1,25) | 1.0 | 72.6 |
| mouse | 90.8 | 93.9 (3,0) | 3.1 | 80.2 | 81.6 (0.1,25) | 1.4 | 78.1 |
| pfalciparum | 92.8 | 93.6 (3,0) | 0.7 | 80.6 | 83.7 (0,5) | 3.1 | 84.8 |
| spombe | 91.2 | 97.2 (3,0) | 6.0 | 78.5 | 81.0 (0,25) | 2.5 | 77.6 |
| tthermophila | 84.0 | 90.6 (3,0) | 6.6 | 76.2 | 76.7 (0.2,25) | 0.5 | 64.8 |

- Complete inputs contain no truly truncated ORFs: among no-refine open ORFs that match a truth
  stop, 41,127 extend into the UTR and 0 are truncated (open_end_mix.tsv). "Always trim" is
  therefore optimal, F1 rises monotonically with alpha, and that optimum is degenerate.
- About 75% of open ORFs in fragments are truncated (53,571 vs 17,922). The optimum is
  alpha 0-0.2 with w ~25.
- Within one input type, at most 6 of 78 cells are within 0.5 points of an organism's best
  (pfalciparum complete: 37).
- Equal weighting gives a moderate plateau. Eight cells are within 0.5 points of the best mean:
  (0.2,25), (0.3,10/25/50), (0.4,10/25/50) and (0.5,25).
- pfalciparum fragments is the one input where refinement costs accuracy (no-refine 84.8,
  default 80.6, best grid cell 83.7). athal fragments are flat (72.4 vs 72.3).

### 2. Weighting complete (lambda) against fragments (mix.tsv)
| lambda | best cell | default gap (pts) | default rank /78 |
|---|---|---|---|
| 0.0 | (0,25) | 1.46 | 16 |
| 0.3 | (0.2,25) | 0.30 | 4 |
| 0.4 | (0.3,25) | 0.10 | 3 |
| 0.5 | (0.4,25) | 0 | 1 |
| 0.6 | (0.4,25) | 0 | 1 |
| 0.7 | (0.5,25) | 0.13 | 4 |
| 0.8 | (0.6,10) | 0.52 | 12 |
| 0.9 | (1,10) | 1.45 | 46 |
| 1.0 | (3,0) | 4.22 | 47 |

### 3. Leave-one-organism-out (loo.tsv)
- Equal-weight objective: 7 of 8 folds choose (0.4,25). The tthermophila fold chooses (0.3,25),
  which scores 0.64 points below the default on tthermophila. The oracle beats the chosen cell by
  0-1.05 points.
- Fragments only: folds choose alpha 0-0.1 with w 25. That gains 0.8-2.6 points on 7 held-out
  organisms and loses 1.1 on tthermophila.
- Complete only: every fold picks the degenerate (3,0).

### 4. GM12878 (used in development; not independent)
"Start precision" = among complete calls compatible with a GENCODE CDS, the fraction with the
exact CDS.

| assembly | cell | MANE exact given assembled | start precision | exact complete calls | wrong-start complete calls | partial calls kept compatible |
|---|---|---|---|---|---|---|
| StringTie | no-refine | 61.03 | 73.27 | 9,309 | 3,396 | 3,486 |
| StringTie | (0.2,25) | 90.69 | 78.58 | 13,504 | 3,680 | 3,603 |
| StringTie | (0.4,25) | 92.70 | 78.43 | 13,805 | 3,796 | 3,580 |
| StringTie | (0.6,25) | 93.79 | 78.07 | 13,980 | 3,928 | 3,531 |
| StringTie | (1.0,25) | 94.92 | 76.98 | 14,178 | 4,239 | 3,330 |
| LRAA | no-refine | 59.06 | 81.67 | 7,562 | 1,697 | 3,274 |
| LRAA | (0.2,25) | 91.80 | 84.86 | 11,768 | 2,100 | 3,173 |
| LRAA | (0.4,25) | 93.45 | 84.34 | 12,018 | 2,231 | 3,112 |
| LRAA | (0.6,25) | 94.47 | 83.50 | 12,199 | 2,410 | 2,968 |
| LRAA | (1.0,25) | 95.47 | 81.46 | 12,381 | 2,817 | 2,625 |

- MANE exact given assembled rises with alpha: it counts only fully assembled MANE CDSs, so
  truncated ORFs cannot hurt it. Start precision falls with alpha.
- From 0.4 to 1.0: about +370 exact complete calls, +443 / +586 wrong-start complete calls, and
  -250 / -487 partial calls correctly left open (StringTie / LRAA).
- From 0.4 to 0.6 roughly breaks even. Dropping to 0.2 costs 1.7-2.0 points of MANE exact.
- GM12878 has about 23-26% truncated ORFs among labelled open ORFs (E09 labels), against about
  48% in the equal-weight benchmark. That corresponds to lambda ~0.75-0.8, where alpha 0.5-0.6
  would be best on the benchmark.

### 5. Kozak weight
- At alpha 0.4, going from w 0 to w 25 raises mean exact F1 from 87.8 to 88.4 (complete) and from
  75.4 to 78.1 (fragments).
- On GM12878 the same change adds 2.8-3.1 points of MANE exact and 0.8-1.1 points of start
  precision.
- w 25 is best or within 0.2 points for alpha 0.2-0.6. w 10-50 are within about 0.5 points; w 100
  is worse.

## Recommendation
1. Keep (0.4, 25) for the release.
2. Document the trade-off in the Methods and in the --start-contrast-alpha help: alpha 0.6-1.0
   for full-length or long-read transcriptomes; alpha 0.1-0.2 for fragmentary assemblies.
3. Optional: alpha 0.5 is defensible. It costs 0.24 points on the benchmark and gains about
   0.5 points of MANE exact on GM12878. The data cannot clearly separate 0.4 from 0.5-0.6.

## Caveats
- The defaults were chosen on this data. Leave-one-out addresses organism-level overfitting, but
  not the 50/50 objective.
- GM12878 is not independent.
- The fragment set is extreme: every transcript is cut inside its CDS.
- Only stranded inputs were run.
- Single deterministic runs: treat differences under about 0.2 points as noise.
- Ribo-seq and MS support were not evaluated.

## Files
- Scripts: run_grid.py, run_gm12878.sh, summarise.py, open_end_mix.py.
- results/: grid_scores, checks, optima, plateau, mix, loo, gm12878 and open_end_mix (.tsv);
  heatmap_{f1,sens}_exact_{complete,fragments}, heatmap_summary and gm12878_heatmap (.png/.pdf).
- Not in git: $P/work/start_param_grid/{runs,gm12878}/ (34 GB), with logs grid.log,
  grid_ext.log, gm12878.log and source_stat_before.txt.
