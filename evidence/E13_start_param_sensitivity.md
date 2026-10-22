# E13. Sensitivity of the start-refinement parameters (alpha, w)

## Claims supported
Manuscript: Methods, "Start-codon refinement" (methods_algorithm.tex, sec:start: how alpha = 0.4
and w = 25 were chosen, and the trade-off). Fills the E09 gap "the grid or summary that selected
alpha and w is not saved".
1. With complete and fragment inputs weighted equally, (0.4, 25) is the best of 78 (alpha, w)
   cells on the 8-organism benchmark (mean exact F1 0.8322). Eight cells (alpha 0.3-0.5,
   w 10-50) lie within 0.5 points.
2. Leave-one-organism-out picks (0.4, 25) in 7 of 8 folds; the T. thermophila fold picks
   (0.3, 25).
3. alpha trades complete-transcript accuracy against fragment accuracy. Complete inputs favour
   always trimming (best at alpha 3); fragments favour alpha 0-0.2. The default is the best cell
   for complete weights of 0.5-0.6 and within 0.13 points for 0.4-0.7.
4. The Kozak term helps. At alpha 0.4, w 25 vs 0 gives +0.6 (complete) and +2.7 (fragments)
   points of mean exact F1, and +2.8-3.1 points of GM12878 MANE exact given assembled.
5. On GM12878, raising alpha above 0.4 gains MANE exact given assembled but lowers the
   exact-start share of complete calls; the default sits between the two.

## Data
- bench_v3 stranded inputs: $P/data/bench_v3/inputs/<org>/{complete,fragments}.{fa,truth.tsv},
  8 organisms (E07).
- LongOrfs output reused from $P/work/bench_v3/runs/<org>/<input>/pytd_dev_v2/input.fa.transdecoder_dir.
- GM12878: $P/work/gm12878_td/{stringtie,lraa}/final_abinitio LongOrfs output; structures
  $P/data/gm12878/star_stringtie/stringtie.gtf and
  $P/data/gm12878/longread_lraa/lraa_denovo/LRAA.gtf; MANE v1.5 and GENCODE v50 (E03).

## Software
- PyTransDecoder PY_dev 48c870c (unmodified); $P/env/td (Python 3.11).
- Scripts: analyses/start_param_sensitivity/{run_grid.py, run_gm12878.sh, summarise.py,
  open_end_mix.py}.

## Commands
    cd $T/analyses/start_param_sensitivity
    python run_grid.py --jobs 6                                         # 55 configs x 16 inputs
    python run_grid.py --jobs 6 --alphas 1.25,1.5,2,3 --no-baseline     # extension
    bash run_gm12878.sh 6                                               # alpha {0.2,0.4,0.6,1.0} x w {0,25,50} + no-refine
    python summarise.py                                                 # -> results/
    python open_end_mix.py > results/open_end_mix.tsv
Runs were on 2026-10-21: grid 14:53-19:42, extension 19:42-21:46, GM12878 21:46-22:13.
Total benchmark predict CPU: 147,507 s.

## Outputs
- results/: grid_scores, checks, optima, plateau, mix, loo, gm12878 and open_end_mix (.tsv);
  heatmap_{f1,sens}_exact_{complete,fragments}, heatmap_summary and gm12878_heatmap (.png/.pdf);
  FINDINGS.md.
- Not in git: $P/work/start_param_grid/runs/<org>/<input>/<config>/ and
  $P/work/start_param_grid/gm12878/<asm>/<config>/ (34 GB), with logs grid.log, grid_ext.log,
  gm12878.log, timing_test.log and source_stat_before.txt.

## Results
See `analyses/start_param_sensitivity/FINDINGS.md` sections 1-5, which give the per-organism
optima table, the complete-weight table, leave-one-out, the GM12878 table and the Kozak weight.
Key numbers:
- Equal-weight mean exact F1: (0.4,25) 0.8322; (0.3,25) 0.8317; (0.5,25) 0.8298;
  (0.4,0) 0.8160.
- Gap between the default and each organism's best: complete 0.7-6.6 points (all optima at
  alpha 3); fragments 0.5-3.1 points (optima at alpha 0-0.2).
- GM12878 at (0.2, 0.4, 0.6, 1.0) x w 25:
  - MANE exact given assembled: StringTie 90.69, 92.70, 93.79, 94.92; LRAA 91.80, 93.45, 94.47,
    95.47.
  - Exact-start share of GENCODE-compatible complete calls: StringTie 78.58, 78.43, 78.07, 76.98;
    LRAA 84.86, 84.34, 83.50, 81.46.
- Open 5' ends among no-refine ORFs matching a truth stop: complete 0 truncated vs 41,127
  UTR-extended; fragments 53,571 vs 17,922.

## Checks
- The (0.4,25) preds.tsv is byte-identical to pytd_dev_v2 on all 16 inputs.
- GM12878 StringTie (0.4,25) gff3 is byte-identical to final_abinitio's. Both GM12878 eval
  summaries are identical to the existing final_abinitio summaries.
- n_pred and tp_stop are invariant across all configs of each input (checks.tsv).
- The source run files were unchanged (stat before and after).
- Spot-checked by the main session: mix.tsv and gm12878.tsv values match the tables above.

## Caveats
- The defaults were chosen during development on these organisms and on GM12878 MANE (E09).
  Leave-one-out tests organism-level overfitting, not the choice of the 50/50 complete/fragments
  objective.
- GM12878 MANE exact given assembled cannot see wrongly trimmed truncated ORFs. The
  prediction-side share is only a partial correction.
- The fragment inputs are cut once inside every CDS, which is more truncation than real
  assemblies have.
- Only stranded inputs were run. Each cell is a single run, with no variance estimate.

## Addendum (2026-10-22): Kozak score in bits, two-score form
At Brian's request the rule was restated for readability, with no change in behaviour beyond rounding:
- coding score = C(u)/mu; UTR score = alpha*|u| + w*K, with K now in bits (log2-odds); trim when the UTR
  score is higher. Algebraically the same as C(u) < (alpha*|u| + w*K)*mu.
- The weight moves from 25 per nat to 17 per bit (25 x ln 2 = 17.33, rounded). All w values in this record
  are per nat; multiply by ln 2 = 0.693 for per-bit values: the grid's 0-100 is 0-69 per bit, and
  the flat range 10-50 is 7-35 per bit.
- Code: PyTD 44e417a, Perl ab38ca6 (see E08 for parity). Documentation:
  docs/start_codon_refinement.md in both code repos.

Check (`compare_bits.py` -> `results/bits_vs_nats.tsv`): the bits code at alpha 0.4 was run on all 16
stranded complete/fragments inputs with w = 17 and with w = 17.32868, in
`$P/work/start_param_grid/runs_bits` (GRID_DIR env var of run_grid.py), and compared with the nats
cell a0.4_w25.
- w = 17.33 (exact equivalent): 0-4 predictions differ per input, from the Kozak model file being written
  to 3 decimals; exact F1 identical to 4 decimals on 13 of 16 inputs, within 0.0001 on the rest.
- w = 17 (new default): 0-26 predictions differ per input (at most 0.13% of predictions); exact F1 within
  0.0003 of the old default on every input.
- Perl/Python regression on the 5 sample sets: identical predictions; vs the previous code, 4 of 4,277
  ORFs change (two duplicated PASA assembly pairs on the decision boundary).
