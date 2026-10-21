# E08. Perl / Python implementation concordance and runtime

## Claims supported
Manuscript: Methods, "Implementations" (`methods_algorithm.tex`, §sec:impl); Supplementary Table S2 (`\todo`).

1. PyTransDecoder reproduces the Perl algorithm, including tie order.
   - On five test data sets the two give identical ORFs, coding scores and pre-refinement selections.
   - At the time of the regression they differed only in start codons chosen by the randomised PWM training.
2. **New, stronger claim:** at the final defaults (frame-contrast + Kozak starts, deterministic), Perl and Python predictions are identical on all 84 benchmark inputs (8 organisms × 10 variants, plus 4 lncRNA inputs).
3. PyTD runs about 2.7–2.9× faster than Perl on the regression sets and 2.3× faster (median) on the benchmark inputs. The manuscript `\todo{2.7}` should be updated; see Results.
4. Before the 2026 fixes (v6.0.0 vs PyTD 43f4293), the implementations differed in ORF selection as well as starts.

## Data
- **Regression sets:** the 5 sample datasets shipped in `TransDecoder/sample_data`. They are extracted by `run_regression.sh` into `$P/work/regression/<run>/inputs/`:
  - Trinity (`simple_transcriptome_target/Trinity.fasta.gz`);
  - PASA (`pasa_example`);
  - StringTie (`stringtie_example`; run with `-S`);
  - SuperTranscripts (`supertranscripts_example`);
  - Cufflinks (`cufflinks_example`, cDNA via `gtf_genome_to_cdna_fasta.pl`).
- **Benchmark inputs:** bench_v3 (E07), `$P/data/bench_v3/inputs/`.

## Software
- **Regression runs** (directories under `$P/work/regression/`; commits inferred from run time against tree commit times):

  | run | Perl tree | PyTD tree | run time | notes |
  |---|---|---|---|---|
  | `baseline_v6.0.0` | TD_base 66d4712 (v6.0.0) | PY_base 43f4293 | 2026-10-17 16:14 | before fixes |
  | `dev_paper_r4` | TD_fix 8cf8d5f (fix/start-refinement) | PY_fix (dev/paper) at about dc7b514, committed 16:43:08 | 2026-10-17 16:43 | after parity fixes |
  | `dev_paper_r6` | TD_fix 694e204 | PY_fix f5b12e7 | 2026-10-17 18:39–18:40 | "final" code of 10-17 (notes/04) |

  - The r6 Perl Predict log shows commands from `work/trees/TD_fix/util/`, and the PyTD predict log has timestamps 2026-10-17 18:40.
  - TD 694e204 and PyTD f5b12e7 were both committed at 18:38:18.
- **Benchmark runs** are the E07 labels:
  - `pytd` (PY_fix 27bfb5e) vs `perl_td` (TD_fix 351634a), PWM start refinement, 10-18;
  - `pytd_dev_v2` (PY_dev 48c870c) vs `perl_td_dev_v2` (TD_fix de3f57e), final defaults, 10-20.
- **Environment:** `$P/env/td` (Python 3.11.17, Biopython 1.88, Perl 5.32.1); `PERL_HASH_SEED=0`.

## Commands
```bash
cd $P
# regression (script: $T/regression/run_regression.sh; comparison: $T/regression/compare_orfs.py)
PYTHON=env/td/bin/python TransDecoderPaper/regression/run_regression.sh work/trees/TD_base work/trees/PY_base work/regression/baseline_v6.0.0
PYTHON=env/td/bin/python TransDecoderPaper/regression/run_regression.sh work/trees/TD_fix  work/trees/PY_fix  work/regression/dev_paper_r6
cp work/regression/<run>/summary.tsv TransDecoderPaper/regression/results/<run>.summary.tsv
# benchmark concordance: E07 commands (launch16.sh), then compare preds.tsv per input
for f in work/bench_v3/runs/*/*/pytd_dev_v2/preds.tsv; do
  cmp -s <(cut -f1-4 $f | sort) <(cut -f1-4 ${f/pytd_dev_v2/perl_td_dev_v2} | sort) || echo DIFF $f; done
```
- **What `run_regression.sh` does:**
  - It runs LongOrfs and Predict at default settings with no homology. Each implementation runs sequentially on the same host, timed by wall clock.
  - `compare_orfs.py` keys ORFs by (transcript, strand, CDS span). It reports exact matches, same-stop-different-start, implementation-only ORFs and Jaccard indices.
- **Reconstruction:** the exact invocations of the regression runs were not saved; the commands above are reconstructed (see Gaps).

## Outputs
- `$T/regression/results/baseline_v6.0.0.summary.tsv`, `dev_paper_r4.summary.tsv` and `dev_paper_r6.summary.tsv`.
- Run directories (not in git): `$P/work/regression/{baseline_v6.0.0,dev_paper_r1..r4,dev_paper_r6,fix_start_refinement}/<dataset>/{perl,py}/`, with `discordant.tsv` and `compare.txt` per dataset.
  - Runs r1–r3, r5 and `fix_start_refinement` were intermediate development iterations.
- Benchmark summaries: `$T/benchmark_v2/results_v3/summary.tsv` (E07).

## Results
**Regression, final ORFs** (`dev_paper_r6.summary.tsv`; r4 is identical except for runtime):

| dataset | ORFs (Perl = Py) | exact match | same stop, different start | Perl-only / Py-only | Perl s | Py s | Perl/Py |
|---|---|---|---|---|---|---|---|
| Trinity | 679 | 667 | 12 | 0 / 0 | 80.2 | 28.6 | 2.80 |
| PASA | 792 | 739 | 53 | 0 / 0 | 77.3 | 27.4 | 2.82 |
| StringTie | 3,037 | 2,940 | 97 | 0 / 0 | 129.2 | 47.2 | 2.74 |
| SuperTranscripts | 117 | 117 | 0 | 0 / 0 | 71.4 | 25.6 | 2.79 |
| Cufflinks | 80 | 76 | 4 | 0 / 0 | 44.8 | 15.5 | 2.90 |

Before fixes (`baseline_v6.0.0.summary.tsv`):

| dataset | Perl ORFs | Py ORFs | Perl-only / Py-only | same stop, different start |
|---|---|---|---|---|
| Trinity | 679 | 682 | 1 / 4 | 15 |
| PASA | 792 | 792 | 0 / 0 | 67 |
| StringTie | 3,037 | 3,045 | 22 / 30 | 71 |
| SuperTranscripts | 117 | 117 | 0 / 0 | 0 |
| Cufflinks | 80 | 80 | 0 / 0 | 4 |

**Intermediate files, r6.** These were re-checked on 2026-10-21 from the run directories:
- **LongOrfs:** the `longest_orfs.pep` IDs, order and sequences are identical in all 5 datasets (845 / 1,067 / 4,353 / 268 / 95 ORFs). `longest_orfs.cds` is identical ignoring case. Case differs only in StringTie, where Perl keeps soft-masked lowercase.
- **Hexamer model:** `hexamer.scores` has the same k-mers, and log-likelihoods agree to ≤ 5.1×10⁻¹⁵.
- **Six-frame scores:** `longest_orfs.cds.scores` is identical after sorting.
- **Pre-refinement selections:** the selected ORF IDs in `best_candidates.gff3` are identical (679 / 792 / 3,037 / 117 / 80).

**Benchmark (E07, `summary.tsv` and `preds.tsv`):**
- **Final defaults** (`pytd_dev_v2` vs `perl_td_dev_v2`):
  - `preds.tsv` coordinates are identical in **84/84** inputs, so every metric in `summary.tsv` is identical.
  - Runtime ratio Perl/Py: median 2.29 (range 1.99–2.51) over 84 runs.
- **PWM-era code** (`pytd` vs `perl_td`, 10-18):
  - tp_exact differs in 78/80 inputs (max |Δ| 3,711, from start choice).
  - tp_stop and n_pred differ in 10/80 (max |Δ tp_stop| 55).
  - Runtime ratio median 2.39 (1.71–2.84).

## Checks
- **Parity of intermediates:** done above with a short script (`$T/regression/parity_check.py`). It parses the FASTA, hexamer and GFF3 files rather than comparing bytes; the formats and float printing differ.
- **Final-default reproducibility:** `pytd_dev_v2` and `perl_td_dev_v2` were rerun at the committed trees (48c870c / de3f57e) on spombe/complete and creinhardtii/complete. Both give `preds.tsv` byte-identical to the benchmark runs (E07).
- **PWM seed noise:** the r6 start differences (12 / 53 / 97 / 4) are within the PyTD seed-to-seed variation of the PWM refinement: PASA 29–61 and StringTie 47–144 starts differ between seeds (E09, `analyses/start_refinement_stability`).

## Gaps
- **Regression not rerun at final code:** the regression has not been run at the final code (PY_dev 48c870c / TD_fix de3f57e). The committed r6 table reflects the 10-17 PWM-refinement code. The "0 Perl-only / 0 Py-only, starts differ" statement in the manuscript describes that code. At the final defaults the stronger evidence is the 84/84 benchmark identity.
- **Invocations not saved:** the exact `run_regression.sh` invocations, and which commits r1–r5 used, were not recorded. Commits are inferred from timestamps. For r4, the PyTD Predict run started 4 s after dc7b514 was committed, so its exact code state is uncertain.
- **Parity claim only in notes:** "the hexamer model and all 6-frame scores of every ORF are identical" (notes/03) had no saved comparison. It was verified here for r6, but the script is not committed.
- **Runtime:** "Python is about 2.7× faster" (notes/03) comes from the summary tables (Perl/Py ratios 2.74–2.90 in r6 and 2.66–2.92 in r4). The manuscript's `\todo{2.7}` should cite these numbers or the benchmark median 2.29. No dedicated, unloaded runtime benchmark exists.
- **Superseded claims in notes/02:** the "10× faster" claim and the "Perl 679 vs Python 682" parity were measured before PWM refinement was ported. Do not use them.

## Caveats
- **Wall-clock timing.** Regression runs were sequential per dataset. Benchmark runs ran 6–16 at a time on a shared host (load ~40). Ratios are more reliable than absolute times.
- **Different modes.** The regression runs Perl via `util/TransDecoder.LongOrfs`/`Predict`, not the v6 wrapper. The benchmark runs either `TD/TransDecoder.LongOrfs` or `TD/util/...`, as `run_tool.sh` resolves them.
- **Unpushed code.** The final-default code (48c870c, de3f57e) is on unpushed local branches (E07).
- **Manuscript out of date.** `methods_algorithm.tex` §Implementations still says the two implementations "differ only in start codons chosen by the randomised PWM training". At the final defaults they do not differ at all on the benchmark.
