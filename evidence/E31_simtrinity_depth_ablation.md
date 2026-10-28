# E31. Read-depth ablation of the simulated GM12878 Trinity benchmark (2026-10-27)

## Question
At full depth the simulated Trinity assembly of E30 is mostly complete CDSs: 83% of the scored truth is complete and only 4–8% falls in each partial class. How do the methods compare when the assembly is more fragmentary?

Approach: subsample the simulated reads to 1/2, 1/4, … of the pairs. Assemble, label, run and score each level exactly as at full depth.

## Data and design
Scripts are in `pipeline/simtrinity/`, tables in `data/simtrinity/depth/`, and the notebook section is `notebooks/02_benchmark.Rmd` › "Simulated GM12878 Trinity assembly" › "Read-depth ablation" (figures `figures/figS_simtrinity_depth.{pdf,png}`: exact F1 per CDS type with negatives counted, overall exact F1 and negatives with a prediction; `figures/figS_simtrinity_depth_sens.{pdf,png}`: exact and stop-level sensitivity per CDS type).

- **Levels.** d1 is the E30 run itself (full depth, 52,633,605 mapped pairs in the STAR BAM). d2, d4, d8, d16 and d32 keep 1/2 … 1/32 of the pairs.
- **Subsampling** (`11_depth_subsample.sh <level> <fraction>`):
  - Command: `samtools view --subsample <f> --subsample-seed 20261026` on the E30 STAR BAM (`$P/data/simtrinity/star/sim_GM12878.STAR.bam`).
  - samtools hashes the read name, so mates stay together and the levels are nested.
  - Nesting check: on chr21, all 133,877 d4 read names are among the 268,101 d2 names.
  - STAR was not re-run. Its two-pass junctions come from the full read set; they are annotated GENCODE junctions in any case, because the reads were simulated from annotated transcripts.
- **Observed fractions of pairs:**

  | level | target | pairs | observed |
  |---|---|---|---|
  | d2 | 0.5 | 26,320,005 | 0.50006 |
  | d4 | 0.25 | 13,160,759 | 0.25004 |
  | d8 | 0.125 | 6,581,285 | 0.12504 |
  | d16 | 0.0625 | 3,291,032 | 0.06253 |
  | d32 | 0.03125 | 1,646,007 | 0.03127 |

- **Per level, as in E30** (`08b_trinity_gg_level.sh`, `12_depth_label_and_run.sh`):
  - **Assembly:** Trinity 2.15.2 genome-guided with `--genome_guided_max_intron 100000 --SS_lib_type RF --min_kmer_cov 2 --max_memory 50G --CPU 4`. E30 used CPU 8; the thread count does not change the assembly.
  - **Alignment:** `minimap2 -ax splice:hq -uf --secondary=no` to GRCh38.
  - **Labelling:** `label_assembly.py` with the final E30 rules and tolerances, including the E29 statistic (`03_utr_gm12878_truth.py`).
  - **Tools:** the same nine methods via the unmodified `run_tool.sh` (input `fragments_simtrinity.fa`; stranded; human models; 4 runs at a time; no CodAn).
  - **Scoring:** `score_simtrinity.py`.
  - **Collection:** `13_depth_collect.py` gathers all levels.
- **Per-type precision** (Brian, 2026-10-27). Two versions are reported, both in `scores_by_type.tsv`.
  - *Truth only* (`precision_*`, `f1_*`): TP of type T / predictions on truth transcripts of type T.
  - *With negatives* (`*_wneg`): the denominator also includes predictions on negatives (lncRNA-derived transcripts) whose own ORF is of type T.
    - Each such ORF is typed from the transcript sequence: a start if it begins with ATG, a stop if it ends with TAA/TAG/TGA. `preds.tsv` coordinates include the stop when the tool reports one. The rule is the same for every tool.
    - Counts are in `neg_pred_types.tsv`; per type they sum to each tool's predictions on negatives.
  - The *with negatives* version applies to d1 (E30) as well. Sensitivity is unchanged.
- **Stop rule:**
  - Initial: stop when scored truth falls below ~3,000, or a partial class below ~300, or the partial counts stop rising.
  - Revised (Brian, 2026-10-27): also stop when every type's share of the scored truth changes by less than 3 points between levels.
  - Actual stop: at d32, by Brian's decision (2026-10-27), before either rule was met.
    - The shares had not stabilised: the d16 → d32 change was complete −12.2, internal +10.9 points.
    - The 5′ and 3′ partial counts had started to fall with the total.
  - d64 was started; its Trinity run was partway through the 15,163 partitions when stopped. It was killed by exact PID, not scored, and its outputs were moved to `$P/_trash/simtrinity_depth_d64_killed_20261027/`.

## Software
As in E30: Trinity 2.15.2 (`$P/env/trinity215`, ParaFly xargs drop-in), STAR index and BAM from E30, samtools 1.24, minimap2 2.31-r1302 and PyTD 3a793ff. The comparators are the same versions.

## Commands
```bash
S=$T/pipeline/simtrinity
# one lane per level, at most two Trinity runs at a time (CPU 4 each); the driver is $P/work/simtrinity/depth/run_lane_v1.sh:
bash $S/11_depth_subsample.sh d4 0.25
CPU=4 bash $S/08b_trinity_gg_level.sh d4
bash $S/12_depth_label_and_run.sh d4          # label, E29 statistic, 9 tools (4 at a time), score
# after all levels:
$P/env/td/bin/python $S/13_depth_collect.py $T/data/simtrinity/depth d1 d2 d4 d8 d16 d32
$T/scripts/knit_notebooks.sh 02_benchmark.Rmd
```
Outputs outside git: `$P/data/simtrinity/depth/<level>/` (subsampled BAM, assembly, labels, inputs) and `$P/work/simtrinity/depth/<level>/` (runs, scores). Logs are in `$P/work/simtrinity/depth/logs/`.

## Results
**Assembly and labels** (`levels.tsv`, `truth_counts.tsv`, `e29_five_prime.tsv`); counts with the % of scored truth in brackets:

| level | pairs | Trinity transcripts | complete | 5′ partial | 3′ partial | internal | truth | negatives | unresolved | E29 share |
|---|---|---|---|---|---|---|---|---|---|---|
| d1 | 52,633,605 | 84,972 | 28,187 (83.1) | 1,454 (4.3) | 1,579 (4.7) | 2,719 (8.0) | 33,939 | 20,815 | 27,666 | 7.0% |
| d2 | 26,320,005 | 70,277 | 22,864 (80.0) | 1,488 (5.2) | 1,572 (5.5) | 2,667 (9.3) | 28,591 | 15,054 | 24,469 | 8.4% |
| d4 | 13,160,759 | 56,018 | 17,188 (72.3) | 1,610 (6.8) | 1,786 (7.5) | 3,205 (13.5) | 23,789 | 9,913 | 20,319 | 12.4% |
| d8 | 6,581,285 | 44,871 | 11,615 (57.6) | 1,926 (9.6) | 1,973 (9.8) | 4,646 (23.0) | 20,160 | 5,885 | 16,963 | 20.8% |
| d16 | 3,291,032 | 36,052 | 7,039 (40.3) | 2,096 (12.0) | 1,920 (11.0) | 6,402 (36.7) | 17,457 | 3,267 | 13,529 | 31.6% |
| d32 | 1,646,007 | 26,744 | 3,889 (28.1) | 1,839 (13.3) | 1,531 (11.1) | 6,587 (47.6) | 13,846 | 1,594 | 9,634 | 43.9% |

- Shares are computed from the counts; `share_*` in `truth_counts.tsv` is written to 6 decimals.
- "E29 share" is the fraction of transcripts carrying a MANE Select stop, with the reference structure to their 5′ end, that begin inside the CDS. The real GM12878 Trinity assembly has 13.6% (E29), between d4 and d8.

**Overall** (`scores_overall.tsv`, `negatives.tsv`), with values for d1 / d2 / d4 / d8 / d16 / d32.

Stop F1:

| method | d1 | d2 | d4 | d8 | d16 | d32 |
|---|---|---|---|---|---|---|
| TransDecoder (pytd_v7) | 0.916 | 0.920 | 0.923 | 0.924 | 0.929 | 0.925 |
| longest ORF | 0.779 | 0.808 | 0.841 | 0.879 | 0.914 | 0.935 |
| Prodigal | 0.630 | 0.674 | 0.713 | 0.771 | 0.844 | 0.892 |
| GeneMarkS-T | 0.873 | 0.886 | 0.891 | 0.882 | 0.857 | 0.821 |
| CPAT | 0.901 | 0.902 | 0.900 | 0.873 | 0.838 | 0.806 |
| RNAsamba | 0.900 | 0.904 | 0.904 | 0.888 | 0.865 | 0.839 |

Exact F1:

| method | d1 | d2 | d4 | d8 | d16 | d32 |
|---|---|---|---|---|---|---|
| TransDecoder | 0.844 | 0.846 | 0.846 | 0.832 | 0.825 | 0.806 |
| TransDecoder v6.0.0 | 0.671 | 0.675 | 0.688 | 0.687 | 0.683 | 0.664 |
| longest ORF | 0.571 | 0.585 | 0.606 | 0.651 | 0.716 | 0.763 |
| Prodigal | 0.453 | 0.478 | 0.506 | 0.566 | 0.656 | 0.725 |
| GeneMarkS-T | 0.687 | 0.697 | 0.697 | 0.706 | 0.706 | 0.694 |
| CPAT | 0.799 | 0.794 | 0.762 | 0.672 | 0.538 | 0.428 |
| RNAsamba | 0.787 | 0.780 | 0.745 | 0.652 | 0.512 | 0.404 |
| CPC2 | 0.781 | 0.776 | 0.742 | 0.653 | 0.512 | 0.389 |
| Borf | 0.723 | 0.729 | 0.718 | 0.682 | 0.611 | 0.547 |

Negatives with a prediction (%):

| method | d1 | d2 | d4 | d8 | d16 | d32 |
|---|---|---|---|---|---|---|
| TransDecoder | 8.1 | 8.9 | 10.0 | 12.0 | 14.3 | 16.6 |
| longest ORF | 70.6 | 68.3 | 67.7 | 66.2 | 65.5 | 69.1 |
| Prodigal | 51.8 | 48.6 | 49.3 | 48.5 | 48.3 | 48.4 |
| GeneMarkS-T | 20.8 | 19.1 | 18.3 | 17.5 | 16.8 | 16.7 |
| CPAT | 11.2 | 11.4 | 11.3 | 11.8 | 12.7 | 13.6 |
| RNAsamba | 11.7 | 12.5 | 13.3 | 14.1 | 17.2 | 20.3 |
| CPC2 | 7.4 | 6.9 | 6.7 | 6.3 | 5.8 | 5.3 |
| Borf | 19.7 | 18.3 | 17.7 | 16.9 | 15.7 | 15.6 |

**By type: exact F1 with negatives, TransDecoder vs the best other method** (`scores_by_type.tsv`):

| type | d1 | d2 | d4 | d8 | d16 | d32 |
|---|---|---|---|---|---|---|
| complete | 0.869 vs CPAT 0.888 | 0.875 vs CPAT 0.899 | 0.884 vs CPAT 0.913 | 0.886 vs CPAT 0.923 | 0.894 vs CPAT 0.932 | 0.891 vs CPAT 0.941 |
| 5′ partial | **0.599** vs GMST 0.455 | **0.587** vs GMST 0.536 | **0.639** vs GMST 0.611 | 0.676 vs GMST 0.707 | 0.724 vs GMST 0.749 | 0.732 vs Prodigal 0.808 |
| 3′ partial | 0.776 vs CPAT 0.796 | 0.816 vs RNAsamba 0.826 | 0.826 vs CPAT 0.858 | 0.837 vs RNAsamba 0.893 | 0.845 vs RNAsamba 0.904 | 0.831 vs RNAsamba 0.904 |
| internal | 0.747 vs Prodigal 0.872 | 0.743 vs Prodigal 0.881 | 0.742 vs Prodigal 0.905 | 0.748 vs longest ORF 0.951 | 0.771 vs longest ORF 0.973 | 0.768 vs longest ORF 0.982 |

- **Truth-only precision** (negatives not counted): TransDecoder's exact F1 at d1 / d32 is complete 0.887 / 0.900, 5′ 0.663 / 0.744, 3′ 0.823 / 0.849, internal 0.765 / 0.773.
  - The best other method on 5′ partials and internals is then the longest-ORF baseline (5′ partials 0.875–0.908; internals 1.000 by construction).
- **Exact sensitivity, TransDecoder, d1 → d32:** complete 0.859 → 0.876, 5′ 0.574 → 0.680, 3′ 0.751 → 0.794, internal 0.689 → 0.719.
- **Predictions on negatives by ORF type** (`neg_pred_types.tsv`), at d1 / d32:
  - TransDecoder: 1,688 / 264 predictions, of which 5′-open 267 / 56.
  - Prodigal: 15,842 / 905, of which 5′-open 7,074 / 414.
  - Longest ORF: 14,687 / 1,102, of which 5′-open 6,651 / 498.

**Observations**
- **Fewer reads make a more fragmentary assembly.** At d4–d8 the 5′-truncation rate (E29 share 12.4–20.8%) brackets the real GM12878 Trinity assembly (13.6%). At d32 internal fragments are the largest class.
- **TransDecoder's overall stop F1 stays at 0.916–0.929 across all levels.** Its overall exact F1 falls from 0.844 to 0.806, because the mix shifts to partials, on which exact matches are harder. Per type its exact F1 rises with fragmentation (5′ partial 0.599 → 0.732).
- **The ATG-requiring classifiers (CPAT, RNAsamba, CPC2) fall steeply** (exact F1 0.80 → 0.43 for CPAT). They cannot report 5′-open ORFs: exact sensitivity ≤ 0.02 on 5′ partials and internals at every level.
- **Prodigal and the longest-ORF baseline improve as partial CDSs dominate** (exact F1 0.45 → 0.73 and 0.57 → 0.76). But they predict on about half and about two-thirds of the negatives at every level. With negatives counted, Prodigal's exact F1 on 5′ partials rises from 0.259 (d1) to 0.808 (d32) and passes TransDecoder from d16 (0.739 vs 0.724). GeneMarkS-T passes TransDecoder on 5′ partials from d8. Prodigal and the longest-ORF baseline lead on internals at every level.
- **CPAT and RNAsamba remain best on complete and 3′ partial CDSs** (start codon present) at every level.

## Checks
- Observed subsampling fractions are within 0.0002 of the target at every level. Nesting was checked on chr21 (above).
- At each level the labels partition the input FASTA; `score_simtrinity.py` asserts this.
- `neg_pred_types.tsv`: the per-type counts sum to `n_pred_on_negatives` of `negatives.tsv` for every level × tool. Example: d1 pytd_v7 1,688; prodigal 15,842.
- d1 rows reproduce the E30 tables (`data/simtrinity/scores_overall.tsv` etc.).
- After the knit, notebook 02 still shows 0 of 133 mismatches.

## Caveats
- **One subsample per level, one simulation.** There are no replicates, so differences of about 0.01 are not interpretable.
- **STAR was not re-run per level.** The junction set and multi-mapping choices come from the full read set.
- **Trinity used CPU 4 for the levels and CPU 8 for d1.** Per-partition assembly does not depend on the thread count.
- **NFS lock stall (d32).** One partition (chr3:10028680–10101908, + strand) stalled for about 1.5 h at the end of its Butterfly step. The ParaFly xargs drop-in was waiting on an NFS `flock` on its `butterfly_commands.completed` file; all its Butterfly commands had finished.
  - The two stalled `flock` processes were killed by PID (1561509, 2758256). The completion records were written and the partition's output (8 transcripts) was collected.
  - No other partition was affected: 18,352 of 18,352 completed.
- **Per-type precision with negatives types each negative-transcript ORF by its own ends.** A tool that trims a 5′-open ORF to an internal ATG is counted under complete or 3′ partial, not under 5′ partial.
- **Stop.** The ablation stopped at d32 by decision, before the share-stability rule was met (see Design). d64 was killed unscored.
