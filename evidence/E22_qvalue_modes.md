# E22. q-value cutoff of the decoy rule: default, --precise and --sensitive (2026-10-24)

## Claims supported
- The decoy rule (E20) has one user setting, the target-decoy q-value cutoff. The default is q = 0.05 (Brian). Named settings only change that cutoff:
  - `--precise` = q 0.01;
  - `--sensitive` = q 0.2;
  - `--max-qvalue Q` sets any value (`--decoy-q` is kept as an alias).
- **The default already cuts calls on lncRNAs about 3.5-fold** relative to the release rule (5.0 / 5.8% vs 18.9 / 18.7%; complete / fragments + lncRNA), with higher sensitivity on fragments. A private comparison with TD2 is in `work/td2_private/qmodes/E22_td2_notes.md`.
- **`--precise` (q 0.01) gives modest further gains.** Against q 0.05:
  - lncRNA calls 4.1 / 4.5% (vs 5.0 / 5.8%);
  - GM12878 nonmatching proteins 193 / 749 (vs 278 / 1,195; LRAA / StringTie);
  - the cost is stop sensitivity: −1.4 points on complete and −6.7 on fragment + lncRNA inputs.
  
  The remaining lncRNA calls are not short marginal ORFs. At q 0.01 in human complete.wLnc, 613 of 1,322 calls are ≥ 200 codons and 236 are ≥ 300 codons, and their number hardly changes with q (236 at q 0.05, 234 at q 0.3). No score cutoff removes them without losing true CDSs. Some may be real (lncRNA annotation is not proof of noncoding).
- **`--sensitive` (q 0.2) is where sensitivity saturates.**
  - Fragments + random: 0.896, against 0.898 at q 0.3 and 0.898 at q 0.5.
  - Calls on random sequences stay below 0.5%: 0.42% / 0.43%.
  - Calls on lncRNAs are 7.7 / 10.5%, still below the release rule (18.9 / 18.7%).
  
  It is far better than taking ORFs on length alone. The longest ORF ≥ 100 codons per transcript calls ORFs on 22% / 11% of random sequences and 26% of lncRNAs, at lower sensitivity (0.954 / 0.759). The longest ORF ≥ 50 codons calls them on 82% / 63% and 81%.
- q 0.5 is not recommended: lncRNA calls reach 15–20%, close to the release rule.

## Code
- PyTD branch `exp/short-orf-taper` (local), commit 8dbe79b: `--precise`, `--sensitive` and `--max-qvalue` are mutually exclusive. The constants are `DEFAULT_DECOY_Q`, `PRECISE_DECOY_Q` and `SENSITIVE_DECOY_Q` in `predict.py`. `tests/test_qvalue_modes.py` was added; 61 tests pass.

## Runs
- Fresh runs, never resumed (`short_orf_rule_run.sh`, setting `decoy_q<Q>_longgc`; GC-table long-ORF accept as in E20):
  - q 0.01, 0.02, 0.05, 0.1, 0.2, 0.3 and 0.5 on all 20 benchmark inputs (`jobs_qmodes.txt`, 226 jobs, none failed);
  - GM12878 at q 0.01, 0.05, 0.1, 0.2 and 0.3 (`jobs_gmq02.txt` for q 0.2).
- Baselines come from the LongOrfs output of the same runs (`cliff_c200/input.fa.transdecoder_dir/longest_orfs.{gff3,cds}`, min 50 aa):
  - `longest_orf_m50`: the longest ORF per transcript (`convert_predictions.py --longest-per-transcript`);
  - `longest_orf_m100`: the same, kept if ≥ 303 nt;
  - `longest_orf_atg`: `longest_orf_with_start.py`.
- Scores: `score_negatives.py` → `data/benchmark/qmodes/qmodes_scores.tsv`; `plot_rules.py` → `qmodes.{pdf,png,medians.tsv}`. GM12878: `score_protein_matches.py` (GENCODE v50 translations) → `gm12878_qmodes_protein_matches.tsv`. The TD2 points are in a private copy only (`work/td2_private/qmodes/`).

## Results: medians (stop sensitivity / exact sensitivity / % negatives with a prediction)
| setting | complete.wRand | fragments.wRand | complete.wLnc | fragments.wLnc |
|---|---|---|---|---|
| release rule (cliff_c200) | 0.957 / 0.897 / 0.48 | 0.815 / 0.708 / 0.11 | 0.976 / 0.906 / 18.9 | 0.858 / 0.753 / 18.7 |
| q 0.01 (`--precise`) | 0.956 / 0.888 / 0.33 | 0.836 / 0.728 / 0.19 | 0.950 / 0.881 / 4.1 | 0.790 / 0.699 / 4.5 |
| q 0.02 | 0.960 / 0.891 / 0.33 | 0.854 / 0.746 / 0.19 | 0.958 / 0.888 / 4.4 | 0.819 / 0.723 / 4.9 |
| **q 0.05 (default)** | **0.965 / 0.899 / 0.36** | **0.872 / 0.773 / 0.22** | **0.964 / 0.894 / 5.0** | **0.857 / 0.755 / 5.8** |
| q 0.1 | 0.968 / 0.903 / 0.37 | 0.888 / 0.788 / 0.31 | 0.969 / 0.898 / 5.9 | 0.884 / 0.777 / 7.5 |
| q 0.2 (`--sensitive`) | 0.970 / 0.905 / 0.42 | 0.896 / 0.799 / 0.43 | 0.975 / 0.904 / 7.7 | 0.909 / 0.797 / 10.5 |
| q 0.3 | 0.973 / 0.906 / 0.47 | 0.898 / 0.800 / 0.46 | 0.978 / 0.907 / 10.0 | 0.920 / 0.806 / 13.5 |
| q 0.5 | 0.976 / 0.908 / 0.59 | 0.898 / 0.800 / 0.60 | 0.982 / 0.910 / 15.3 | 0.929 / 0.813 / 19.7 |
| longest ORF ≥ 100 codons | 0.954 / 0.677 / 22.3 | 0.759 / 0.617 / 11.3 | 0.972 / 0.621 / 25.7 | 0.786 / 0.632 / 25.7 |
| longest ORF ≥ 50 codons | 0.988 / 0.701 / 82.5 | 0.907 / 0.760 / 62.6 | 0.988 / 0.634 / 80.9 | 0.919 / 0.736 / 80.9 |
| Perl v5.7.1 | 0.946 / 0.758 / 0.30 | 0.772 / 0.623 / 0.09 | 0.969 / 0.764 / 17.9 | 0.792 / 0.653 / 17.6 |

## Results: GM12878 (TP = complete + 30-mer match / nonmatching; LRAA | StringTie)
| setting | LRAA | StringTie |
|---|---|---|
| release rule (cliff_c200) | 33,587 / 1,005 | 36,547 / 4,119 |
| q 0.01 (`--precise`) | 32,011 / 193 | 33,869 / 749 |
| q 0.05 (default) | 32,763 / 278 | 35,391 / 1,195 |
| q 0.1 | 33,135 / 389 | 36,169 / 1,735 |
| q 0.2 (`--sensitive`) | 33,508 / 579 | 36,903 / 2,832 |
| q 0.3 | 33,657 / 706 | 37,294 / 4,186 |

## Caveats
- The q levels were chosen on this benchmark. They are rounded values, not fitted optima.
- With `--preset`, a q between stored levels uses the largest stored level ≤ q (E21). 0.01, 0.05 and 0.2 are all stored.
