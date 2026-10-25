# E18. Stringency sweep for a --precise mode; TD2-preprint benchmarks revisited (2026-10-23)

## Claims supported (planned; not yet in the manuscript)
- PyTD's coding-evidence threshold can be applied to ORFs of every length, not only to short ones. Doing this gives a precision/sensitivity trade-off curve, and one point on that curve is a candidate `--precise` mode. Brian asked for this on 2026-10-23: "a precision/recall plot for adding additional score thresholding on top of our current approach that will enable a --precise parameterization option".
- With the contrast threshold applied to all ORFs (`sall_c200`), the share of lncRNAs that get a prediction drops from about 19% to 6%. This costs 2 points of stop sensitivity on complete transcripts and 6 points on fragments.
- The TD2 preprint's lncRNA table measures false positives at three minimum ORF lengths; we reproduce that measurement. TD2 itself was run privately only (see `work/td2_private`), because it is not yet a paper comparator.

## Code state
All runs used the pytd_dev_v2 checkpoints at PyTD 48c870c. The final commits are PyTD 44e417a and Perl ab38ca6; they differ only in the Kozak weight (10 of 18,976 human starts change). Rerun before final numbers.

## Sources (benchmark_v2/)
- `scripts/stringency_sweep_pytd.sh`: reruns selection and refinement from the training and scoring checkpoints of the stored pytd_dev_v2 runs. Each setting is named `s<L|all>_c<C>[_noauto]`:
  - `--short-orf-length L` (`all` = 10^9, so the contrast test applies to every ORF);
  - `--short-orf-min-contrast C`;
  - `_noauto` turns off length auto-acceptance (`--retain-long-orfs-mode strict --retain-long-orfs-length 1e9`).

  The sweep has 16 settings over 20 inputs: 8 organisms × complete/fragments with random negatives, plus human and mouse × complete/fragments with lncRNA negatives. Outputs are in `work/bench_v3/stringency/`.
- `scripts/score_negatives.py`: per input and setting, reports sensitivity (stop and exact), precision, and the fraction of negatives with a prediction. Its output matches `summary.tsv` for the tools that are in both.
- `scripts/plot_stringency.py`: medians over organisms. Writes `data/benchmark/stringency/stringency_sweep.{pdf,png,medians.tsv}` and takes scores from `data/benchmark/stringency/stringency_scores.tsv`.
- `scripts/prediction_features.py`: lists features for each prediction. `data/benchmark/stringency/postfilter_rules.tsv` holds the post-filter rules that were tested on those features.
- `scripts/minlen_run.sh`: runs pytd, perl v5.7.1 or TD2 at LongOrfs `-m` 50, 98 and 150 on human and mouse wLnc inputs. Scores are in `data/benchmark/stringency/minlen_scores.no_td2.tsv`.
- `scripts/score_protein_matches.py`: scores proteins the way TD2 Tables 2 and 3 do:
  - complete = identical to a GENCODE v50 protein;
  - match30 = shares an exact 30-aa window with one;
  - everything else is nonmatching.

  The GM12878 LRAA and StringTie assemblies were scored this way; results are in `data/benchmark/stringency/gm12878_protein_matches.no_td2.tsv`. The `pytd_sall_c200` rows come from a checkpoint rerun in `work/bench_v3/stringency/gm12878/<asm>/sall_c200`.
- LNCipedia (TD2 Table 4) could not be downloaded: every lncipedia.org URL returned 404 on 2026-10-23. We used the GENCODE lncRNA negatives in the wLnc inputs instead.

## Results: medians (stop sensitivity / % negatives with a prediction)
| setting | complete.wLnc | fragments.wLnc | complete.wRand | fragments.wRand |
|---|---|---|---|---|
| defaults (s100_c200) | 0.976 / 18.9 | 0.858 / 18.7 | 0.957 / 0.48 | 0.815 / 0.11 |
| s100_c400–c1600 | 0.970 / 17.9 | 0.79 / 17.6 | | |
| **sall_c200 (candidate --precise)** | **0.957 / 6.0** | **0.795 / 5.8** | **0.935 / 0.33** | **0.767 / 0.10** |
| sall_c300 | 0.931 / 3.4 | 0.679 / 3.3 | 0.904 / 0.06 | 0.644 / 0.03 |
| sall_c400 | 0.905 / 2.7 | 0.573 / 2.6 | 0.872 / 0.03 | 0.556 / 0.02 |
| sall ≥ c1200 (only the length rule keeps ORFs) | 0.80 / 1.4 | 0.34 / 1.4 | | |
| perl v5.7.1 | 0.969 / 17.9 | 0.792 / 17.6 | 0.946 / 0.30 | 0.772 / 0.09 |

Exact sensitivity at sall_c200 is 0.887 on complete.wLnc and 0.704 on fragments.wLnc; at the defaults it is 0.906 and 0.753.

Raising the threshold for short ORFs alone saturates: lncRNA false positives are mostly ORFs of 100–300 codons that pass the length gate. Post-filters on the selected ORF did not beat simply extending the contrast test to all ORFs. The filters tried were per-base contrast rate ≥ 0.3–0.5, total contrast ≥ 200·mu, and combinations of the two (`postfilter_rules.tsv`).

## Results: minimum ORF length (TD2 Table 4 design; GENCODE lncRNAs; human / mouse)
| tool | -m 50 | -m 98 | -m 150 |
|---|---|---|---|
| PyTD final, lncRNAs with prediction % | 20.6 / 17.2 | 19.7 / 16.2 | 6.1 / 5.2 |
| PyTD complete stop sensitivity | 0.980 / 0.973 | 0.973 / 0.967 | 0.926 / 0.922 |
| perl v5.7.1, lncRNAs with prediction % | 63.6 / 62.1 | 20.6 / 17.0 | 6.1 / 5.2 |
| perl v5.7.1 complete stop sensitivity | 0.989 / 0.983 | 0.974 / 0.968 | 0.926 / 0.922 |

With its minimum lowered to 50 aa, v5.7.1 calls ORFs on about 63% of lncRNAs (the TD2 preprint reports 47.6% on LNCipedia). PyTD's contrast test for short ORFs holds this at about 20%. The remaining false positives come from ORFs of 100 codons or more.

## Results: GM12878, TD2-style protein matches (LRAA / StringTie)
| tool | predictions | complete | true positives (complete + match30) | nonmatching |
|---|---|---|---|---|
| PyTD defaults | 34,592 / 40,666 | 20,692 / 17,940 | 33,587 / 36,547 | 1,005 / 4,119 |
| PyTD sall_c200 | 32,811 / 36,364 | 20,095 / 17,500 | 32,507 / 35,221 | 304 / 1,143 |
| perl v5.7.1 | 36,861 / 49,443 | 16,808 / 16,075 | 34,699 / 39,180 | 2,162 / 10,263 |
| GMST | 34,405 / 39,980 | 16,966 / 14,048 | 33,276 / 35,060 | 1,129 / 4,920 |
| CPAT | 33,521 / 37,424 | 21,708 / 18,530 | 33,114 / 35,292 | 407 / 2,132 |
| RNAsamba | 34,531 / 38,350 | 22,107 / 18,544 | 33,746 / 35,216 | 785 / 3,134 |
| Prodigal | 51,599 / 97,217 | 14,121 / 13,766 | 38,486 / 46,440 | 13,113 / 50,777 |

On real data, sall_c200 cuts nonmatching predictions by 70% (LRAA) and 72% (StringTie). True positives fall by 3.2% and 3.6%.

The GM12878 Ribo-seq and MS support for the ORFs that a precise rule would drop is in `work/bench_v3/stringency/gm12878_orf_features_support.tsv`. Using the total-contrast ≥ 200·mu rule on the 48c870c final_abinitio runs, among ORFs with at least 10 P-sites:
- dropped ORFs have frame-supported Ribo-seq at a rate of about 0.33–0.36, against 0.62–0.66 for kept ORFs;
- some dropped ORFs are real.
