# E26. Frame contrast against the forward frames only vs all five other frames (2026-10-25)

Question (Brian): does adding the three reverse-strand frames to the frame contrast used for ORF acceptance (as before v5.0.0) change accuracy, now that acceptance uses the decoy q-value test? v5.0.0 dropped them on the strength of a smaller benchmark.

Code: PyTD branch `exp/contrast-frames` (83a4806, local worktree `work/trees/PY_frames`, not merged), = release main 3a793ff plus `--contrast-frames {forward,all}`. `all` uses C = S1 − max(S2,…,S6) instead of S1 − max(S2,S3) for both the C > 0 condition and the q-value contrast, for targets and decoys alike (decoy contrast scores the reverse complement's frames too). Start refinement and its μ unchanged. Check: `forward` reproduces release `pytd_v7` exactly on S. pombe complete.wRand (4,984 predictions); 67 tests pass.

Runs: fresh, `benchmark_v2/scripts/frames_run.sh` (same settings as `run_tool.sh pytd_v7`), launchers `data/benchmark/frames/launch2.sh` (84 inputs, `all`, q 0.05) and `launch_q.sh` (36 inputs with negatives, `all`, q 0.1 and 0.2). 156 runs, none failed. Release comparison: E25 `pytd_v7`, `pytd_v7_precise` (q 0.01), `pytd_v7_sensitive` (q 0.2).

## Claims supported
- **At the same q (0.05), all-frames is worse on 75 of 84 inputs** (exact F1 median −0.33 points, stop F1 −0.35; better on 9, at most +0.27). Medians by input class, all-frames vs release (`all_scores.tsv`):

  | inputs | stop F1 | exact F1 | stop sensitivity | stop precision |
  |---|---|---|---|---|
  | complete (stranded) | 0.976 vs 0.977 | 0.905 vs 0.907 | 0.959 vs 0.963 | 0.996 vs 0.995 |
  | complete (unstranded) | 0.975 vs 0.976 | 0.905 vs 0.906 | 0.957 vs 0.961 | 0.993 vs 0.992 |
  | complete, full UTR (stranded) | 0.969 vs 0.973 | 0.902 vs 0.904 | 0.945 vs 0.952 | 0.995 vs 0.995 |
  | complete, full UTR (unstranded) | 0.959 vs 0.962 | 0.887 vs 0.890 | 0.936 vs 0.941 | 0.988 vs 0.984 |
  | fragments (stranded) | 0.919 vs 0.928 | 0.800 vs 0.817 | 0.855 vs 0.872 | 0.994 vs 0.993 |
  | fragments (unstranded) | 0.908 vs 0.917 | 0.789 vs 0.806 | 0.838 vs 0.863 | 0.988 vs 0.983 |

  Exact F1 change by strandedness: stranded median −0.39 (44 inputs), unstranded −0.23 (40). Largest losses: T. thermophila fragments −7.7, T. thermophila fragments.wRand −3.0, mouse fragments.unstranded −2.1 points.
- **At matched rates of calls on negatives, the two rules are close; the release rule is ahead on stranded fragments.** Medians over organisms, stop sensitivity at % negatives called (`negatives_curve.tsv`):

  | inputs | release q 0.01 / 0.05 / 0.2 | all-frames q 0.05 / 0.1 / 0.2 |
  |---|---|---|
  | complete + random | 0.955@0.29 / 0.964@0.29 / 0.971@0.41 | 0.960@0.21 / 0.964@0.27 / 0.966@0.31 |
  | fragments + random | 0.837@0.17 / 0.872@0.23 / 0.895@0.45 | 0.856@0.20 / 0.869@0.24 / 0.878@0.38 |
  | complete + lncRNA | 0.950@4.3 / 0.963@5.0 / 0.974@7.8 | 0.958@4.5 / 0.965@5.4 / 0.972@7.1 |
  | fragments + lncRNA | 0.794@4.6 / 0.858@5.9 / 0.906@10.2 | 0.816@5.1 / 0.850@6.4 / 0.881@8.8 |
  | complete unstranded + random | — / 0.961@0.39 / — | 0.957@0.25 / 0.960@0.27 / 0.964@0.33 |
  | fragments unstranded + random | — / 0.859@0.32 / — | 0.836@0.28 / 0.858@0.32 / 0.866@0.49 |

  Linear interpolation of the release curve at the all-frames points: complete inputs equal within ~0.1 point; stranded fragments release ahead by ~1 point (random 0.888 vs 0.878 at 0.38%; lncRNA 0.863 vs 0.850 at 6.4%, 0.890 vs 0.881 at 8.8%); unstranded inputs equal (fragments 0.858 vs 0.859 at 0.32%) or all-frames slightly ahead (complete 0.964 at 0.33% vs 0.961 at 0.39%). The release rule needs no more than its default to reach these points; all-frames needs q 0.1 to match the release default.
- **Mechanism of the largest loss**: in stranded T. thermophila fragments, requiring S1 > all reverse frames removes 18% of candidates before the q-value test (eligible candidates with C > 0: 23,448 vs 28,726; `predict.log`). Coding fragments in this AT-rich genome with a single stop codon (TGA) often have a reverse-strand frame that scores higher than their own frame.
- **Run time**: all-frames takes a median 1.41× the release run time (both frames sets are scored for every target and decoy ORF).

Conclusion (benchmark and GM12878): no reason to restore the reverse frames. On the full benchmark with the q-value test, they lower sensitivity at the same q, are at best equal at matched negative rates (slightly better only on unstranded complete transcripts), cost ~1 point on stranded fragments and up to 7.7 points in T. thermophila, and add 40% run time. This replaces the v5 rationale in Methods.

## GM12878 (ab initio, q 0.05; all-frames vs release E25 `final_v7_abinitio`)
`genome_annotation/scripts/gm12878_frames.sh`, launcher `data/benchmark/frames/launch_gm12878_frames.sh`. Eval `eval_genome_cds.py`; proteins `score_protein_matches.py` (GENCODE v50 translations; `gm12878_frames_protein_matches.tsv`); catalogue `match_riboseq.py` (stop matches, counted as in E25).

| assembly | distinct CDSs | MANE exact (given assembled) | no GENCODE CDS overlap (all preds) | matching / non-matching proteins | catalogued Ribo-seq ORF stop matches |
|---|---|---|---|---|---|
| LRAA | 22,559 vs 22,769 | 5,810 vs 5,848 (0.923 vs 0.929) | 1.54% vs 1.50% | 32,445 / 280 vs 32,751 / 277 | 99 vs 101 |
| StringTie | 30,173 vs 30,846 | 8,124 vs 8,160 (0.911 vs 0.915) | 5.07% vs 5.47% | 34,709 / 1,017 vs 35,294 / 1,185 | 97 vs 95 |
| Trinity | 57,553 vs 59,514 | 6,814 vs 6,857 (0.888 vs 0.893) | 3.94% vs 4.25% | 75,955 / 2,317 vs 78,181 / 2,616 | 184 vs 186 |

All-frames loses matching proteins in every assembly (−306, −585, −2,226) and MANE exact models (−38, −36, −43). Non-matching proteins: +3 (LRAA), −168 (StringTie), −299 (Trinity). Lost matching per non-matching removed: StringTie 3.5, Trinity 7.4, against 3.4 for moving the release rule from q 0.05 to q 0.01 on StringTie (−1,647 / −487, E25), i.e. no better than tightening q. Catalogued Ribo-seq ORFs unchanged.
`--sensitive` (q 0.2) all-frames: LRAA 23,232 distinct CDSs, MANE exact 5,866 (0.932), no overlap 2.18%, proteins 33,176 / 492, catalogue 132; StringTie 32,750, 8,232 (0.923), 8.24%, 36,298 / 2,267, catalogue 139 (release `--sensitive`: LRAA 33,443 / 560, catalogue 144; StringTie 36,903 / 2,879, catalogue 165). At matched non-matching counts (linear interpolation between release default and `--sensitive`), the release rule gives LRAA 33,277 vs 33,176 matching proteins at 492 non-matching, and StringTie 36,322 vs 36,298 at 2,267: equal or slightly better.

## Files
`data/benchmark/frames/`: `all_scores.tsv` (score_predictions.py, release rows from E25 plus `all`), `negatives_scores.tsv`, `negatives_curve.tsv` (score_negatives.py), `score_all.sh`, launchers and job lists.
