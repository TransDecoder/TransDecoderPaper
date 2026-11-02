# E40. One code version throughout: pinned release commits and reruns (2026-10-31 to 2026-11-02)

Brian (2026-10-31): "we want the version of transdecoder to be consistent throughout and our latest version, so please
rerun what's needed and add Trinity". Every TransDecoder/PyTD number, figure and table now comes from the pinned commits
below at their defaults or at stated options. The only exceptions are the named earlier releases and rejected
alternatives that the released code cannot run (listed below).

## Pinned code
- PyTransDecoder `main` **d3892b5** (v0.2.0). Its parent a5b24c3 is the algorithm commit; d3892b5 only adds
  `docs/homology_te_screen.md` and a README section. `git ls-tree` gives the same `pytransdecoder/` tree for both
  (4dde030…), so runs made at a5b24c3 are runs of d3892b5. 3a793ff differs (tree 9d9a84a…): it lacks the homology
  start assist and `--no-homology-start-assist`. 76/76 unit tests pass at d3892b5.
- TransDecoder (Perl) `devel` **41fc5ee** (v7.0.0), unchanged since E25. `origin/devel` is still 41fc5ee; `master`
  (2fc944d) is the v6 line.
- Clean worktrees: `$P/work/trees/PY_rel_d3892b5`, `$P/work/trees/TD_rel_41fc5ee` (`git worktree add --detach`). No
  tags were created, and no code was changed.

## Byte identity of the earlier release runs (`pipeline/release_pin/10_byte_identity_v1.py` -> `data/release_pin/byte_identity.tsv`)
Fresh runs at d3892b5 (full LongOrfs + Predict, new directories) against the runs the tables were made from
(3a793ff, E25). Files: final GFF3, BED, CDS, pep, `preds.tsv` and LongOrfs GFF3 (benchmark); GFF3, BED, CDS, pep and
genome GFF3 (GM12878). All identical:
- benchmark default (`rel` vs `pytd_v7`): 108/108 inputs (84 benchmark + 24 mixed);
- `--precise`, `--sensitive`: 20/20 each;
- retention only (`--no-homology-start-assist` with the same DIAMOND/Pfam hits) vs `pytd_v7_homology` (3a793ff, homology
  for selection only): 84/84;
- GM12878 default (3 assemblies), `--precise` and `--sensitive` (LRAA, StringTie), Trinity retention-only vs
  `trinity/final_v7_homology`: all identical;
- the start grid's default cell `a0.4_w17` vs `rel`: 16/16 (`preds.tsv` and GFF3).
So every table built on the E25 runs (`pytd_v7` labels, `final_v7_*`) is a table of d3892b5. The run labels are kept.

## Perl/PyTD parity at the pinned commits
`regression/run_regression.sh $P/work/trees/TD_rel_41fc5ee $P/work/trees/PY_rel_d3892b5 $P/work/release_pin/regression`:
identical final ORFs on all five sample sets (737, 786, 2,584, 69, 60), now `data/regression/release_d3892b5_41fc5ee.summary.tsv`.
`regression/parity_check.py` on the same run (`data/release_pin/parity_intermediates.txt`): LongOrfs ORFs, peptides and
CDS identical; hexamer log-likelihoods within 5.1e-15. Benchmark parity (84/84) stands from E25: the Perl code is the
same commit, and the PyTD ab initio output is byte-identical (above).

## Reruns (all fresh, never resumed; launcher `pipeline/release_pin/launch_v1.sh`, 4 processes, normal priority, local)
Stage A 2026-10-31 17:40 to stage END 2026-11-02 01:58; 2,528 jobs (A-D) plus 1,264 grid runs; 0 failed. Workers:
`run_bench_v1.sh`, `run_gm_v1.sh`, `run_preset_v1.sh`, `analyses/start_param_sensitivity/run_grid_v2.py`. Outputs under
`$P/work/release_pin/` (about 80 GB): `bench/<org>/<input>/<label>`, `gm/<asm>/<cfg>`, `presets/`, `start_grid/runs/`,
`regression/`, logs in `logs/` and `launch_v1.log`.

| set | runs | options of the release code |
|---|---|---|
| benchmark default | 108 | none |
| `--precise`, `--sensitive` | 2 x 20 | (q 0.01, 0.2) |
| q sweep | 4 x 20 | `--max-qvalue` 0.02, 0.1, 0.3, 0.5 |
| minimum length without the q-value test | 6 x 20 | LongOrfs `-m` 100, 90, ..., 50; `--max-qvalue 1000000` (estimated FDRs can exceed 1) |
| long-ORF rule off | 2 x 20 | `--retain-long-orfs-mode strict --retain-long-orfs-length 1000000`, q 0.05 and 0.1 |
| homology retention only | 84 | hits of `pytd_v7_homology`, `--no-homology-start-assist` (`longest_orfs.pep` checked with `cmp`) |
| GM12878 x 3 assemblies | 55 | default, precise, sensitive, `m100_noq`, `m50_noq`, `--no-refine-starts`, alpha {0.2,0.4,0.6,1} x w {0,17,35}; Trinity retention only |
| held-out presets (E21 design) | 8 references + 728 tests | presets built with `build-preset` from fresh runs on the E21 half-A references; tests on the E21 half-B sets; trained on the input, same-species, same-species "mix", other-species |
| alpha x w grid | 1,264 | 16 inputs x (13 alpha x w {0, 3.5, 7, 17, 35, 69} per bit + `--no-refine-starts`) |

The grid and the GM12878 configs other than `rel` and `m*_noq` reuse the LongOrfs files of the `rel` run of the same
input (same commit and settings), symlinked; Predict runs in full.

**Weight mapping (old grid in nats -> release grid in bits):** 0 -> 0, 5 -> 3.5, 10 -> 7, 25 -> 17, 50 -> 35, 100 -> 69.
The default cell is exactly 17 per bit. The old grid's 25/nat was 17.33/bit.

## Scoring (scripts in `pipeline/release_pin/`, `analyses/start_param_sensitivity/`)
- `11_score_sweeps_v1.sh` -> `data/benchmark/release_pin/`: `negatives_scores.tsv` / `negatives_medians.tsv`
  (score_negatives.py; includes the longest-ORF baselines from the `rel` LongOrfs and the historical `cliff_c200`),
  `by_length.tsv`, `short_orf_summary.tsv`, `sens_by_length.tsv`, `lnc_call_lengths.human.tsv`, `lnc_riboseq_orfs.tsv`.
- `12_gm_tables_v1.py` -> `data/gm12878/release_pin/gm_runs.tsv` (eval summary, prediction-side start counts as E13,
  catalogued Ribo-seq ORF stop matches), `no_overlap_by_length.tsv`.
- `13_gm_support_v1.sh` and `14_rule_support_v1.py` -> `protein_matches.tsv`, `riboseq_summary.tsv`, `ms_summary.tsv`,
  `rule_support.tsv` (E01/E02 methods; removed ORFs = predictions with the q-value test off whose stop is not predicted
  at the default). The default rows reproduce E25's release tables exactly (protein matches 32,751/277, 35,294/1,185,
  78,181/2,616; Ribo-seq 0.657/0.615/0.612; MS unannotated 167/748/765).
- `15_score_presets_v2.sh` -> `data/benchmark/presets_release/` (v1 kept one stray row per extra organism through an
  `NR==1` filter; v2 differs only in that filter).
- `summarise_v2.py`, `open_end_mix_v2.py`, `kozak_score_distribution_v2.py` -> `data/start_param_sensitivity/release/`.
- Longest-ORF baseline: `longest_orf` rows of the benchmark (`pytd` run, 100-codon LongOrfs) equal the longest ORF
  >= 100 codons of the release LongOrfs (`m100_noq`) on all 16 wRand inputs checked (`data/release_pin/longest_orf_baseline_check.tsv`).

## Results
**Benchmark alpha x w grid (stranded complete + fragments, 8 organisms).** (0.4, 17) is still the best of the 78 cells
with the two input types weighted equally (mean exact F1 0.8477; then (0.3, 17) 0.8474, (0.5, 17) 0.8453, (0.4, 7)
0.8451). Eight cells within 0.5 points (alpha 0.2-0.5, w 7-35 per bit). Leave-one-organism-out picks (0.4, 17) in 7 of 8
folds (*T. thermophila* fold: (0.3, 17), held-out loss 0.67 points). The default is the best cell for complete weights
0.5-0.6 and within 0.13 points for 0.4-0.7 (`mix.tsv`). No default changes; none was made.

| quantity | development grid (E13; 48c870c, w per nat, fixed short-ORF test) | release grid |
|---|---|---|
| best cell, equal weights | (0.4, 25/nat) 0.8322 | (0.4, 17/bit) 0.8477 |
| LOO folds choosing default | 7/8 (Tth -> 0.3) | 7/8 (Tth -> 0.3) |
| complete-transcript optima | alpha 3 | alpha 3 |
| fragment optima | alpha 0-0.2 | alpha 0-0.3 |
| open 5' ends without refinement, complete | 41,670 of 130,066 | 41,318 of 129,778 |
| open 5' ends, fragments | 72,056 of 111,104 | 75,125 of 116,162 (75% truncated; Tth 65%, Pfa 96%) |
| refinement gain, complete | 3.1-38.7 (median 20.1) | 3.0-38.7 (median 20.1) |
| Tth / human complete, no refine -> default | 0.452->0.840 / 0.563->0.922 | 0.446->0.833 / 0.563->0.920 |
| refinement effect, fragments | -4.3 to +11.4; + in 6, Ath -0.1, Pfa 0.848->0.806 | -4.9 to +11.3; + in 7, Pfa 0.917->0.868 |
| Kozak term, fragments (vs w 0) | +1.3 to +4.3 | +1.2 to +4.3 |
| Kozak term, complete | -0.5 to +2.9 (mean +0.6) | -0.5 to +2.8 (mean +0.6) |
| Kozak at alpha 0.4, complete / fragments (means) | +0.6 / +2.7 | +0.6 / +2.8 |

**GM12878 grid, now with Trinity** (LRAA / StringTie / Trinity; `gm_runs.tsv`). The default cell equals the release run,
so the "start refinement" paragraph now quotes the release values. This removes the earlier mismatch, where the grid
gave 92.7 / 93.5 from development code with the fixed short-ORF test and the release gave 91.5 / 92.9.

| quantity | old text (E13 grid, LRAA/StringTie only) | release |
|---|---|---|
| MANE exact given assembled, no refinement -> default | ST 61.0->92.7, LRAA 59.1->93.5 | 58.6->92.9, 60.3->91.5, 63.3->89.3 |
| exact complete calls | ST 9,309->13,805; LRAA 7,562->12,018 | 7,420->11,828; 9,044->13,395; 9,622->13,380 |
| 5'-partial calls | ST 10,544->3,030; LRAA 10,693->2,631 | 10,348->2,502; 9,483->2,646; 16,006->5,931 |
| exact-start share | ST 73.3->78.4; LRAA 81.7->84.3 | 82.2->84.8; 73.9->79.0; 53.1->58.7 |
| wrong-start complete calls | ST 3,396->3,796; LRAA 1,697->2,231 | 1,609->2,113; 3,190->3,556; 8,510->9,413 |
| Kozak term, MANE exact / start share (points) | 3.1, 2.8 / 0.8, 1.1 | 2.9, 3.0, 2.1 / 1.1, 0.7, 0.8 |
| alpha 1: MANE exact / start share | ST 94.9 / 77.0; LRAA 95.5 / 81.5 | 94.9, 93.7, 91.1 / 81.9, 77.5, 56.6 |

Trinity shows the same trade-off as the other two (alpha 1 raises MANE exact to 91.1 and lowers the start share to 56.6).

**Minimum length and the q-value test (GM12878).** The comparison baseline is now the release code with the q-value test
off, not the development fixed short-ORF test. No-GENCODE-overlap share at -m 100 -> -m 50 without the test: LRAA 3.6->10.7%,
StringTie 11.8->34.7% (as the development runs), Trinity 15.4->53.8% (new). Protein matches without the test
(matching / non-matching): 34,382/3,690, 38,575/18,278, 88,821/97,397. At the default, non-matching proteins fall 13- to
37-fold for a loss of 4.7-12.0% of matching ones. The old text had 3.5-3.6-fold and 2.5-3.4% against the fixed test,
a different baseline. Removed ORFs: Ribo-seq frame test 40/37/40% vs 66/62/61% of all predictions; MS peptide 21/9/6%
vs 64/69/65%. The old text had 39/38/40% and 35/20/17% against the fixed-test runs. Catalogued Ribo-seq ORF stop
matches: 458/578/1,195 without the test, 184/221/451 at -m 100 without the test, 101/95/186 at the default, 144/165/295
with `--sensitive`. MS: unannotated predictions with a peptide 2-9% without the test (8.7/6.6/2.4%), against the
old "9-24% with the fixed test".

**Benchmark short-ORF sweeps** (`release_pin/negatives_medians.tsv`, `short_orf_summary.tsv`). The minimum-length sweep
without the q-value test reproduces the E04 development values (human complete.wLnc: 50-99 aa sensitivity 0.31 -> 0.86,
lncRNAs 6,608 -> 21,416 = 19.6 -> 63.6%, precision 0.74 -> 0.47; random calls 0.013-0.159 per sequence). The one
difference is 14,878 rather than 14,880 lncRNA calls of 50-99 aa. The q sweep: fragments + random stop sensitivity
0.897 at q 0.3 (dev 0.898), 0.898 at 0.5; human lncRNA calls >= 300 codons 236 (q 0.01), 236 (0.05), 235 (0.2),
234 (0.3). Observed false fraction with the long-ORF rule off, q 0.1: complete + lncRNA 0.090 (dev 0.089), fragments +
lncRNA 0.119 (dev 0.121), random 0.006 / 0.009 (unchanged). Fragments 50-74 codons: 0.37 / 0.53 / 0.71 at
q 0.01 / 0.05 / 0.2. Long-ORF rule effect, >= 400 codons, complete + random: 0.998 with, 0.995 without (unchanged).

**Held-out presets** (release code; `presets_release/preset_scores.tsv`). Trained on the input, <= 500 sequences: 8-82%
of random sequences called (unchanged); 1,000: 0.9% / 3.8% (was 0.8 / 3.7); 5,000: 0.4% / 0.2% (was 0.6 / 0.3).
Same-species presets: random calls 0-0.3% (was 0-0.4); mix preset stop sensitivity 0.954-0.960 complete, 0.840-0.870
fragments (was 0.948-0.960, 0.854-0.875); trained on 5,000: 0.967 / 0.878 (was 0.967 / 0.880). Complete-only preset on
fragments 0.797 vs mix 0.867 (was 0.806 vs 0.873); on complete transcripts 0.956 vs 0.957 (was 0.958 vs 0.956).
Human/mouse cross-species 0.79-0.80 (was 0.78-0.81); human -> Ath/Dmel 0.43 / 0.66 (unchanged). Conclusions unchanged.

**Homology, retention alone.** The retention-only numbers in Results (median exact-F1 gain on fragments 0.7 points;
5'-partial gain of the start rules 11.4 points over retention alone; GM12878 Trinity 7,063 MANE exact with retention
alone) were made at 3a793ff. They are now runs of the release code with `--no-homology-start-assist`, byte-identical
(above), so the numbers are unchanged. The simulated-assembly retention-only runs (E33, condition C) were made at a5b24c3
with `--no-homology-start-assist` (same code tree as d3892b5).

**Supplementary Table S6.** The "first ATG" column is now the release code with w = 0 (grid cell `a0.4_w0`; median
exact F1 0.895, was 0.893 for the development column). The PWM-choice column stays an unreleased development version,
labelled. Against the release first-ATG rule the PWM choice has fewer exact matches in 7 of 8 organisms (more in
*T. thermophila*, +296). The within-version comparison (+69 *C. reinhardtii*, +1 *P. falciparum*, lower in six) is
also stated, labelled.

## Exceptions that remain (stated in the text)
1. Earlier releases: v6.0.0 (benchmark, start comparison, Supplementary Note 1); r2012-v5.7.1 (Note 1).
2. The fixed short-ORF test (`cliff_c200`, an unreleased development version of v7): the release code has no option for
   it (`--short-orf-rule` exists only on the experimental branch). It is used as the historical comparator in the
   q-value paragraph, Fig. 4a-c,e, Supplementary Table S11 / Fig. S3, and the 50-74-codon sensitivity (0.27). Labelled
   in prose and captions.
3. The prototype with ORF-length strata only (E20; "5-6 points less" on fragments): no release option.
4. The PWM choice among later ATGs and the scale-corrected PWM (Supplementary Table S6, Results "0.728", Methods):
   development versions; no release option.
5. Experimental branches of the release code, not adopted and labelled as such: all-frame contrast (E26), UTR-length
   prior (E29).
6. Run-time ratio (Abstract "2.7 times faster", Results 2.71): from the E25 timings of Perl 41fc5ee and PyTD 3a793ff,
   whose ab initio code path equals d3892b5 (outputs byte-identical). Not re-timed.

## Checks
- Notebooks (all knitted 2026-11-02): 01 60/60, 02 173/173, 03 118/118, 04 108/108, 05 89/89, 06 56/56, 07 25/25,
  08 151/151, S1 97/97 check rows pass.
- `data/supplementary_tables/benchmark_per_input.tsv`: all non-TransDecoder rows byte-identical before/after (md5 of
  the rows without TD labels equal). The 84 homology rows are now the released homology mode (`pytd_v7h_homology`)
  instead of the 3a793ff retention-only runs. `data/benchmark/summary.release.tsv` = the release and comparator rows of
  `summary.interim_20261025.tsv`, unchanged (0 lines differ).
- The manuscript copy of Fig. 5 (`manuscript/figures/fig5_failures.pdf`) was stale (the 2026-10-26 3a793ff version while
  the repo figure and caption were the a5b24c3 version). It was synced.
- main.tex and supplement.tex compile with tectonic, no undefined references.
