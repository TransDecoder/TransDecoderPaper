# E33. Keep-open vs trim for 5'-open ORFs: how much better can a classifier decide? (2026-10-28)

Status: **exploratory, not in the paper.** Step 1 (oracle-label ceiling test) and step 2 (self-trained
prototypes, approved by Brian 2026-10-28) done. No TransDecoder/PyTD code was changed.

## Question
The release start refinement (methods_algorithm.tex "Start-codon refinement"; E09, E13, E29) takes an
ORF that runs off the 5' end, finds its first in-frame ATG, lets u = the segment before it, and trims
to the ATG when 0.4|u| + 17K > C(u)/mu. On simulated Trinity assemblies (E30/E31) it finds the stop of
~88% of 5' partials but the exact ORF of only ~59%. Can rules or a classifier using |u|, C(u)/mu, the
Kozak score K, transcript length, ORF length, 3' UTR length etc. decide better? Do the shuffled decoys
give a usable null for the UTR class?

## Data
`$P = /broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T = $P/TransDecoderPaper`.
All runs are existing, finished PyTD release runs (PyTD 3a793ff, `pytd_version.txt` checked by the
driver); nothing was resumed or rerun.
- Benchmark (E07/E28), stranded, 8 organisms x {complete, complete_fullutr, fragments,
  fragments_mixed_t1/t2/t3} = 48 inputs: runs `$P/work/bench_v3/runs/<org>/<input>/pytd_v7/`,
  truth `$P/data/bench_v3/inputs/<org>/<input>.truth.tsv`.
- Simulated human Trinity assemblies (E30/E31): d1 run `$P/work/simtrinity/runs/pytd_v7/`, truth
  `$P/data/simtrinity/inputs/fragments_simtrinity.truth.tsv`; d2, d4, d8, d16, d32 runs
  `$P/work/simtrinity/depth/<d>/runs/pytd_v7/`, truth `$P/data/simtrinity/depth/<d>/inputs/`.
- Simulated Arabidopsis Trinity assembly d1 (other agent's E32 work, finished pytd_v7 run, read only):
  `$P/work/simtrinity_athal/depth/d1/runs/pytd_v7/`, truth
  `$P/data/simtrinity_athal/depth/d1/inputs/fragments_simtrinity.truth.tsv`.
- Each run's `input.fa` was checked identical (`cmp`) to the input FASTA for the simulated assemblies
  and one mixed input.

## Software
- PyTD 3a793ff (`$P/work/trees/PY_v7`), imported on PYTHONPATH for every score (frame contrast, mu,
  Kozak model and window, dinucleotide shuffle, Xoshiro128ss PRNG).
- New env `$P/env/startclf` (`micromamba create -p $P/env/startclf -c conda-forge python=3.11
  scikit-learn pandas numpy scipy biopython tqdm matplotlib`; log `$P/env/startclf.install.log`):
  Python 3.11.18, scikit-learn 1.9.1, pandas 3.1.0, numpy 2.4.6, scipy 1.17.1, biopython 1.88.
  (`$P/env/td` has no scikit-learn; it was not modified.)
- Paper repo at 0483a09 (+ uncommitted E33 files).

## Method (scripts `$T/pipeline/start_classifier/`, tables `$T/data/start_classifier/`)
1. `01_features.py` (driver `01_features_all_v1.sh`): for every ORF selected before refinement
   (`best_candidates.gff3` in the run's work dir; refinement only moves starts, so this is the
   `--no_refine_starts` selection) that is 5'-open and has an in-frame ATG, recomputes with PyTD's
   code |u| (`seg`), C(u)/mu (`coding`), per-base `coding_pb`, K (`kozak`), Kozak window positions off
   the transcript (`kozak_off`), the rule's UTR score and decision, and reads the run's decision from
   `best_candidates.gff3.revised_starts.gff3`. Context: transcript length, ORF length, ORF length /
   median selected ORF length of the input (`orf_len_rel`), ORF lacks a stop (`internal`), 3' UTR
   length, GC of u.
   Labels (only candidates whose 3' end matches a truth CDS within 2 nt on the same strand):
   `truncated` = truth 5prime_partial/internal (correct: keep open); `utr` = truth complete/3prime_partial
   (correct: trim), with `true_atg_rank` 1 if the true start is the first in-frame ATG, >= 2 if later
   (then neither decision is exact), 0 if inconsistent (dropped).
   Decoy null `dec_coding`: the transcript is dinucleotide-shuffled exactly as PyTD's decoy k = 0
   (`decoys.dinucleotide_shuffle`, seed `42:<tid>:0`); a random stop-free in-frame stretch of length |u|
   (random frame and position, Python `random.Random('E33:<dataset>')`) is scored as C/mu.
2. `02_separation.py` -> `step1/separation.tsv`, `candidates_by_dataset.tsv`, `lr_coefficients.tsv`,
   `coding_score_by_seg.tsv`, `coding_score_scaling_fit.tsv`; stdout `step1.stdout.txt`.
   Models (p = P(truncated), keep open when p >= 0.5): `rule` (release rule, fixed);
   `rule_lr` (logistic regression on |u|, C, K); `lr_intr` (LR on intrinsic u/ATG features: |u|,
   log|u|, C, C/|u|, C/sqrt|u|, K, kozak_off, GC(u)); `lr_full` (+ log tx length, log ORF length,
   log ORF length rel., internal, log 3' UTR); `hgb_intr`, `hgb_full` (sklearn
   HistGradientBoosting, 300 iterations, lr 0.05, 31 leaves, l2 1). Schemes: `cv_bench` (pooled 48
   inputs, 5-fold GroupKFold by organism+transcript, since complete/fragments/mixed share
   transcripts); `looo` (leave one organism out); `bench2sim` (train on all benchmark, test each
   simulated assembly; `bench2sim_oracle_prior` additionally shifts the logit to the test set's true
   truncated share); `cv_sim` (5-fold CV within each simulated assembly, grouped by contig).
   Metrics: AUC; accuracy at the decision; best accuracy over thresholds on the test set (oracle
   threshold); exact = decisions that yield the exact ORF; `exact_max` = truncated + utr with rank 1.
   |u| bins 0-19, 20-49, 50-99, 100-199, 200-499, >= 500.
3. `03_ablation.py` -> `step1/lr_ablation.tsv` (stdout `step1_ablation.stdout.txt`): LR feature-set
   ablation under cv_bench, cv_sim and bench2sim.

## Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; cd $T
bash pipeline/start_classifier/01_features_all_v1.sh          # ~25 min serial; log data/start_classifier/features/extract.log
$P/env/startclf/bin/python pipeline/start_classifier/02_separation.py data/start_classifier/features data/start_classifier/step1 > data/start_classifier/step1.stdout.txt
OMP_NUM_THREADS=3 $P/env/startclf/bin/python pipeline/start_classifier/03_ablation.py data/start_classifier/features data/start_classifier/step1/lr_ablation.tsv
```
(02 ran HistGradientBoosting with default OpenMP threads, 6.6 min wall / 43 min CPU; later runs set
OMP_NUM_THREADS=3 to stay within the shared-machine budget.)

## Checks
- Re-derived rule decisions equal the runs' own decisions in all 55 runs (`trim_mismatch_vs_run` = 0
  everywhere; `features/extract.log`, `step1/candidates_by_dataset.tsv`).
- Inconsistent utr labels (true start upstream of the first ATG / not an in-frame ATG) are rare:
  0-74 per input (tthermophila fragments 74 of 15,074), dropped.

## Results (step 1)
### Candidates (`candidates_by_dataset.tsv`)
- 378,486 labelled benchmark candidates; 63,693 in the 7 simulated assemblies.
- Truncated share pi: benchmark complete 0; fragments and mixed 0.63-0.95 (pfalciparum 0.95);
  simulated human 0.27 (d1), 0.29, 0.37, 0.50, 0.65, 0.72 (d32); simulated athal d1 0.33.
- UTR candidates whose true start is a later ATG (no exact answer possible): 1-3% except
  creinhardtii (~30% of utr candidates; E09) and human sim d1 5%.

### Separation (AUC; accuracy at the decision), `separation.tsv`
| test | rule | rule_lr | lr_intr | lr_full | hgb_intr | hgb_full |
|---|---|---|---|---|---|---|
| benchmark pooled, cv | 0.899; 0.818 | 0.901; 0.822 | 0.905; 0.828 | 0.936; 0.862 | 0.933; 0.853 | 0.974; 0.913 |
| sim human d1, cv within | 0.930; 0.868 | 0.930; 0.880 | 0.952; 0.905 | 0.972; 0.931 | 0.965; 0.924 | 0.980; 0.942 |
| sim human d32, cv within | 0.911; 0.822 | 0.918; 0.864 | 0.942; 0.891 | 0.959; 0.914 | 0.952; 0.902 | 0.962; 0.918 |
| sim athal d1, cv within | 0.817; 0.763 | 0.858; 0.809 | 0.869; 0.816 | 0.947; 0.908 | 0.866; 0.817 | 0.949; 0.914 |
| sim human d1, trained on benchmark | 0.930; 0.868 | 0.928; 0.853 | 0.921; 0.829 | 0.926; 0.876 | 0.942; 0.864 | 0.898; 0.897 |
| sim human d32, trained on benchmark | 0.911; 0.822 | 0.917; 0.850 | 0.906; 0.839 | 0.914; 0.790 | 0.933; 0.873 | 0.806; 0.583 |
| sim athal d1, trained on benchmark | 0.817; 0.763 | 0.840; 0.756 | 0.842; 0.759 | 0.854; 0.800 | 0.863; 0.729 | 0.900; 0.874 |

- Leave-one-organism-out (benchmark), AUC rule vs lr_intr vs hgb_full: athal 0.838/0.852/0.926,
  creinhardtii 0.908/0.912/0.939, dmel 0.919/0.929/0.975, human 0.918/0.902/0.977, mouse
  0.907/0.898/0.977, pfalciparum 0.906/0.914/0.972, spombe 0.882/0.883/0.966, tthermophila
  0.906/0.898/0.920. But hgb_full's LOOO accuracy on complete inputs (median 0.986 vs rule 0.820) shows
  it recognises the benchmark input type from context (complete inputs: long 3' UTRs, full-length
  transcripts), a construction artefact, not a transferable signal.
- Exact ORFs over these candidates, sim human d1: rule 8,506 of a possible 9,710; in-domain CV
  hgb_full 9,150 (+644), lr_full 9,039 (+533), lr_intr 8,809 (+303); athal d1: rule 5,747 / 7,532,
  hgb_full 6,888 (+1,141), lr_full 6,843; human d32: 6,893 / 8,369, hgb_full 7,696 (+803).
  Benchmark-trained models do not reliably gain on the simulated assemblies (human d1: lr_full 8,553,
  hgb_full 8,679; d32: hgb_full 4,873, lr_full 6,650 vs rule 6,893).
- Oracle prior shift (true pi of the test set) changes bench2sim accuracy by <= 0.04 (largest for
  athal); it does not rescue the context models on d32 (hgb_full 0.643).
- By |u| (cv_sim d1, AUC rule -> hgb_full): 0-19 0.797 -> 0.918; 20-49 0.900 -> 0.970; 50-99 0.932 ->
  0.980; 100-199 0.952 -> 0.989; 200-499 0.951 -> 0.985. The short-|u| bins gain most, through
  context (intrinsic hgb_intr 0.872 at 0-19).

### Which features carry it (`lr_ablation.tsv`, LR, AUC cv_sim human d1 / athal d1 / bench2sim human d1)
- Rule features alone (|u|, C, K): 0.930 / 0.858 / 0.928. Adding C/sqrt|u| or log|u|: no change
  (<= 0.003). The rule's linear form already uses the u-intrinsic coding and Kozak information.
- + GC of u: 0.951 / 0.868 in domain, but 0.920 when transferred (human UTRs GC-richer than CDS here,
  not on the benchmark mix): an input-specific signal.
- Single context features on top of intrinsic: + log transcript length 0.968 / 0.929 / 0.948 (also
  transfers: bench2sim athal 0.842 -> 0.895); + 3' UTR length 0.969 / 0.941 / 0.934; + internal
  0.969 / 0.936 / 0.930; + log ORF length 0.965 / 0.920 / 0.942. All together (lr_full) 0.972 / 0.947
  / 0.926: combined, they overfit the benchmark's construction when transferred.

### Coding score vs |u|, and the decoy null (`coding_score_by_seg.tsv`, `coding_score_scaling_fit.tsv`)
- Mean C(u)/mu per base: truncated 0.86 (benchmark), 0.92 (human d1), 0.77-0.88 (human d2-d32),
  0.63 (athal d1); UTR 0.007 (benchmark), -0.02 to -0.04 (human sims), 0.09 (athal d1).
- Spread grows faster than sqrt|u|: sd ~ |u|^0.73 (truncated), |u|^0.60-0.64 (UTR) on every set
  (benchmark sd at 200-499 nt: truncated 147, UTR 102).
- Decoy null (PyTD dinucleotide shuffles, stop-free, same length): median and IQR close to real
  UTRs at |u| < 200 (benchmark KS D = 0.026-0.044 per bin; median -5.0 vs -5.8 at 50-99 nt; human d1
  KS 0.03-0.08 at < 100 nt). It is narrower in the tail: sd exponent 0.52-0.59 vs 0.60-0.64, benchmark
  200-499 nt sd 78 vs 102; KS at 100-499 nt 0.03-0.08 on the benchmark but 0.07-0.15 on the
  simulated assemblies (0.03-0.12 below 100 nt); the real UTR tail toward positive coding scores
  (UTR-located uORFs/coding-like UTRs) is missing. Stop-free decoy stretches of >= 500 nt are rare
  (28 of 103 benchmark candidates). For truncated candidates the decoy is far from the real segment
  (KS 0.21-0.93), as expected.
  -> The decoys approximate the UTR class well where most decisions fall (|u| < 200 nt) and could
  serve as the UTR null there; they understate UTR variance at long |u|.

## Conclusion (step 1)
- With u-intrinsic features only (coding score, Kozak, |u|, any transform), a fitted classifier does
  not beat the release rule outside its training domain (bench -> sim AUC 0.92-0.94 vs 0.93), and in
  domain gains ~+0.02 AUC, mostly from GC(u). The rule is near the ceiling of what u itself says.
- The usable extra separation comes from transcript context (transcript length, 3' UTR length, whether
  the ORF is internal, ORF length): fitted in domain it raises accuracy on the simulated assemblies
  from 0.87 to 0.93-0.94 (human d1), 0.82 to 0.91-0.92 (d32) and 0.76 to 0.91 (athal), recovering
  about half of the exact-ORF gap (+533 to +1,141 exact ORFs per assembly). Context relations differ
  between inputs (benchmark-trained context models fail on d32 and athal), so they must be learned
  from each input itself: that is what step 2 would test, and its main risk is that synthetic cuts
  reproduce the context distributions of real assemblies poorly.

## Caveats
- Labels come from simulation truth; only stop-matched candidates are studied (decisions for wrong
  stops are irrelevant to exact accuracy).
- cv_sim results are a ceiling for an input-specific model with true labels; a self-trained model
  (step 2) cannot see true labels.
- Benchmark context features are confounded with input type by construction (cut model, UTR clipping).
- Decoy segment sampling uses one shuffle per transcript and one random stretch per candidate.


## Step 2: self-trained decisions (no truth labels), 2026-10-28

### Method (scripts `pipeline/start_classifier/04_selftrain.py`, `04_selftrain_all_v1.sh`, `05_evaluate.py`, `06_summarize.py`)
Per input, using only the finished release run (hexamer model, mu, Kozak model, selected ORFs) and the
features of its real 5'-open candidates; truth is read only by `05_evaluate.py`.
- **(a) Synthetic self-training.** Confident ORFs = selected ORFs (before refinement) beginning with ATG
  and >= 100 codons (complete or 3' partial; the benchmark `fragments` inputs contain no complete
  CDSs), at most 4,000 per input (seeded subset). Per ORF 6 cuts inside the CDS (truncated examples, u
  = real CDS) and 6 cuts in the open 5' UTR after the last upstream in-frame stop (UTR examples, u =
  real UTR, the ATG is the true start); with probability = the share of the input's real candidates
  lacking a stop, the 3' end is also cut inside the CDS. Remaining ORF >= 5th percentile of the real
  candidates' ORF length. Every synthetic candidate is scored with PyTD's code exactly as a real one.
  Decoy UTR fallback below 300 UTR examples was never needed (0 decoys in all 55 inputs; 15,605-23,400
  examples per class except spombe/pfalciparum fragments and mixed, 7,347-13,350).
  Weights: `nat` (classes balanced) or `match` (each class reweighted to the real candidates' joint
  histogram of |u| decile x transcript-length quintile, clipped at 20).
  Models: logistic regression on `u` (|u|, log|u|, C/mu, C/mu/|u|, C/mu/sqrt|u|, K, Kozak bases off,
  GC(u)), `uc` (+ log 3' UTR, lacks stop, log ORF length), `uct` (+ log transcript length); boosted trees
  on `uct` for comparison. Each with Saerens EM prior adjustment (default) and without (`_noem`, prior 0.5).
- **(b) Unsupervised mixture on the real candidates.** Intrinsic log-odds = release rule score
  calibrated on the synthetic examples (`b_rule_*`) or the (a) `u` model (`b_au_*`); EM over pi and
  naive-Bayes histograms of context (3' state = lacks stop or 3' UTR quintile; ORF-length quintile;
  `_tx` + transcript-length quintile), initialised from the intrinsic posteriors. `b_rule_pi` = EM over
  pi only (prior-adapted calibrated rule).
- **Evaluation** (`05_evaluate.py`): each method's predictions are re-derived from the run
  (`preds.tsv` = refined `best_candidates.gff3`; verified identical, e.g. 42,582/42,582 rows on sim d1)
  by applying its keep/trim decision to every 5'-open candidate; stops never change. Scored as
  `score_predictions.py` / `13_depth_collect.py` (2 nt; excluded contigs removed on the simulated
  assemblies; per-type F1 = 2 tp_exact / (n_truth + predictions on that type)). Check: `release`
  reproduces the official exact F1 (e.g. spombe fragments 0.8295 = `utr_prior_eval.bench.tsv`) with
  0 changed predictions. `perfect` = oracle decisions (ceiling).

### Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
# manifest data/start_classifier/step2_manifest.tsv (55 inputs: 48 benchmark + human d1-d32 + athal d1)
nohup bash $T/pipeline/start_classifier/04_selftrain_all_v1.sh &   # 08:37-08:57, serial, OMP_NUM_THREADS=3; log data/start_classifier/selftrain/run.log (55/55, ALLDONE, no errors)
$P/env/startclf/bin/python $T/pipeline/start_classifier/06_summarize.py $T/data/start_classifier/step2 > $T/data/start_classifier/step2.summary.stdout.txt
```
Outputs: `data/start_classifier/selftrain/<dataset>.{decisions,synthetic}.tsv.gz, .log.tsv`;
`data/start_classifier/step2/{decision_metrics,exact_f1,selftrain_log,summary_by_group,summary_sim}.tsv`.

### Results (change in exact F1, points, vs release; `summary_by_group.tsv`, `summary_sim.tsv`, `exact_f1.tsv`)
Release exact F1: sim human d1-d32 0.806-0.846, athal d1 0.882; perfect start decisions would add
+3.6 (d1) to +11.2 (d32) and +8.7 (athal); benchmark complete +0.7..+6.7.

| variant | sim human d1/d2/d4/d8/d16/d32 | sim athal d1 | bench complete (16 inputs) | bench fragments (8) | bench mixed (24) |
|---|---|---|---|---|---|
| a_u_match_noem | +0.83/+0.86/+0.90/+1.15/+0.83/+1.58 | +0.56 | mean +0.24; spombe -3.0/-3.5, dmel -0.26/-0.23 | +1.17 (1 down, -0.48) | +1.41 (0 down) |
| a_hgb_uct | +0.07/+0.41/+0.33/+1.02/+2.26/+3.63 | -0.05 | -0.14 / +2.44 (4 of 8 complete down, to -5.3) | +1.33 | +1.83 (3 down) |
| a_u_match (EM) | mean -1.16 (4 of 6 down, to -4.4) | -0.15 | +0.6..+6.7 (0 down) | +2.95 | +3.50 (2 down) |
| b_rule_pi | -2.09/-2.56/-3.56/-3.33/-1.44/+1.16 | -0.55 | +0.6..+6.7 (0 down; = perfect within 0.1 except creinhardtii 0.3, tthermophila 0.4) | +1.81 | +1.87 (3 down) |
| b_rule_ctx | +0.76/+0.73/+0.38/+0.69/+0.85/+2.69 | +4.38 | +0.6..+6.7 (0 down) | +8.91 (0 down) * | +1.89 (1 down, -0.016) |
| b_au_ctx | +0.44/+0.33/+0.48/+1.32/+1.67/+3.19 | +3.66 | +1.94 / +3.11 (2 complete down, to -0.9) | +7.92 * | +0.55 (10 down) |

\* Benchmark-artefact: in `fragments` inputs every 5'-open UTR-extension candidate is a 3' partial (no
stop) and every truncated one ends in a stop, so "ORF lacks a stop" separates the classes perfectly
(context AUC 0.98-1.00); the mixture finds this. These fragment gains say nothing about real data.

Per type (simulated assemblies), the overall gains of the mixture come from complete and internal
CDSs while 5' and 3' partials lose: b_rule_ctx 5'-partial F1 -14.1 to -31.7 (human), -33.9 (athal),
3'-partial -9.4 to -20.9 / -16.3, complete +2.1..+2.9 / +5.8, internal +9.0..+12.9 / +17.9. The mixture
learns "lacks stop -> truncated" and "has stop -> UTR", which is right for internal and complete CDSs
and wrong for 3' and 5' partials. a_u_match_noem is near neutral per type (5' partial -2.3..+0.4 human,
+1.8 athal; complete +0.9..+1.3; internal -1.5..+2.7).

Decision accuracy (keep/trim, labelled candidates; release in parentheses): a_u_match_noem human sims
0.847-0.903 (0.822-0.868), athal 0.778 (0.763); AUC 0.936-0.947 (0.911-0.930), athal 0.843 (0.817).
Step-1 in-domain oracle ceiling (hgb_full, cv_sim) was 0.918-0.942 / 0.914.

Prior estimates (`selftrain_log.tsv`), true truncated share in parentheses: b_rule_ctx human d1-d32
0.25/0.28/0.33/0.46/0.60/0.70 (0.27/0.29/0.37/0.50/0.65/0.72), athal 0.27 (0.33); complete inputs
0.000-0.034 (0). Saerens on the synthetic-trained models is biased: unmatched `a_u` overestimates
(human d1 0.38, athal 0.53, complete inputs up to 0.29), matched `a_u_match` underestimates on the
assemblies (d1 0.09, athal 0.23) but finds ~0 on most complete inputs (spombe 0.17-0.20). The
synthetic examples are not calibrated for the real candidates, so prior adaptation through them is
unreliable; it is the main reason the EM variants lose on the assemblies.

Does |u|/transcript-length matching matter? Yes, mainly through the prior: with prior 0.5 (`_noem`),
matching turns the u-only model from a loss (`a_u_noem`: sims -1.3, complete -2.5/-1.2) into a gain
(`a_u_match_noem`: sims +0.6..+1.6, complete +0.2/+0.2); with EM it flips the bias of the pi estimate
(see above). Adding context features to the synthetic models (`uc`, `uct`) did not help on the
assemblies (`a_uc_match_noem` +0.27..+1.25 human, -0.64 athal) because synthetic 5' cuts do not
reproduce how real assembly truncation co-varies with 3' completeness and transcript length.

### Conclusion (step 2)
**No variant meets the success criteria** (5' partials better without loss on complete, and overall
exact F1 >= release on every input type):
- `b_rule_ctx` (calibrated rule + EM mixture over context) has the best overall exact F1 and is
  >= release on 55/55 inputs within 0.02 points (human mixed t1 -0.016): +0.4..+2.7 human assemblies,
  +4.4 athal, +0.6..+6.7 complete benchmark (= the oracle's gain within 0.1 except creinhardtii and
  tthermophila, within 0.3), +1.9 mixed. But it wins by trimming more: 5'-partial F1 falls 14-34 points on every
  simulated assembly, and 3' partials 9-21 points. It fails the 5'-partial criterion, and its fragment
  gains are the stop-flag artefact.
- `a_u_match_noem` (u-only logistic regression on matched synthetic cuts, no prior adaptation) is the
  only variant that is per-type near neutral and gains on every assembly (+0.6..+1.6), but it loses on
  spombe (-3.0/-3.5) and slightly on dmel (-0.25) complete inputs.
- Prior adaptation alone (`b_rule_pi`) recovers the oracle gain on complete inputs (pi -> 0) but
  underestimates pi on assemblies (-1.4..-3.6 points).

**Recommendation:** do not change the release rule. What self-training can learn honestly is (i) the
input's truncated share (works well when it is ~0, i.e. complete inputs, and from the mixture
`b_rule_ctx` within ~0.05 on the assemblies), and (ii) a u-only reweighting worth ~+1 point on
assemblies. Context features learned without labels help overall F1 only by trading 5'/3' partial
accuracy for complete/internal accuracy. If pursued, the candidate to test next (not tested here, it
would be post hoc) is the u-only matched model combined with the mixture's pi estimate (the better
calibrated prior), evaluated on the assemblies and GM12878.

### Caveats (step 2)
- Only 7 simulated assemblies (6 human depth levels share genes and reads) and one athal assembly are
  honest tests; benchmark fragments and mixed inputs are artefact-prone for any context feature.
- Synthetic training ORFs are selected ORFs of the same run (the hexamer model was trained on the
  longest of them), so their coding scores are slightly optimistic.
- One seed per input; differences below ~0.3 points on a single input were not tested for stability.

## Storage note (2026-10-28)
The large intermediate outputs were moved out of the repository after the analysis: `data/start_classifier/features/` (33 MB) and `data/start_classifier/selftrain/` (144 MB) now live at `$P/data/start_classifier/features/` and `$P/data/start_classifier/selftrain/`. Paths in the commands above refer to the original location in `$T`. The summary tables (`step1/`, `step2/`, `*.stdout.txt`) remain in the repository.

## Outcome
Not adopted; Brian decided (2026-10-28) that this stays out of the manuscript and is recorded in `notes/EXPERIMENT_LOG.md`. Step 3 (homology override, below) was added afterwards at Brian's request; its outcome is open (a method change needs his decision).

## Step 3: homology (BLAST / Pfam) evidence for keep-open vs trim (2026-10-28)

### Question
PyTD uses homology only for ORF selection, never in start refinement. Can a DIAMOND/Swiss-Prot or
Pfam hit of the 5'-open ORF (the longest-ORF peptide, before trimming) decide keep-open vs trim?

### Data and searches
- Benchmark: the existing `pytd_v7_homology` runs (`$P/work/bench_v3/runs/<org>/<input>/pytd_v7_homology/`)
  for complete, complete_fullutr and fragments of the 8 organisms (24 inputs). Their `longest_orfs.gff3`
  is byte-identical to the `pytd_v7` runs' (checked with `cmp` for all 24), so ORF ids join directly to
  the step-1 candidates (release `pytd_v7` decisions). There are no homology runs for the mixed inputs;
  they were skipped.
- Simulated assemblies human d1, d8, d32 and athal d1: new searches
  (`pipeline/start_classifier/07_homology_candidates_v2.sh`, driver `07_homology_all_v2.sh`;
  outputs `$P/work/start_classifier/homology/<dataset>/`; log `$P/work/start_classifier/homology_v2.log`,
  09:27-10:25, one dataset at a time). The parameters are those of `run_tool.sh` `pytd_v7_homology`:
  `diamond blastp --max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads 3` against
  `$P/data/homology/sprot_minus_<taxid>.dmnd` (human 9606, athal 3702), and `hmmsearch --cpu 1 --cut_ga
  --noali` against `$P/data/homology/Pfam-A.hmm` in 3 chunks run in parallel. DIAMOND 2.2.8 and
  HMMER 3.4 from `$P/env/td`.
  - Queries are the longest-ORF peptides of the 5'-open candidates only (13,205 / 11,230 / 9,710 /
    7,898), not all longest ORFs. Per-query results do not depend on the other queries (`-k 1`; GA
    thresholds are bit scores). Check: the same subset search on bench.human.fragments reproduced
    the full `pytd_v7_homology` run exactly for those ORFs: 10,682 / 10,682 best hits identical
    (subject, identity, qstart..send) and 19,267 / 19,267 Pfam domains identical.
  - v1 of the driver failed at once because the feature tables had been moved to
    `$P/data/start_classifier/features/`. Its empty output was moved to
    `$P/_trash/e33_homology_v1_failed_20261028/`; v2 differs only in that path.

### Method (`pipeline/start_classifier/08_homology_rules.py`, then `05_evaluate.py`)
- For each candidate:
  - a = peptide position of the first in-frame ATG (seg/3 + 1);
  - BLAST: the first (best) line of the query, giving qstart and sstart;
  - Pfam: the smallest envelope start (env_from) over the ORF's domains.
- Rules:

  | rule | call | fires when |
  |---|---|---|
  | R1_k[_s] | keep open | qstart <= a - k (k = 1, 5, 10, 20); `_s` also requires sstart > 1 |
  | R2_s_d | trim | sstart <= s and \|qstart - (sstart-1) - a\| <= d (s = 1, 3, 5, 10; d = 2, 5, 10) |
  | R3_k | keep open | a Pfam envelope starts at <= a - k (k = 1, 5) |

- Policies: the homology call where the rule fires, otherwise the release decision.
  - Each rule alone.
  - The combination chosen before evaluation: `H` = R2_3_5 (trim) + R1_5_s (keep) + R3_1 (keep).
  - `H_noR3` = the same without R3.
  - Conflicting calls fall back to the release decision.
  - `oracle_on_H` = the true decision wherever H fires, i.e. the ceiling for H's coverage.
- Scoring as in step 2.

### Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; cd $T
nohup bash pipeline/start_classifier/07_homology_all_v2.sh > $P/work/start_classifier/homology_v2.log 2>&1 &
OMP_NUM_THREADS=3 $P/env/startclf/bin/python pipeline/start_classifier/08_homology_rules.py $P/data/start_classifier/features data/start_classifier/step3_manifest.tsv $P/data/start_classifier/step3
$P/env/startclf/bin/python pipeline/start_classifier/05_evaluate.py $P/data/start_classifier/features $P/data/start_classifier/step3/decisions data/start_classifier/step3_manifest.tsv data/start_classifier/step3
cp $P/data/start_classifier/step3/rule_metrics.tsv data/start_classifier/step3/
```
Outputs:
- `data/start_classifier/step3/{rule_metrics,exact_f1,decision_metrics}.tsv`;
  `selftrain_log.tsv` there is only a placeholder.
- Decision columns are in `$P/data/start_classifier/step3/decisions/`.

### Results
**Coverage** (`rule_metrics.tsv`, rows `(coverage)`). A BLAST hit exists for:
- human assemblies: 92-96% of truncated and 99% of UTR-extension candidates;
- athal d1: 48% / 65%.

A Pfam hit exists for:
- human: 64-71% / 90-93%;
- athal: 73% / 89%.

**Rules on the simulated assemblies** (precision of the call on labelled candidates where the rule fires; recall of its class):

| rule | human d1 | human d8 | human d32 | athal d1 |
|---|---|---|---|---|
| R1_5_s keep: prec / recall | 0.909 / 0.795 | 0.972 / 0.824 | 0.990 / 0.858 | 0.961 / 0.339 |
| R1_5_s: changed release decisions, exact fixed / broken | 638, 323 / 113 | 969, 719 / 56 | 1141, 956 / 30 | 262, 226 / 22 |
| R2_3_5 trim: prec / recall | 0.997 / 0.776 | 0.994 / 0.821 | 0.994 / 0.866 | 0.993 / 0.132 |
| R2_3_5: changed, fixed / broken | 558, 438 / 4 | 417, 369 / 7 | 266, 226 / 7 | 107, 104 / 2 |
| R3_1 keep: prec / recall | 0.863 / 0.519 | 0.938 / 0.474 | 0.978 / 0.485 | 0.930 / 0.493 |

- Larger k raises R1's precision and lowers its recall: R1_20_s reaches 0.93-0.99.
- R2's precision is 0.99-1.00 for every s and d.
- On the benchmark complete inputs every 5'-open end is a UTR extension:
  - R2_3_5 trims 23,063 candidates (16 inputs), fixing 2,508 exact ORFs and breaking none;
  - the keep rules fire wrongly by construction: R1_5_s 682 times (415 exact broken), R3_1 1,597 (1,026 broken).
  - These are candidates where the homolog's alignment begins upstream of the annotated ATG: an
    N-terminal extension in the homolog or a later annotated start.

**Exact F1, change in points vs release** (`exact_f1.tsv`):

| policy | human d1 | human d8 | human d32 | athal d1 | bench complete (16) | bench fragments (8) * |
|---|---|---|---|---|---|---|
| H | +2.02 | +5.44 | +9.00 | +2.02 | mean +0.53 (-0.57..+2.94; 8/16 down) | +4.30 (+0.31..+8.53) |
| H_noR3 | +1.93 | +5.25 | +8.70 | +1.49 | mean +0.75 (-0.15..+3.08; 3/16 down) | +3.18 |
| R2_3_5 alone | +1.30 | +1.86 | +1.66 | +0.50 | mean +0.89 (+0.02..+3.41; 0/16 down) | +0.43 |
| R1_5_s alone | +0.63 | +3.40 | +7.04 | +0.99 | mean -0.14 (14/16 down, to -0.32) | +2.75 |
| oracle_on_H | +2.82 | +6.31 | +9.62 | +2.54 | +1.06 | +4.56 |
| perfect (all candidates) | +3.62 | +7.61 | +11.22 | +8.66 | +4.05 | +9.27 |

Release exact F1 for comparison: d1 0.844, d8 0.832, d32 0.806, athal 0.882.

\* The homology rules do not use the "ORF lacks a stop" flag, so the step-2 benchmark artefact does
not apply to them. The benchmark fragments are still synthetic cuts.

- Per type on the assemblies, H improves every type:
  - complete +1.0 / +1.9 / +2.9 / +0.2;
  - 5' partial +10.6 / +13.1 / +12.7 / +14.2;
  - 3' partial +1.2 / +3.7 / +5.1 / +0.7;
  - internal +9.7 / +12.9 / +12.9 / +12.3 (human d1 / d8 / d32 / athal).
- H_noR3 also improves every type (5' partial +8.0 to +12.1).
- H recovers 56% (d1), 72% (d8), 80% (d32) and 23% (athal) of the gap to perfect decisions, and
  71-94% of what perfect decisions on its own firings would give (`oracle_on_H`).
- Athal recovers less because BLAST coverage is lower (Swiss-Prot minus Arabidopsis is sparser for
  plants). R3 matters more there (H +2.02 vs H_noR3 +1.49).
- The complete inputs that lose under H are athal, creinhardtii, dmel and tthermophila, at most
  -0.57 (tthermophila): R3 keeps the end open when a Pfam envelope reaches upstream of the annotated
  start. Under H_noR3 the losses are creinhardtii -0.15 / -0.09 and tthermophila -0.003.

### Conclusion (step 3)
Homology evidence decides keep-open vs trim far better than any self-trained classifier (step 2), and
it improves every CDS type on the simulated assemblies.
- **Trimming when the subject's N-terminus aligns to the first ATG (R2) is almost never wrong**
  (precision >= 0.99). It never lowers exact F1 on any input (0 of 16 complete inputs down, +0.02 to
  +3.41).
- **Keeping the end open when the alignment starts upstream of the ATG (R1, R3)** gives most of the
  assembly gain: R1_5_s alone adds +0.6 to +7.0. It costs up to 0.3 points (R1) or 0.6 (R3) on
  complete inputs.
- **Policies:**
  - `H` gives +2.0 / +5.4 / +9.0 (human d1 / d8 / d32) and +2.0 (athal), with a mean +0.5 on complete
    inputs but 8 of 16 slightly down.
  - `H_noR3` gives +1.9 / +5.3 / +8.7 / +1.5, with a mean +0.75 on complete inputs and 3 of 16 down
    by at most 0.15.
- **A rule with no loss on any complete input** would need a stricter keep rule. The R1_20_s results
  above (precision 0.93-0.99; -0.03 on complete) suggest one, but no combination with it was tested,
  and choosing one now would be post hoc.

**Recommendation:** a homology override for start refinement is worth implementing as an option used
when homology results are supplied (`--retain-blastp-hits` / `--retain-pfam-hits` are already
inputs). Candidates are R2 (trim) alone as the zero-risk version, or R2 + R1 (k = 5, sstart > 1) for
the larger assembly gain. It needs Brian's decision and fresh runs. Before adoption:
- validate on GM12878 with MANE truth;
- validate with the full Swiss-Prot, i.e. without excluding the organism's own entries, as users
  would run it;
- test a stricter R1 (k = 20) on complete inputs.

### Caveats (step 3)
- Swiss-Prot minus the organism's own entries (as in the benchmark homology runs). With the own
  species included, coverage and agreement would be higher, especially for human.
- Only the best DIAMOND hit (`-k 1`) and its first HSP are used, as in the TransDecoder wrapper.
- The rule parameter grids are reported in full, but the combined policy was fixed before
  evaluation. Simulated human assemblies at different depths share genes, so they are not
  independent tests.
- Benchmark mixed inputs were not evaluated (no homology runs).
