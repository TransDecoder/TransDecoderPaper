# E33. Keep-open vs trim for 5'-open ORFs: how much better can a classifier decide? (2026-10-28)

Status: steps 1-2 exploratory, not adopted. Step 3 (homology evidence) ADOPTED in PyTransDecoder only
(homology mode, a5b24c3, 2026-10-28; see "Adoption in PyTD"). Perl TransDecoder unchanged.
Step 5 (UniRef50 instead of Swiss-Prot, 2026-10-29): not adopted; Swiss-Prot stays the recommended database.

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
Not adopted; Brian decided (2026-10-28) that this stays out of the manuscript and is recorded in `notes/EXPERIMENT_LOG.md`. Step 3 (homology override) was added afterwards at Brian's request and was ADOPTED in PyTransDecoder only (homology mode; commit a5b24c3, section "Adoption in PyTD").

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

## Adoption in PyTD: homology-assisted start refinement (2026-10-28)

### Decision
Brian (2026-10-28) chose to implement policy H from step 3 in PyTransDecoder only, not in Perl:
- the BLAST rules (R1, R2) apply only if `--retain-blastp-hits` is given;
- the Pfam rule (R3) applies only if `--retain-pfam-hits` is given;
- the feature is on by default whenever hits are supplied, with the opt-out `--no-homology-start-assist`.

### Code
**Commit:** PyTD `main` **a5b24c3**, "Homology-assisted start refinement when BLASTP/Pfam hits are supplied (v0.2.0)".
- Parent 3a793ff, which was `origin/main` before the push.
- Pushed to `git@github.com:TransDecoder/PyTransDecoder.git` as `main`; not tagged.
- Worktree: `$P/work/trees/PY_hstart`.

**Rules** (`pytransdecoder/core/homology_starts.py`), with a = peptide position of the first in-frame ATG:

| rule | call | fires when |
|---|---|---|
| R2 | trim | sstart <= 3 and \|qstart - (sstart-1) - a\| <= 5 |
| R1 | keep open | qstart <= a - 5 and sstart > 1 |
| R3 | keep open | Pfam env_from <= a - 1 |

- A trim call together with a keep call (a conflict), or no call at all, falls back to the ab initio rule.
- R1 and R2 cannot fire together on one hit, so conflicts are always R2 + R3.

**Changes in `predict.py`:**
- The hit parsers are refactored into `_read_blastp_first_hits` and `_read_pfam_domains`. The existing selection parsers `_parse_blastp_hits` and `_parse_pfam_hits` are built on them, with unchanged output.
- `_homology_start_evidence` reads the best hit's qstart and sstart, and each ORF's smallest envelope start.
- `_revise_start_sites` applies the rules. It writes `start_refinement.homology.tsv` (orf id, ATG offset in nt and aa, BLAST qstart and sstart, Pfam env_from, rules fired, coding and UTR scores, ab initio decision, final decision, source), logs the counts of R1/R2/R3 firings, homology trims and keeps, conflicts and changed decisions, and tags R2 trims with `start_homology=R2`.

**CLI and documentation:**
- `--no-homology-start-assist` added to `predict` and to the `pyTransdecoder` pipeline.
- README (start refinement and options), `docs/start_codon_refinement.md` (new section and comparison-table row), new `CHANGELOG.md` (v0.2.0); version set to 0.2.0.

### Tests
- **Unit tests:** `tests/test_homology_starts.py`, 9 tests: each rule and its boundaries, the conflict fallback, Pfam-only and opt-out behaviour, and the per-ORF table. All 76 PyTD tests pass at a5b24c3.
- **Ab initio regression** (`pipeline/start_classifier/09_pytd_hstart_abinitio_regression_v1.sh`; output `$P/work/start_classifier/hstart_regression/`, summary `data/start_classifier/adoption/abinitio_regression.txt`):
  - Run on the working tree that was then committed unchanged as a5b24c3.
  - Five sample datasets (trinity, pasa, stringtie, supertx, cufflinks, as in `regression/run_regression.sh`): PY_v7 (3a793ff) and PY_hstart run fresh.
  - Four benchmark inputs (spombe complete and fragments, dmel complete, pfalciparum fragments_mixed_t1): PY_hstart run fresh, compared with the existing 3a793ff `pytd_v7` runs.
  - Files compared: final gff3/bed/cds/pep, `best_candidates.gff3` and its `revised_starts`. All 54 files are byte-identical.
- **In-program vs offline decisions** (`12_check_inprogram_vs_offline.py`; `adoption/inprogram_vs_offline_{bench,sim}.tsv`):
  - Scope: on the 24 benchmark homology inputs (complete, complete_fullutr and fragments of 8 organisms; the new `pytd_v7h_homology` runs), every 5'-open ORF that both the program and the offline step-3 policy H evaluate.
  - Benchmark: **150,356 ORFs compared, 0 decision mismatches**, and 0 ab initio-decision mismatches. Homology selection adds some ORFs that the ab initio selection lacks; those are not compared.
  - Simulated assemblies: 42,043 compared, 0 mismatches.
  - Rule counts in the simulated runs' logs equal the step-3 offline counts. Human d1, for example: R1 3,199, R2 6,793, R3 2,308, 1,317 changed decisions.

### Evaluation runs (fresh, with the pushed code a5b24c3)
- **Driver:** `pipeline/start_classifier/15_eval_jobs_v1.sh` runs `11_run_pytd_v7h_homology_v1.sh` and `13_gm12878_trinity_v7h_v1.sh`, 3 at a time, 11:19-14:36. Log `$P/work/start_classifier/eval_jobs_v1.log`: 92 jobs, 0 failed.
- **Each run:** a new run directory, a fresh LongOrfs (`-S` for stranded inputs, the organism's genetic code), then Predict with `--retain-blastp-hits` / `--retain-pfam-hits`.
- **Hits are reused, as a user would supply them:**
  - benchmark: the existing release homology run `pytd_v7_homology` (3a793ff; DIAMOND to Swiss-Prot minus own organism, Pfam-A, on all LongOrfs peptides). Its LongOrfs is identical to the fresh one.
  - simulated assemblies: the step-3 candidate searches.
  - GM12878 Trinity: the release run `trinity/final_v7_homology`. Its `longest_orfs.pep` was checked identical (`cmp`).
- **Not run:** GM12878 LRAA and StringTie have no homology searches on the v7 LongOrfs. Their homology runs predate v7 and have different LongOrfs, so they were skipped.
- **Scoring:** `16_score_adoption.py` (as `score_predictions.py` / E25; per type as `02_partial_types.py`) writes `adoption/{bench_scores,bench_summary,sim_scores}.tsv`. GM12878 uses `eval_genome_cds.py` (as E03/E25; `adoption/gm12878_trinity_summary.tsv`) and the E29 statistic (`14_gm12878_e29_stat.py`; `adoption/gm12878_trinity_e29_stat.tsv`).

### Results
**Benchmark, 84 homology inputs.** The new code is `pytd_v7h_homology`; the comparison is with `pytd_v7_homology`, i.e. the same hits used for selection only.
- Predictions and stop-level matches are identical on all 84 inputs; only start codons change.

Change in exact F1 (points; median, range, number of inputs down):

| inputs | vs pytd_v7_homology | vs ab initio pytd_v7 |
|---|---|---|
| all 84 | +1.25 (-0.58..+9.18; 24 down) | +1.65 (-0.44..+12.90; 13 down) |
| fragments* (34) | +5.92 (+0.29..+9.18; 0 down) | +6.52 (+0.53..+12.90; 0 down) |
| complete* (34) | +0.02 (-0.58..+2.92; 16 down) | +0.67 (-0.44..+4.21; 8 down) |
| complete_fullutr* (16) | -0.01 (-0.54..+2.85; 8 down) | +0.68 (-0.42..+4.78; 5 down) |

Change in exact F1 by truth type in the fragments inputs, vs `pytd_v7_homology`:
- 5' partial: median +11.4 points (+0.9..+15.7);
- 3' partial: median +0.02 (-0.51..+3.25).

Per organism, the complete inputs are unchanged across the stranded, unstranded and wRand variants:

| organism | complete | complete_fullutr |
|---|---|---|
| human | +2.9 | +2.85 |
| mouse | +1.5 | +1.45 |
| spombe | +0.7 | +0.7 |
| pfalciparum | +0.02 | +0.02 |
| athal | -0.08 | -0.04 |
| dmel | -0.07 | -0.09 |
| creinhardtii | -0.35 | -0.26 |
| tthermophila | -0.58 | -0.54 |

The per-organism changes are the step-3 offline values. The losses come from the keep rules on transcripts whose homolog's alignment starts upstream of the annotated ATG.

**Simulated Trinity assemblies** (`sim_scores.tsv`; excluded contigs removed):
- Arms: ab initio `pytd_v7`; assist off (the same hits used for selection only); assist on.
- Assist off is identical to ab initio here, because the candidate hits do not change selection.
- Predictions and stops are identical in all three arms.

| assembly | exact F1 off -> on | complete | 5' partial | 3' partial | internal |
|---|---|---|---|---|---|
| human d1 | 0.8440 -> 0.8641 (+2.0) | 0.887 -> 0.897 | 0.663 -> 0.769 | 0.823 -> 0.835 | 0.765 -> 0.862 |
| human d8 | 0.8315 -> 0.8859 (+5.4) | 0.901 -> 0.920 | 0.706 -> 0.837 | 0.862 -> 0.899 | 0.754 -> 0.883 |
| human d32 | 0.8063 -> 0.8963 (+9.0) | 0.900 -> 0.929 | 0.745 -> 0.871 | 0.849 -> 0.900 | 0.773 -> 0.902 |
| athal d1 | 0.8815 -> 0.9017 (+2.0) | 0.919 -> 0.920 | 0.622 -> 0.764 | 0.784 -> 0.791 | 0.738 -> 0.861 |

**GM12878 Trinity GG** (MANE / GENCODE evaluation as E03/E25, `gm12878_trinity_summary.tsv`):

| metric | ab initio v7 | homology v7 (3a793ff) | homology v7h (a5b24c3) |
|---|---|---|---|
| MANE exact | 6,857 | 7,063 | 7,190 |
| MANE exact given assembled | 0.8933 | 0.9201 | 0.9367 |
| MANE compatible | 10,261 | 10,960 | 10,941 |
| GENCODE complete exact | 51,807 | 54,537 | 54,889 |
| GENCODE complete exact given assembled | 0.8563 | 0.9014 | 0.9072 |
| predicted CDS models | 77,481 | 96,969 | 96,421 |
| distinct predicted complete | 44,494 | 54,116 | 51,351 |
| distinct predicted 5' partial | 5,931 | 9,206 | 11,717 |

E29 statistic: transcripts carrying a MANE stop with the reference structure to their 5' end; `gm12878_trinity_e29_stat.tsv`.

5'-partial references (1,915):
- stop predicted: 1,188 (ab initio) and 1,471 (homology);
- left open (exact) given the stop: 0.791 (ab initio), 0.777 (homology v7), **0.965 (v7h)**;
- trimmed to an internal ATG: 248 / 328 / **52**.

Complete references (12,172):
- exact given the stop: 0.927 / 0.928 / **0.944**;
- 5' end upstream of the true start: 818 / 829 / 642.

### Conclusion
- The homology start assist is in PyTD main (a5b24c3, v0.2.0). Ab initio output is unchanged.
- With homology hits:
  - the simulated Trinity assemblies gain 2-9 exact-F1 points, every CDS type improving;
  - benchmark fragments gain +5.9 points (median);
  - GM12878 Trinity gains +127 MANE exact CDSs, and its 5' partials kept open given the stop rise from 78% to 96%.
- Cost on complete benchmark inputs: up to -0.58 points (tthermophila), with a median of about 0.
- Perl and Python now agree only on ab initio runs. Homology runs differ in start codons only; selection and stops are identical.

### Caveats (adoption)
- GM12878 LRAA and StringTie are not evaluated: there are no homology searches on the v7 LongOrfs.
- The simulated-assembly searches cover only the 5'-open candidate ORFs, so selection there is unaffected by homology. The comparison isolates the start assist.
- Swiss-Prot excludes the organism's own entries, as in E25. With the species included, rule coverage would rise.

## Arabidopsis depth series (E32) with the homology start assist (2026-10-28)

### Question
Brian asked how the homology start assist performs on the simulated Arabidopsis Trinity depth series (E32,
levels d1-d32), run as a user would, and how it separates from homology-based retention.

### Data and runs
- **E32 inputs (read only):** the assembly and labels of each level, `$P/data/simtrinity_athal/depth/<level>/{inputs,labels}`.
  - Levels d1-d16 had finished.
  - d32 was started only after E32's own d32 label/run/score had finished: `d32: LEVEL_DONE` in its `d32.run.log` and `scores/scores_overall.tsv` present, at 18:14. My lane started d32 at 18:15.
- **Code:** PyTD a5b24c3 (`$P/work/trees/PY_hstart`).
- **Scripts:** `pipeline/start_classifier/17_athal_depth_level_v1.sh` (one level) and `18_athal_depth_lane_v1.sh` (the lane, 15:54-18:28, 0 failures; log `$P/work/start_classifier/athal_depth/lane_v1.log`).
- **One LongOrfs per level:** a fresh `-S -G Universal` run. Its work directory is copied into each condition before a fresh Predict, so all three conditions start from the same LongOrfs.
- **Conditions:**
  - A: ab initio;
  - B: full homology mode, start assist on (the default);
  - C: the same hits as B with `--no-homology-start-assist`, i.e. homology for retention only.
- **Searches for B and C** run on all LongOrfs peptides with `run_tool.sh pytd_v7_homology`'s commands:
  - DIAMOND 2.2.8 `--max-target-seqs 1 --evalue 1e-5` against `sprot_minus_3702.dmnd`;
  - hmmsearch 3.4 `--cut_ga --noali` against Pfam-A, in 3 chunks;
  - 3 threads; outputs in `$P/work/start_classifier/athal_depth/searches/<level>/`.
- **The step-3 d1 searches** covered only the 5'-open candidates (step 3, "Adoption"). They were not reused for B and C here; d1 was rerun fully like the other levels.
- **Scoring with E32's scorer and labels:**
  - each condition's run sits at `$P/work/start_classifier/athal_depth/<A|B|C>/depth/<level>/runs/pytd_v7/`;
  - scored with `pipeline/simtrinity/score_simtrinity.py`, then collected with E32's `13b_depth_collect.py --data-root $P/data/simtrinity_athal`;
  - `19_athal_depth_summary.py` runs the same collector over E32's own runs and writes `data/start_classifier/adoption/athal_depth_scores.tsv` (stdout `athal_depth_scores.stdout.txt`).
  - The collected E32 table is identical to E32's committed `data/simtrinity_athal/depth/scores_by_type.tsv`.
  - "best_other" is the best E32 method other than TransDecoder (pytd_v7, perl_td_v600) and longest ORF: overall by exact F1, and per type by exact F1 with negatives.

### Checks
- Condition A reproduces E32's TransDecoder (pytd_v7, 3a793ff) at every level:
  - `preds.tsv` byte-identical;
  - `longest_orf.preds.tsv` identical at d1;
  - overall and per-type score rows identical.
- B and C have identical predictions and stops; only start codons differ.

Start-assist rule counts in B (from the Predict logs):

| level | R1 | R2 | R3 | homology trims | homology keeps | conflicts | changed vs ab initio |
|---|---|---|---|---|---|---|---|
| d1 | 999 | 708 | 1,525 | 685 | 1,769 | 23 | 648 |
| d32 | 3,746 | 442 | 4,612 | 421 | 5,604 | 21 | 1,132 |

### Results (`athal_depth_scores.tsv`)
Exact F1 overall. Stop F1 is in parentheses where the column holds both.

| level | A ab initio | C homology, no assist | B homology + assist | best other method | longest ORF |
|---|---|---|---|---|---|
| d1 | 0.8815 (0.9719) | 0.8910 (0.9841) | **0.9138** (0.9841) | 0.8798 (borf) | 0.7308 |
| d2 | 0.8727 (0.9699) | 0.8822 (0.9841) | **0.9148** (0.9841) | 0.8432 (borf) | 0.7473 |
| d4 | 0.8589 (0.9731) | 0.8683 (0.9863) | **0.9091** (0.9863) | 0.7878 (borf) | 0.7733 |
| d8 | 0.8379 (0.9684) | 0.8480 (0.9834) | **0.9008** (0.9834) | 0.7794 (prodigal) | 0.8152 |
| d16 | 0.8218 (0.9693) | 0.8321 (0.9845) | **0.8989** (0.9845) | 0.8235 (prodigal) | 0.8515 |
| d32 | 0.8281 (0.9710) | 0.8364 (0.9862) | **0.9083** (0.9862) | 0.8540 (prodigal) | 0.8760 |

Per-type exact F1 with negatives, B / C / A:

| level | complete | 5' partial | 3' partial | internal |
|---|---|---|---|---|
| d1 | 0.928 / 0.926 / 0.918 | 0.787 / 0.634 / 0.617 | 0.808 / 0.801 / 0.774 | 0.893 / 0.750 / 0.736 |
| d8 | 0.929 / 0.926 / 0.916 | 0.836 / 0.706 / 0.688 | 0.868 / 0.870 / 0.852 | 0.893 / 0.791 / 0.784 |
| d32 | 0.929 / 0.919 / 0.907 | 0.869 / 0.756 / 0.743 | 0.872 / 0.865 / 0.847 | 0.914 / 0.819 / 0.815 |

Other levels are in the table file. The best other method per type is the same at every level: RNAsamba for complete and 3' partial, Prodigal for 5' partial and internal.

| level | RNAsamba complete | Prodigal 5' partial | RNAsamba 3' partial | Prodigal internal |
|---|---|---|---|---|
| d1 | 0.973 | 0.893 | 0.858 | 0.975 |
| d32 | 0.988 | 0.950 | 0.944 | 0.987 |

### Summary
- Homology for retention (C vs A) adds +0.8 to +1.0 exact-F1 points at every level, and +1.2 to +1.5 points stop F1.
- The start assist (B vs C) adds a further +2.3 (d1), +3.3, +4.1, +5.3, +6.7 and +7.2 (d32) points.
- So B's gain over ab initio grows from +3.2 (d1) to +8.0 (d32) points. B's exact F1 stays at 0.899-0.915 across depths, while A falls from 0.882 to 0.822-0.828.
- Most of the gain is in 5' partials (+11 to +16 points vs C) and internal CDSs (+9 to +15). Complete and 3'-partial CDSs change by -0.2 to +0.9 points.
- B has the best overall exact F1 of all methods at every level: +3.4 (d1) to +12.1 (d8) points over the best non-TransDecoder method, and +3.2 to +18.3 over longest ORF.
  - Per type, the specialist methods still lead: RNAsamba on complete (0.97-0.99) and Prodigal on 5' partial and internal.
  - These per-type F1 values count only predictions on transcripts of that type, so a method that predicts few ORFs on other types is not penalised here. E32 discusses how to read them.

### Caveats
- Swiss-Prot without Arabidopsis entries (as E25/E32): BLAST coverage of plant CDSs is lower than for human (step 3: 48-65% of candidates). With the species' own proteins included, B would gain more.
- Simulated reads; one replicate per level.

## Step 4: hit depth and e-value (2026-10-28)

### Question
Brian asked which 5'-partial and 3'-partial CDSs the homology start assist (PyTD a5b24c3) still gets wrong.
He also asked whether using only the top DIAMOND hit (`-k 1`, first HSP) at e-value 1e-5 is too stringent
to collect useful evidence. The alternatives are the top n hits checked for consistency, looser e-values,
all HSPs, more-sensitive DIAMOND, or a Pfam E-value threshold instead of `--cut_ga`.

### Data and searches
- **Error anatomy:** the assist-on runs, read only.
  - Simulated human d1 / d8 / d32 and athal d1: `$P/work/start_classifier/hstart_runs/<ds>/on`, with the step-3 candidate searches. Homology there does not change selection.
  - Athal d1, full homology mode: condition B of the Arabidopsis depth series, `$P/work/start_classifier/athal_depth/B/depth/d1/runs/pytd_v7`. It is another agent's output and was only read. Its searches ran on all LongOrfs, so selection changes too.
  - Truth, excluded contigs and scoring are as in `16_score_adoption.py`. Per-type non-exact counts equal `adoption/sim_scores.tsv` (e.g. human d1 complete 3,692, 5' partial 486).
- **Deep searches** (`pipeline/start_classifier/step4_01_deep_search_v3.sh`; outputs `$P/data/start_classifier/step4/searches/<dataset>/`; log copied to `data/start_classifier/step4/deep_search_v3.log`; 16:36-19:19).
  - Queries: the longest-ORF peptides of the 5'-open candidates (`01_features.py` tables), i.e. the step-3 candidate peptides.
    - Simulated: human d1 / d8 / d32 and athal d1.
    - Benchmark: complete and fragments inputs of all 8 organisms. The task asked for at least human, athal, dmel and tthermophila; there are 2,596-15,222 candidates per input.
  - DIAMOND 2.2.8 against the same `sprot_minus_<taxid>.dmnd` with `--max-target-seqs 25 --max-hsps 0 --evalue 10`, default (fast) mode and `--more-sensitive`; 8 threads, `nice -n 19`, one search at a time.
  - hmmsearch 3.4 `--domE 1 -Z 1 --domZ 1 --noali` against Pfam-A.
    - With `-Z 1 --domZ 1` every E-value is a P-value. The results then do not depend on the query set or on chunking, and `--domE 1` reports every domain that passes the filter pipeline.
    - The search ran as 7 serial chunks (`--cpu 0`, 7 threads in all).
    - Step 4 applies "domE" thresholds as i-Evalue at a reference domZ of 10,000. In the step-3 `--cut_ga` searches the i-Evalue / P ratio was 4,267-4,546, the number of sequences per chunk.
    - So `dE1` = P <= 1e-4, `dE1e-1` = 1e-5, `dE1e-3` = 1e-7, `dE1e-5` = 1e-9. For comparison, 90% of GA domains have P <= 1.1e-8.
  - Two earlier versions were stopped by me:
    - v1: a single `hmmsearch --cpu 7` used under 2 CPUs, because the master thread re-reads the targets for every HMM;
    - v2: `--cpu 1` chunks are 2 threads each.
    - Their partial outputs are in `$P/_trash/e33_step4_hmmsearch_v1_killed/`. The DIAMOND results of v1 (sim.human.d1) were kept, because the commands are identical.
- **Check:** the deep fast search's top hit at 1e-5 against the a5b24c3 `-k 1` hit (`top1_vs_k1.tsv`).
  - The same queries have hits, plus 0-20 extra per input.
  - Same subject: 96-100%.
  - Same keep/trim call on 99.9-100% of candidates.

### Method
- **Error anatomy** (`step4_02_error_anatomy_v1.py`). Every truth CDS that is not exact (stop and start within 2 nt) gets one class:
  - `selection`: no prediction with the true stop. Split into "no prediction on the transcript" and "wrong ORF / stop".
  - `outside`: right stop, but the predicted ORF is not a 5'-open candidate with an in-frame ATG, or the true start is neither the open end nor the first ATG. Neither keep nor trim can be exact.
  - `no hit`: wrong keep/trim decision and the ORF has no DIAMOND hit and no Pfam domain.
  - `hit, no rule`: wrong decision, a hit exists but no rule fired, so the ab initio rule decided. The reason is recorded (a = peptide position of the first ATG):
    - BLAST alignment starts at sstart > 3 (bins 4-10, 11-50, > 50) without reaching >= 5 aa upstream of the ATG;
    - sstart <= 3 but the subject N-terminus maps > 5 aa from the ATG, or sstart = 1 with qstart <= a-5 (R1 needs sstart > 1);
    - Pfam only, with the envelope starting at or after the ATG.
  - `rule wrong`: a homology rule (R1, R2, R3) made the wrong call.
  - `conflict`: R2 + R3 conflict, and the ab initio fallback was wrong.
- **Offline variants** (`step4_03_rule_variants_v1.py`). The a5b24c3 rules R1 / R2 / R3, with R2 + R3 conflict and no call both falling back to the release decision, are applied with:
  - search: `k1` (the a5b24c3 inputs), `fast`, `msens`;
  - HSP e-value cutoff: 1e-10, 1e-5, 1e-3, 1e-1, 1, 10;
  - HSP per subject: `best` (highest bitscore) or `all` (the HSP reaching furthest to the subject N-terminus);
  - depth: `top1`, or the top n = 3, 5, 10, 25 subjects (by best-HSP bitscore) combined by one of five schemes:
    - `maj`: majority of the firing hits;
    - `unan`: all firing hits agree;
    - `bits`: larger summed bitscore of trim vs keep hits;
    - `anyup`: any R2 hit trims, otherwise any R1 hit keeps;
    - `first`: the best-ranked firing hit decides;
  - Pfam R3: `ga`, `dE1`, `dE1e-1`, `dE1e-3`, `dE1e-5`, or `none`.
  - This gives 3,030 variants per dataset.
  - Scoring: the release run's predictions with each candidate's decision substituted, as `05_evaluate.py`. Exact F1 = 2 tp / (n_truth + n_pred), overall and per type with negatives, as `16_score_adoption.py`.
  - The baseline `k1|1e-05|best|top1|ga` reproduces the adoption runs: human d1 / d8 / d32 0.8641 / 0.8859 / 0.8963, athal d1 0.9017. It also reproduces the step-3 `H` column.
- **Selection** (`step4_04_summarize_v1.py`):
  - tuning score = mean change in exact F1 over simulated human d1 and d8 only;
  - held out: human d32, athal d1, the 8 benchmark fragments inputs and the 8 benchmark complete inputs.
- **Evidence check** (`step4_05_evidence_v1.py`):
  - For labelled candidates, i.e. truncated (keep is right) and UTR extensions with the true start at the first ATG (trim is right), whether any of the 25 hits supports the right call or the wrong call. Every HSP is tested with R1 / R2; R3 uses Pfam P <= 1e-4.
  - The counts are split by whether a5b24c3 decides the candidate correctly.

### Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; B=$P/data/start_classifier/step4; S=$T/data/start_classifier/step4; cd $T
nohup bash pipeline/start_classifier/step4_01_deep_search_v3.sh > $B/deep_search_v3.log 2>&1 &
PY=$P/env/startclf/bin/python
$PY pipeline/start_classifier/step4_02_error_anatomy_v1.py $B/anatomy
$PY pipeline/start_classifier/step4_03_rule_variants_v1.py $S/step4_manifest.tsv $B/variants
$PY pipeline/start_classifier/step4_05_evidence_v1.py $S/step4_manifest.tsv $B/evidence.tsv
$PY pipeline/start_classifier/step4_04_summarize_v1.py $B/variants/variants_f1.tsv.gz $S
$PY pipeline/start_classifier/step4_06_tables_v1.py $B $S
```
Outputs:
- in `data/start_classifier/step4/`:
  - `anatomy_summary.tsv`, `anatomy_counts.tsv`, `anatomy_complete_left_open.tsv`;
  - `evidence_summary.tsv`;
  - `tuning_ranking.tsv` (all variants), `selected_variants.tsv`, `family_by_dataset.tsv`, `family_marginals.tsv`;
  - `top1_vs_k1.tsv`;
- large files in `$B`: `anatomy/anatomy_per_cds.tsv.gz`, `variants/variants_f1.tsv.gz` (every dataset x variant, overall and per type).

### Results: error anatomy (`anatomy_summary.tsv`)
CDSs that are not exact with the assist on, by class. "athal B" = depth-series condition B (full homology mode).
- sel: selection error, given as no prediction / wrong stop or ORF.
- outside: right stop but no keep/trim choice can be exact.
- hit, no rule: given as total (BLAST sstart > 3 / subject N-terminus mismatch / Pfam only downstream).

| type | dataset | not exact | sel: none / wrong | outside | no hit | hit, no rule (sstart>3 / N-term / Pfam-down) | rule wrong | conflict |
|---|---|---|---|---|---|---|---|---|
| 5' partial | human d1 | 486 | 391 / 21 | 0 | 12 | 50 (45 / 1 / 4) | 12 (R2) | 0 |
| | human d8 | 491 | 422 / 18 | 0 | 9 | 31 (28 / 0 / 3) | 10 (R2) | 1 |
| | human d32 | 376 | 320 / 13 | 0 | 3 | 33 (28 / 1 / 4) | 7 (R2) | 0 |
| | athal d1 | 299 | 124 / 4 | 0 | 39 | 128 (77 / 4 / 47) | 4 (R2) | 0 |
| | athal B d1 | 245 | 65 / 10 | 0 | 38 | 128 (77 / 4 / 47) | 4 (R2) | 0 |
| 3' partial | human d1 | 375 | 275 / 11 | 61 | 1 | 9 (6 / 1 / 2) | 16 (R1, R1+R3, R3) | 2 |
| | human d8 | 304 | 231 / 8 | 27 | 3 | 13 (12 / 1 / 0) | 20 | 2 |
| | human d32 | 242 | 196 / 2 | 15 | 1 | 9 (6 / 2 / 1) | 18 | 1 |
| | athal d1 | 157 | 84 / 7 | 5 | 19 | 33 (17 / 0 / 16) | 9 | 0 |
| | athal B d1 | 130 | 55 / 7 | 7 | 19 | 33 (17 / 0 / 16) | 9 | 0 |
| internal | human d1 | 609 | 544 / 0 | 0 | 25 | 34 (28 / 0 / 6) | 6 (R2) | 0 |
| | human d8 | 858 | 715 / 0 | 1 | 47 | 83 (71 / 3 / 9) | 11 | 1 |
| | human d32 | 1,060 | 916 / 3 | 0 | 39 | 97 (90 / 0 / 7) | 5 | 0 |
| | athal d1 | 384 | 220 / 1 | 0 | 90 | 72 (35 / 0 / 37) | 1 | 0 |
| | athal B d1 | 255 | 84 / 1 | 0 | 90 | 79 (37 / 0 / 42) | 1 | 0 |
| complete | human d1 | 3,692 | 1,753 / 435 | 1,138 | 6 | 127 (90 / 32 / 5) | 230 (R1, R1+R3, R3) | 3 |
| | human d8 | 1,221 | 633 / 100 | 297 | 1 | 66 (55 / 9 / 2) | 112 | 12 |
| | human d32 | 373 | 205 / 16 | 73 | 0 | 28 (20 / 7 / 1) | 47 | 4 |
| | athal d1 | 1,632 | 549 / 37 | 72 | 96 | 784 (507 / 38 / 239) | 88 | 6 |
| | athal B d1 | 1,333 | 215 / 59 | 75 | 95 | 789 (508 / 38 / 243) | 94 | 6 |

- **5' partials and internal CDSs: most remaining errors are selection errors, nearly all "no prediction on the transcript".** Human: 85-90% of 5' partials and 83-89% of internal CDSs. Start decisions cannot fix these.
  - The decision errors that remain are mostly "hit, no rule": the best hit aligns from deep inside the subject (sstart > 50 for most internal ones), so neither "N-terminus at the ATG" (R2) nor "alignment >= 5 aa upstream of the ATG" (R1) holds, and the ab initio rule trims wrongly.
  - Wrong R2 trims are rare: 4-12 per assembly.
  - Athal has more `no hit` (38-90) and Pfam-only cases, because Swiss-Prot minus Arabidopsis is sparse.
- **3' partials enter the start decision only through their 5' end.** By definition a 3' partial has no stop, so its ORF is 3'-open. It becomes a 5'-open candidate only when its 5' UTR has no in-frame stop: the longest ORF then runs open at both ends, and the right decision is to trim to the ATG.
  - This holds for 449 of 1,579 3' partials in human d1, 741 / 1,973 (d8), 657 / 1,531 (d32) and 188 / 593 (athal); the others are 5'-complete ORFs that the assist never sees.
  - Their decision errors are "kept open when the trim was right": R1 / R3 firing on a homolog whose alignment starts upstream of the annotated ATG (9-20), or the ab initio rule keeping it with a non-firing hit (9-33).
  - `outside` (5-61) are predicted ORFs whose start is not a keep/trim choice, or whose true start is a later ATG.
  - Most 3'-partial errors (73-81% in human) are again "no prediction".
- **Complete CDSs wrongly left open** (right stop, trim needed, kept open; `anatomy_complete_left_open.tsv`):

  | dataset | rule wrong | hit, no rule | no hit | conflict |
  |---|---|---|---|---|
  | human d1 | 230 | 127 | 6 | 3 |
  | human d8 | 112 | 66 | 1 | 12 |
  | human d32 | 47 | 28 | 0 | 4 |
  | athal d1 | 88 | 784 | 96 | 6 |
  | athal B d1 | 94 | 789 | 95 | 6 |

  - In human, R1 or R3 firing wrongly dominates, i.e. homologs with a longer N-terminus.
  - In athal, the dominant case is a hit that does not fire, with the ab initio rule keeping the end open: the alignment starts at subject position > 3 (507), or there is only a Pfam domain downstream of the ATG (239).
  - A deeper search cannot change these: the subject does not align to the region that the decision is about.
- Hit-no-rule reasons over all types:
  - human d1: sstart 4-10 / 11-50 / > 50 = 32 / 53 / 84; N-terminus mismatch 34; Pfam only 17;
  - athal d1: 126 / 235 / 275, 42, 339.

### Results: is the evidence there in deeper hits? (`evidence_summary.tsv`)
Labelled candidates, split by whether a5b24c3 decides them right or wrong. Counts are those with support (from any of the 25 hits, every HSP, e <= 10, more-sensitive, Pfam P <= 1e-4) for the right call only, the wrong call only, or both.

| dataset | a5b24c3 wrong: n / right only / both | a5b24c3 right: n / wrong only / both |
|---|---|---|
| human d1 | 533 / 37 / 44 | 9,177 / 69 / 1,174 |
| human d8 | 422 / 73 / 42 | 8,799 / 50 / 776 |
| human d32 | 292 / 51 / 28 | 8,077 / 16 / 460 |
| athal d1 | 1,369 / 178 / 36 | 6,163 / 373 / 250 |
| bench fragments (8) | 6,642 / 1,446 / 127 | 62,693 / 2,529 / 1,664 |
| bench complete (8) | 5,572 / 226 / 229 | 32,903 / 3,675 / 2,824 |

- Only 7-17% of the wrong decisions on the assemblies, 22% on the benchmark fragments and 4% on the complete inputs gain uncontested support for the right call from deeper hits.
- Deeper hits put contradicting evidence on many more of the right decisions: 6-14% on the assemblies, 7% on the fragments and 20% on the complete inputs.
- With the top hit at 1e-5 this is 0.1-0.8% (all of it "both").
- The `-k 1` hit and `-k 25` at 1e-5 give nearly the same picture as e <= 10 (the table file has every source). Looser e-values add few hits, because Swiss-Prot hits of these ORFs are either strong or absent.

### Results: variants (`family_by_dataset.tsv`, `selected_variants.tsv`)
Change in exact F1 (points) vs a5b24c3.
- Tuning = mean of human d1 and d8.
- Benchmark columns: median over 8 inputs, with (min, number down) for complete.

| variant | tuning | d32 | athal d1 | bench fragments | bench complete (min, n down) |
|---|---|---|---|---|---|
| **best on tuning:** msens, e<=10, all HSPs, top3 bitscore-weighted, Pfam dE1e-5 | +0.093 | +0.167 | +0.097 | +0.183 (2 down, min -0.35) | -0.027 (-1.17, 4/8 down) |
| best with Pfam GA: msens, e<=10, all HSPs, top3 bits, GA | +0.087 | +0.152 | +0.175 | +0.167 (1 down) | -0.019 (-0.24, 4/8 down) |
| best e-value x mode (top1, best HSP, GA): msens, e<=10 | +0.039 | +0.106 | +0.136 | +0.128 (0 down) | -0.015 (-0.20, 4/8 down) |
| msens, e<=1e-5 (top1) | +0.025 | +0.068 | +0.083 | +0.084 | +0.007 (-0.017, 1/8 down) |
| fast, e<=10 (top1) | -0.001 | +0.008 | +0.005 | +0.046 | 0.000 (-0.036, 3/8 down) |
| fast, e<=1e-10 (top1) | -0.014 | -0.023 | -0.019 | -0.032 | -0.004 (-0.012, 4/8 down) |
| best depth/scheme (fast, 1e-5, best HSP, GA): top5 anyup | +0.042 | 0.000 | +0.029 | +0.010 (3 down) | +0.021 (-0.066, 2/8 down) |
| top25 maj / unan (fast, 1e-5) | -0.097 / -0.133 | -0.015 / +0.008 | +0.044 / +0.039 | -0.008 / 0.000 | -0.065 / -0.050 |
| all HSPs instead of best (fast, 1e-5, top1) | +0.007 | +0.008 | 0.000 | 0.000 | 0.000 (-0.008, 3/8 down) |
| Pfam dE1e-5 instead of GA (k1) | -0.003 | +0.008 | -0.087 | +0.011 | -0.027 (-0.95, 7/8 down) |
| Pfam dE1 instead of GA (k1) | -0.216 | -0.038 | -1.257 | -0.309 | -1.613 (-8.86, 8/8 down) |
| no R3 (k1) | -0.135 | -0.296 | -0.534 | -1.009 | +0.190 (0/8 down) |
| oracle on all candidates (ceiling) | +1.60 / +2.16 | +2.22 | +6.64 | +4.79 | +3.29 |

- **Top hit vs top n.**
  - Taking the top 3-5 hits gains at most +0.04 (tuning) on the a5b24c3 search. Every scheme loses at n = 10-25: majority and unanimity are diluted by weaker hits that start elsewhere in the subject.
  - The best consistency scheme is "R2 from any hit trims, otherwise any R1 keeps" (`anyup`). It is still neutral on the held-out sets: +0.00 on d32, +0.03 athal, +0.01 bench fragments, +0.02 complete.
- **E-value.**
  - With the default (fast) DIAMOND mode, loosening from 1e-5 to 10 changes nothing measurable (|change| <= 0.05 points), because the top hit is the same.
  - 1e-10 is slightly worse everywhere (-0.01 to -0.03).
  - More-sensitive mode adds hits for athal and the benchmark fragments: +0.08 to +0.14 points at the top hit. It costs on complete inputs at e <= 10 (tthermophila -0.20, creinhardtii -0.08), less at e <= 1 (min -0.06) and little at 1e-5 (min -0.017).
- **All HSPs** instead of the best HSP of the top hit: no effect (<= 0.01 points difference).
- **Pfam.**
  - Any E-value threshold looser than about GA hurts, sharply on complete inputs: dE1 -1.6 median, -8.9 tthermophila. R3 then keeps UTR extensions open on weak domain hits.
  - dE1e-5 (P <= 1e-9, about as strict as GA) is about neutral on the assemblies, but 7 of 8 complete inputs go down (to -0.95).
  - GA is the best Pfam setting. Dropping R3 costs 0.3-1.0 points on the assemblies and fragments, while gaining +0.19 on complete inputs (the step-3 trade-off).
- **The best variant chosen on the tuning sets** (more-sensitive, e <= 10, all HSPs, top 3 bitscore-weighted) is consistently positive on the held-out assemblies and fragments: +0.10 to +0.18. Per type it mostly helps 5' partials and internal CDSs, e.g. athal internal +1.4 to +1.6, 5' partial +0.9 to +1.4.
  - It lowers 4 of the 8 complete inputs: tthermophila -0.24 (GA) or -1.17 (dE1e-5), creinhardtii -0.18, dmel -0.06, athal -0.06.
  - 3'-partial F1 on the fragments falls in 5 of 8 inputs (to -0.63, tthermophila; GA version).
  - It needs a more-sensitive search with 25 targets and all HSPs: 2-5x the DIAMOND time of default mode at the same `-k 25` (113 s vs 58 s on human d1, 82 s vs 17 s on tthermophila fragments, 8 threads), and ~17x the output lines of `-k 1`.

### Conclusion (step 4)
- **Most of what the assist still gets wrong is not a start decision.**
  - On the human assemblies, 83-90% of the non-exact 5'-partial and internal CDSs, and 76-82% of 3' partials, have no prediction with the true stop (nearly always no prediction on the transcript at all).
  - Of the decision errors that remain, the main class is "hit, no rule": the best homolog aligns from inside the subject (sstart > 3), so it says nothing about the 5' end.
  - Complete CDSs wrongly left open come mainly from R1 / R3 in human (homologs with longer N-termini) and from non-firing hits in athal.
- **3' partials** are 5'-open candidates only when the 5' UTR has no in-frame stop (28-43% of them). There the right call is always "trim", and errors come from keep calls.
- **The top hit at 1e-5 is not too stringent.**
  - Deeper or looser hits rarely supply new uncontested evidence for the wrong decisions (4-22%).
  - They add contradicting evidence to 6-20% of the right ones.
  - Varying the e-value (1e-10 to 10), the number of hits (1-25, five consistency schemes), the HSPs used, and Pfam GA vs E-value thresholds changes exact F1 by at most +0.18 points on the held-out assemblies and fragments.
  - Every variant that gains there lowers some complete benchmark inputs, by up to 0.24 points (Pfam GA) or more (looser Pfam).
- **Recommendation: keep the a5b24c3 settings** (DIAMOND `-k 1 --evalue 1e-5`, first HSP; hmmsearch `--cut_ga`). No variant beats them on the held-out sets without making complete inputs worse.
  - The available gain is about +0.1 to +0.2 exact-F1 points: 20-36 CDSs per held-out assembly, out of an oracle ceiling of +1.6 to +6.6. It needs a 2-5x slower search.
  - The remaining ceiling sits in hits that do not reach the subject's N-terminus and in selection, not in hit depth.

### Caveats (step 4)
- Pfam "domE" thresholds are defined at a reference domZ of 10,000. A user's `hmmsearch --domE x` on all LongOrfs has a domZ set by the number of significant sequences, so their `--domE 1` is not exactly `dE1` here. The result (looser than GA hurts; dE1 is far looser than GA) does not depend on this.
- Variant selection used human d1 + d8 only. The held-out human d32 shares genes with them.
- The benchmark scoring substitutes decisions in the `pytd_v7` ab initio predictions, as step 3. Selection with homology hits is not re-run, so these are start-decision effects only.
- Swiss-Prot without the organism's own entries, as before.
- Rule geometry (R2 window sstart <= 3, R1 offset 5) was not varied here. Step 3's grid covers it.

## Step 5: UniRef50 vs Swiss-Prot (2026-10-29)

### Question
Brian asked whether UniRef50 (and possibly UniRef90) works better than Swiss-Prot as the DIAMOND database of
homology mode, in particular for the homology start assist. UniRef50 has ~67x more sequences and more diversity, so
more 5'-open ORFs may get evidence; but its unreviewed and fragmentary entries may bring false positives and
unreliable N-termini.

### Database (UniRef50 release 2026_03, 38,840,027 clusters)
- **Build** (`pipeline/start_classifier/step5_01_uniref50_build_v3.sh`, log `data/start_classifier/step5/build_v3.log`):
  streamed from the UniProt FTP, no FASTA or XML on disk.
  - Pass 1, `uniref50.xml.gz`: per cluster the representative's NCBI taxid and the benchmark organisms among
    all members (`step5_01_uniref50_xml_v2.awk`). 38,840,027 rows, as the release note.
  - Pass 2, `uniref50.fasta.gz`: headers cut to the cluster id, piped into `diamond makedb` with
    `--taxonmap` (representative taxid) and NCBI taxdump `nodes.dmp` / `names.dmp` (downloaded 2026-10-29);
    each sequence's length, first residue, `n=` and `TaxID=` recorded on the way.
  - DB: 12,182,780,325 letters, all 38,840,027 sequences mapped to a taxid, 13.9 GB.
- **Why pass 1:** the FASTA header's `TaxID=` is the cluster's lowest common taxon, not the representative's
  (e.g. `UniRef50_Q8WZ42`, titin, `n=1909 TaxID=7742`, representative human). Excluding by it would keep
  clusters whose representative is the query organism's own protein.
- **Versions:**
  - v1 (regex awk) was stopped for speed; same output.
  - v2 left 55,666 representative taxids that NCBI has merged unknown to `nodes.dmp`. DIAMOND
    `--taxon-exclude` silently drops sequences with unknown taxids: 55,666 extra sequences (0.14%) in every
    search. Found in a pilot ("Filtered database contains 38,732,414 sequences" for human).
  - v3 resolves them through `merged.dmp` (55,657; 4 more get the common taxid, 5 the root) and overwrote the
    v2 db. With v3 the human filter removes exactly the 51,947 human-represented clusters.
- **Exclusion, mirroring `sprot_minus_<taxid>`:** those databases (no builder saved, E07 Gaps) remove exactly
  the entries with the organism's own `OX=` taxid. The Swiss-Prot counts match: 575,748 entries; minus 9606
  leaves 555,317 (20,431 `OX=9606`); 3702 and 7227 likewise exact. Descendant taxa were kept, e.g. one
  Neanderthal entry.
  - `diamond --taxon-exclude <taxid>` removes the taxon and its NCBI descendants (`step5_00_own_taxsets_v1.py`).
    Human has 2 descendant taxa, mouse 17 subspecies, C. reinhardtii 1, the others none.
  - It works on the representative's taxid. `hits/uniref50/subject_check.tsv`: 0 hits with an excluded
    representative taxid in all 16 searches.
- **Clusters removed** (`own_cluster_stats.tsv`):

  | organism | clusters removed | their members | clusters with a member of the organism, kept (other representative) |
  |---|---|---|---|
  | human | 51,947 | 7,236,815 | 44,111 |
  | mouse | 25,138 | 2,694,683 | 26,379 |
  | athal | 19,677 | 2,106,696 | 9,898 |
  | dmel | 12,139 | 845,641 | 6,602 |
  | tthermophila | 18,790 | 48,632 | 5,044 |
  | creinhardtii | 9,547 | 73,569 | 7,093 |
  | spombe | 4,796 | 178,025 | 205 |
  | pfalciparum | 2,145 | 51,318 | 2,997 |

  - **This is not the same experiment as Swiss-Prot minus taxon.** Removing a human-represented cluster
    removes every member within 50% identity of it, i.e. the human protein's orthologs in other mammals.
  - What remains near a human query is:
    - clusters with a non-human representative (some contain human members);
    - and, because UniRef requires 80% overlap with the cluster seed, separate clusters of length variants
      of the same protein, e.g. N-terminally extended gene models.
  - Swiss-Prot minus human, in contrast, keeps the reviewed chimp or mouse ortholog.
  - **Near-identical representatives:** as with sprot_minus (chimp kept for human), representatives from
    other taxa can be near-identical to the organism's own protein. That is realistic for a user.
  - **How often the top -k 1 hit is >= 95% identical** (5'-open candidates, `hit_stats.tsv`):

    | input | Swiss-Prot | UniRef50 |
    |---|---|---|
    | human simulated assemblies | 0.33-0.46 | 0.30-0.49 |
    | human benchmark | 0.33-0.37 | 0.25-0.38 |
    | mouse benchmark | 0.32-0.35 | 0.22-0.28 |
    | athal | 0.011-0.014 | 0.08-0.16 |
    | dmel | 0.02-0.04 | 0.08-0.16 |
    | other organisms | 0.005-0.07 | 0.007-0.17 |

  - The UniRef50 top hit's cluster contains a member of the query organism for:
    - human 21-27%, mouse 24-34%;
    - athal 9-15%, dmel 27-34%, tthermophila 31-33%, pfalciparum 42-46%, creinhardtii 29-30%, spombe ~1%.
  - Top hit >= 95% identical and own member in its cluster: human 10-15%, athal 1-2%.

### Searches
- **Queries** (`step5_02_queries_v2.py`): identical peptides were searched once per organism; ORF ids were
  mapped back afterwards (`step5_05_expand_hits_v1.py`).
  - Candidates: the 5'-open candidate peptides of 21 inputs: simulated human d1 / d8 / d32, athal d1 / d16,
    and the benchmark complete and fragments inputs of the 8 organisms. 167,190 peptides, 162,709 unique.
    - These are the step-4 `cand.pep` files.
    - Athal d16 is new: candidates from `01_features.py` on E32's d16 `pytd_v7` run
      (`$P/data/start_classifier/step5/features/`, 14,242 5'-open, rule re-derivation 0 mismatches).
  - Full mode: all LongOrfs peptides of athal d1 / d16 and human d8 (the a5b24c3 LongOrfs; checked `cmp`
    identical to the release runs'), and of bench human / athal / dmel / tthermophila complete, fragments,
    `.wRand`, plus human `.wLnc`. The `pytd_v7_homology` LongOrfs, which are identical to `pytd_v7h_homology`.
  - 2,142,066 peptides, 1,167,893 unique per organism.
- **UniRef50 searches** (`step5_04_uniref50_searches_v2.sh`; log `uniref50_searches_v2.log`; 12:16-14:22):
  - DIAMOND 2.2.8, default (fast) mode, `--taxon-exclude <taxid>`, 8 threads, nice 19, one search at a time;
  - `-k 25 --max-hsps 0 --evalue 1e-5` on the candidates, and `-k 1 --evalue 1e-5` (as `run_tool.sh`) on
    full-mode LongOrfs plus all candidates;
  - v1 of the driver looked for the wrong query file name and was stopped at once.
- **Swiss-Prot additions** (`step5_03_sprot_searches_v2.sh`):
  - athal d16 candidates: fast `-k 25 --max-hsps 0 --evalue 10`, as step 4;
  - human d8 full mode: `-k 1 1e-5` on all 100,141 LongOrfs peptides, and hmmsearch `--cut_ga` (4 chunks);
  - v1 was stopped before its first search ended, to keep within 8 threads while the db was being built.
- **Check, query-set independence:** the pooled human search ran DIAMOND's double-indexed mode; a single-input
  search of bench.human.complete (query-indexed) gave 56 of ~25,850 lines different, nearly all weak hits.
  28 extra queries had a hit, with e-values 1e-6 to 1e-15.

### Method
- **Start-decision isolate** (`step5_06_rule_variants_v1.py`; `step5_07_summarize_v1.py`): as step 4.
  - The release run's predictions, with each candidate's keep/trim decision substituted.
  - Exact F1 overall and per type with negatives.
  - Same Pfam GA hits for both databases.
  - Baseline: a5b24c3 rules on the Swiss-Prot `-k 1` hits. It reproduces the adoption values, e.g. human d1
    0.8641, R1 3,199 / R2 6,793.
  - Variants per database (Swiss-Prot `-k 25` = the step-4 fast search filtered to 1e-5):
    - a5b24c3 rules on the top hit (`k1`, and `k25` top1);
    - top-n R1/R2 schemes (n = 3, 5, 10, 25; maj / unan / bits / anyup / first, as step 4);
    - **consensus N-terminus:**
      - each hit projects subject residue 1 to p = qstart - (sstart - 1);
      - class trim if |p - a| <= 5, keep if p < a - 5 (p = a - 5 counts as trim);
      - bitscore-weighted majority of the top n (1-25), denominator all hits or firing hits;
      - optionally only hits with sstart <= 30;
    - **M filter:**
      - `r2`: a trim class needs a subject beginning with M;
      - `all`: non-M subjects dropped.
      - The first residue was recorded while streaming the FASTA (Swiss-Prot from `uniprot_sprot.fasta.gz`).
  - Selection on simulated human d1 + d8 only. Held out: d32, athal d1 / d16, bench fragments and complete.
- **Full homology mode** (`step5_08_full_predict_v1.sh`; driver `step5_08_full_predict_all_v1.sh`, 14:23-14:57, 21/21 rc 0):
  - fresh PyTD a5b24c3 Predict with the UniRef50 `-k 1` hits on all LongOrfs plus the same Pfam hits;
  - the LongOrfs files are linked, not copied (Predict only reads them);
  - Swiss-Prot arm: the existing `pytd_v7h_homology` (bench) and athal depth-series B runs; human d8 was run
    fresh with the new Swiss-Prot full searches;
  - scored with `step5_09_full_scores_v2.py`:
    - as `16_score_adoption.py`, excluded contigs removed on the assemblies;
    - stop sensitivity = tp_stop / n_truth, stop precision = tp_stop / n_pred;
    - negatives: `rand_*` / `lnc_*` transcripts, or the assemblies' lncRNA-only contigs.
  - The athal Swiss-Prot runs reproduce the depth-series exact F1 (d1 0.91383, d16 0.89887).
- **Post hoc (after seeing the results, not tuned):** the "hybrid" input, Swiss-Prot hits plus UniRef50 hits for
  ORFs with no Swiss-Prot hit. A user gets it by concatenating the two outfmt6 files, Swiss-Prot first, because
  PyTD reads the first line per query.
  - Isolate: `step5_13_hybrid_v1.py`.
  - Full mode: `step5_08_full_predict_v2.sh` / `_all_v2.sh`, 15:18-15:48, 21/21 rc 0.
- **Diagnosis:** `step5_10_diagnose_v1.py`, the candidates whose a5b24c3 call differs between databases.
- **Tables:** `step5_12_tables_v2.py` writes `isolate_compact.tsv` and `full_mode_compact.tsv`.

### Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; B=$P/data/start_classifier/step5; S=$T/data/start_classifier/step5; D=$T/pipeline/start_classifier; PY=$P/env/startclf/bin/python
bash $D/step5_01_uniref50_build_v3.sh > $B/build_v3.log 2>&1      # pass 1 reused from v2 (same output)
PYTHONPATH=$P/work/trees/PY_v7 $PY $D/01_features.py $P/work/simtrinity_athal/depth/d16/runs/pytd_v7 Universal \
  $P/data/simtrinity_athal/depth/d16/inputs/fragments_simtrinity.truth.tsv sim.athal.d16 $B/features/sim.athal.d16.tsv.gz
python3 -I $D/step5_02_queries_v2.py $B/queries
bash $D/step5_03_sprot_searches_v2.sh 4 > $B/sprot_searches_v2.log 2>&1
bash $D/step5_04_uniref50_searches_v2.sh human athal dmel mouse spombe pfalciparum creinhardtii tthermophila > $B/uniref50_searches_v2.log 2>&1
$PY -I $D/step5_05_expand_hits_v1.py $B
$PY $D/step5_06_rule_variants_v1.py $S/step5_manifest.tsv $B/variants
$PY $D/step5_07_summarize_v1.py $B/variants/variants_f1.tsv.gz $B/variants/hit_stats.tsv $S > $S/summary.stdout.txt
bash $D/step5_08_full_predict_v1.sh sim.human.d8 sprot; bash $D/step5_08_full_predict_all_v1.sh; bash $D/step5_08_full_predict_all_v2.sh
$PY $D/step5_10_diagnose_v1.py $S/step5_manifest.tsv $B/taxonomy/names.dmp $S   # names.dmp now gzipped
$PY $D/step5_13_hybrid_v1.py $S/step5_manifest.tsv $S
$PY $D/step5_09_full_scores_v2.py $S; $PY $D/step5_12_tables_v2.py $S > $S/tables.stdout.txt
bash $D/step5_11_timing_v1.sh > $B/timing_v1.log 2>&1
```
- **Outputs in the repo** (`data/start_classifier/step5/`):
  - isolate: `main_by_dataset.tsv`, `tuning_ranking.tsv`, `selected_variants.tsv`, `selected_by_dataset.tsv`,
    `hit_stats.tsv`, `hybrid_f1.tsv`, `isolate_compact.tsv`;
  - full mode: `full_scores.tsv`, `full_mode_compact.tsv`;
  - diagnosis: `diag_*.tsv`;
  - logs.
- **Large files** in `$B`:
  - `db/`: the dmnd, taxonmap, seqinfo, XML table;
  - `queries/`, `searches/`, `hits/` (per-dataset hit files, gzipped after use; the `hits_used.tsv` of the
    runs names the plain paths), `variants/variants_f1.tsv.gz`.
- **Runs:** `$P/work/start_classifier/step5_full/{uniref50,hybrid,sprot}/<dataset>/`.

### Results: coverage and rule firings (start-decision isolate, `hit_stats.tsv`, `main_by_dataset.tsv`)
**5'-open candidates with any hit at 1e-5, Swiss-Prot -> UniRef50:**
- human 0.93-0.97 -> 0.98-1.00;
- mouse 0.95 -> 0.99;
- athal 0.49-0.59 -> 0.98-0.99;
- dmel 0.51-0.54 -> 0.98;
- spombe 0.61-0.77 -> 0.88-0.96;
- tthermophila 0.18-0.27 -> 0.72-0.81;
- pfalciparum 0.29-0.37 -> 0.94-0.99;
- creinhardtii 0.29-0.37 -> 0.87-0.91.

**What the top hit looks like:**
- The top hit's subject begins at residue <= 3 less often with UniRef50 in human (d1 0.61 -> 0.44), and more often
  where Swiss-Prot was sparse (athal d1 0.19 -> 0.37).
- The top hit's subject begins with M in 0.86-0.98 of UniRef50 hits, against 0.97-1.00 for Swiss-Prot.

**Rule firings, a5b24c3 rules on the `-k 1` hit, Swiss-Prot -> UniRef50:**

| dataset | R1 keep | R2 trim | wrong keeps | wrong trims |
|---|---|---|---|---|
| human d1 | 3,199 -> 5,080 | 6,793 -> 4,960 | 394 -> 1,659 | 139 -> 107 |
| human d8 | 4,829 -> 6,087 | 4,255 -> 3,089 | 229 -> 1,068 | 193 -> 155 |
| human d32 | 6,077 -> 6,844 | 2,200 -> 1,556 | 108 -> 587 | 184 -> 161 |
| athal d1 | 928 -> 2,664 | 686 -> 2,531 | 1,035 -> 1,034 | 334 -> 184 |
| athal d16 | 3,708 -> 9,111 | 533 -> 1,491 | 434 -> 494 | 1,148 -> 378 |
| bench human complete | 89 -> 1,408 | 6,263 -> 4,525 | 316 -> 1,622 | 0 -> 0 |
| bench athal fragments | 4,040 -> 8,761 | 330 -> 1,802 | 729 -> 643 | 891 -> 344 |
| bench tthermophila fragments | 1,228 -> 5,838 | 104 -> 1,741 | 839 -> 734 | 1,480 -> 813 |

- Wrong keep = kept open where the truth is a UTR extension with the true start at the first ATG.
- Wrong trim = trimmed where the truth is truncated.
- R3 (Pfam) is identical by design.

### Results: start-decision isolate, exact F1 (`isolate_compact.tsv`)
Change in points vs the a5b24c3 baseline (Swiss-Prot `-k 1`).
- **UniRef50 a5b24c3**: the a5b24c3 rules on the UniRef50 `-k 1` hit.
- **best UniRef50 variant**: `uniref50|k25|top25_maj|all`, the top 25 hits with non-M subjects dropped, majority
  of R1/R2 calls; chosen on human d1 + d8.
- **hybrid** (post hoc): Swiss-Prot hit if any, else the UniRef50 hit.

| dataset | UniRef50 a5b24c3 | best UniRef50 variant | hybrid | hybrid, fallback top25 maj | perfect decisions |
|---|---|---|---|---|---|
| sim human d1 (tune) | -3.70 | -0.31 | +0.08 | +0.09 | +1.60 |
| sim human d8 (tune) | -4.11 | -0.32 | +0.22 | +0.23 | +2.16 |
| sim human d32 | -3.46 | -0.30 | +0.27 | +0.30 | +2.22 |
| sim athal d1 | +0.73 | +2.90 | +0.87 | +1.64 | +6.64 |
| sim athal d16 | +3.81 | +4.83 | +3.68 | +3.85 | +8.48 |
| bench human complete / fragments | -6.90 / -3.52 | -0.67 / -0.26 | -0.11 / +0.18 | +0.03 / +0.25 | +1.67 / +1.88 |
| bench mouse | -2.26 / -1.10 | -0.32 / -0.13 | +0.08 / +0.32 | +0.16 / +0.36 | +1.51 / +1.62 |
| bench dmel | -0.01 / +1.69 | +0.63 / +2.25 | +0.04 / +1.50 | +0.46 / +1.96 | +2.32 / +4.48 |
| bench athal | +0.16 / +2.61 | +2.37 / +3.39 | +0.34 / +2.36 | +1.18 / +2.84 | +6.75 / +6.68 |
| bench spombe | +0.85 / +1.10 | +0.91 / +1.27 | +0.48 / +1.10 | +0.54 / +1.16 | +5.34 / +5.02 |
| bench tthermophila | +0.59 / +3.22 | +0.26 / +3.34 | +0.38 / +2.96 | +0.31 / +3.18 | +7.15 / +9.67 |
| bench pfalciparum | +0.25 / +2.54 | +0.40 / +2.58 | +0.19 / +2.48 | +0.32 / +2.52 | +0.72 / +4.55 |
| bench creinhardtii | -0.66 / +0.83 | -1.33 / +0.71 | -0.57 / +0.74 | -0.91 / +0.69 | +4.26 / +5.91 |

Benchmark rows give complete / fragments.

**Per type, UniRef50 a5b24c3 vs Swiss-Prot** (`selected_variants.tsv`):
- human: complete -4.3 / -6.4 / -9.3, 3' partial -7.0 / -6.1 / -8.9, 5' partial +0.1..+0.6 and internal
  +0.4..+1.0 (d1 / d8 / d32). The losses are trims that became keeps.
- athal: 5' partial +6.1 / +6.3, internal +4.5 / +6.8, 3' partial +2.9 / +1.7, complete -0.1 / -1.5 (d1 / d16).

**Top-n and the other variants on UniRef50** (`tuning_ranking.tsv`; tuning = mean of human d1 + d8):
- No UniRef50 variant reaches the Swiss-Prot baseline on the human tuning sets; the best is -0.31.
- Majority over the top 25 hits dilutes the extended-N-terminus hits. The best variant, chosen on human, is
  then the best UniRef50 variant elsewhere too:
  - held out: d32 -0.30, athal +2.90 / +4.83;
  - bench fragments median +1.76 (2 of 8 down, min -0.26);
  - bench complete median +0.33 (3 of 8 down: creinhardtii -1.33, human -0.67, mouse -0.32).
- Other variants:
  - top-5 anyup -0.47;
  - top-25 unanimous -1.34;
  - top-25 majority without the M filter -0.38; with `r2` -0.44.
- **Consensus N-terminus rule: worse.**
  - Top 1: -6.2 on tuning. Projecting subject residue 1 through an unaligned subject N-terminus is unreliable.
  - The best form (top 25, firing hits, sstart <= 30): tuning -1.37, human d32 -4.26, athal +2.57 / +0.88,
    bench fragments median -0.95 (6 of 8 down), complete median +0.37.
  - It trims truncated ORFs (human d32: 743 wrong trims vs 184).
  - With Swiss-Prot hits it also loses: best `sprot|k25|cons1_all` tuning -0.17, bench complete median -1.76.
- **M filter: almost no effect.**
  - On the top hit: human d1 -3.70 -> -3.74 (`r2`) / -3.52 (`all`); bench complete median +0.08 -> +0.07 / +0.15.
  - The misleading UniRef50 subjects mostly do begin with M (diagnosis below). They are full-length entries
    with longer N-termini, not fragments.
- **Swiss-Prot top-n**, for reference: the best Swiss-Prot variant is top-5 anyup, +0.04 on tuning and within
  +-0.2 elsewhere, as step 4.

### Diagnosis: why UniRef50 loses on human and mouse (`diag_r1_wrong.tsv`, `diag_r1_species.tsv`, `diag_flips.tsv`)
**UTR extensions that UniRef50 keeps open by R1** (true start at the first ATG):

| dataset | n | Swiss-Prot called R2 (trim, right) on the same candidate | median identity | subject begins with M | median cluster size | subject's residue 1 maps upstream of the ATG (median) | own member in cluster |
|---|---|---|---|---|---|---|---|
| human d1 | 1,456 | 78% | 90% | 91% | 11 | 101 aa | 19% |
| human d8 | 915 | 80% | 91% | 90% | 11 | 101 aa | 18% |
| human d32 | 518 | 86% | 93% | 89% | 9 | 93 aa | 18% |
| bench human complete | 1,382 | 81% | 90% | 92% | 9 | 115 aa | 14% |

- Only 10-14% of these clusters are singletons.
- Their representatives are unreviewed mammalian entries: macaque, horse, marmoset, dog, crab-eating macaque,
  rabbit, orangutan and others, each 2-5%.
- These are orthologs with N-terminal extensions of about 100 aa: alternative upstream starts in other species'
  gene models, or isoforms. UniRef keeps them out of the human ortholog's cluster because of the 80%-overlap
  rule, so excluding human-represented clusters leaves them as the best hit.
- **Candidates flipped from right (Swiss-Prot) to the opposite call (UniRef50):** 1,141 in human d1, 1,122 in
  bench human complete, against 31 the other way.
- **Athal and the sparse organisms:**
  - Most changes are candidates that had no Swiss-Prot hit and now get a call: 2,860-5,980 per input
    (athal, dmel, tthermophila, pfalciparum, creinhardtii).
  - Right-to-wrong flips are 0-104 per input.
- **The cost on complete inputs where it remains** (creinhardtii -0.66, dmel -0.01): more R1 keeps.

### Results: full homology mode (`full_scores.tsv`, `full_mode_compact.tsv`)
Fresh a5b24c3 Predict, all LongOrfs searched, same Pfam hits. Differences UniRef50 minus Swiss-Prot.

| input | exact F1 Swiss-Prot -> UniRef50 (hybrid) | stop sens. | stop prec. | predictions on negatives (ab initio / Swiss-Prot / UniRef50 / hybrid) |
|---|---|---|---|---|
| bench human complete | 0.9591 -> 0.8883 (0.9580) | 0.9947 -> 0.9973 | 0.9986 -> 0.9984 | - |
| bench human fragments | 0.9428 -> 0.9013 (0.9424) | 0.9474 -> 0.9540 | 0.9977 -> 0.9920 | - |
| bench human complete.wRand | 0.9563 -> 0.8858 (0.9552) | same as complete | 0.9928 -> 0.9926 | 113 / 113 / 113 / 113 of 19,167 |
| bench human complete.wLnc | 0.8832 -> **0.7310** (0.7884) | 0.9943 -> 0.9972 | 0.8527 -> **0.6985** | 1,639 / 3,268 / **8,220** / 8,231 of 33,668 |
| bench human fragments.wLnc | 0.8612 -> **0.7366** (0.7699) | 0.9460 -> 0.9534 | 0.8360 -> 0.6814 | 1,928 / 3,522 / **8,389** / 8,400 |
| sim human d8 | 0.9209 -> 0.8614 (0.9055) | 0.9808 -> 0.9884 | 0.9397 -> 0.9006 | 704 / 1,129 / 2,026 / 2,028 of 5,885 lncRNA contigs |
| bench athal complete / fragments | 0.9066 -> 0.9135 (0.9156) / 0.8475 -> 0.8963 (0.8945) | 0.974 -> 0.987 / 0.853 -> 0.902 | -0.1 / -0.2 pt | wRand: 39 / 58 unchanged |
| bench dmel complete / fragments | 0.9310 -> 0.9385 (0.9394) / 0.8794 -> 0.9056 (0.9037) | 0.969 -> 0.986 / 0.888 -> 0.909 | -0.01 / -0.1 pt | wRand: 44 / 53 unchanged |
| bench tthermophila complete / fragments | 0.8290 -> 0.8359 (0.8338) / 0.7785 -> 0.8233 (0.8208) | 0.894 -> 0.900 / 0.820 -> 0.848 | -0.5 / -0.4 pt | wRand: 61 / 126 unchanged |
| sim athal d1 | 0.9138 -> 0.9249 (0.9267) | 0.9764 -> 0.9935 | 0.9918 -> 0.9836 | 83 / 93 / 265 / 265 of 1,372 |
| sim athal d16 | 0.8989 -> 0.9496 (0.9483) | 0.9702 -> 0.9966 | 0.9992 -> 0.9983 | 4 / 4 / 18 / 18 of 154 |

**Per type, exact F1 with negatives, UniRef50 - Swiss-Prot:**
- human d8: complete -6.8, 5' partial +0.5, 3' partial -6.4, internal +1.3;
- athal d1 / d16: 5' partial +8.3 / +8.9, 3' partial +7.0 / +5.2, internal +6.3 / +7.9, complete +0.3 / -1.0.

**Wall time for 21 Predict runs:** about the same with either database (114-1,005 s each).

**Retention:**
- UniRef50 hits retain more ORFs: stop sensitivity rises everywhere, +0.3 to +4.8 points.
- But lncRNA ORFs also gain hits, from unreviewed predicted proteins of other species:
  - lncRNA calls rise 2.5x on the human benchmark (3,268 -> 8,220);
  - 1.8x on the human d8 assembly's lncRNA contigs (1,129 -> 2,026);
  - 2.8x on athal d1 (93 -> 265; 1,372 lncRNA-only contigs).
- Random GC-matched negatives (`.wRand`) gain no calls at all.
- The hybrid keeps the start decisions of Swiss-Prot where Swiss-Prot has a hit, but it adds the same UniRef50
  hits for retention. So it removes the human start-decision loss (bench complete -0.11 instead of -7.08) but
  not the lncRNA false positives (8,231 lncRNA calls; wLnc exact F1 -9.5 / -9.1; human d8 -1.5).

### DIAMOND run time and disk
**Same input, as a user would run it:** bench.human.complete, 63,551 LongOrfs peptides, `-k 1 --evalue 1e-5`,
8 threads, back to back on the shared machine (load ~30-50 on 28 cores) (`step5_11_timing_v1.sh`, `timing_v1.log`):

| database | wall | CPU | peak RAM |
|---|---|---|---|
| Swiss-Prot minus human | 37.5 s | 125 s | 0.85 GB |
| UniRef50 minus human | 13 min 22 s (21x) | 2,757 s (22x) | 6.9 GB |

- About 3 min of each UniRef50 search is fixed cost: opening the 13.9 GB db and building the reference seed
  index for its 7 blocks. A 2,000-peptide pilot took 3 min.
- Pooled per-organism searches: 3-20 min each; 2 h 06 min for all 16 (167k + 1.17M unique peptides).

**Disk:**
- UniRef50 db 13.9 GB, plus 1.1 GB of tables (taxon map, per-sequence first residue and length, XML-derived
  representative taxids).
- Swiss-Prot db: 0.28 GB.

### Conclusion (step 5)
**UniRef50 minus taxon is not uniformly better than Swiss-Prot minus taxon. It depends on how well Swiss-Prot
covers the organism's clade.**
- **Sparse in Swiss-Prot** (Arabidopsis, Drosophila, Tetrahymena, Plasmodium, S. pombe fragments):
  - UniRef50 lifts the candidates with a hit from 18-77% to 72-99%;
  - start-decision exact F1 rises by +0.7 to +3.8 points on the assemblies and fragments with the a5b24c3 rules;
  - full homology mode gains +0.7 to +5.1 points, every type except complete CDSs improving.
- **Human and mouse:**
  - UniRef50 is clearly worse: start-decision exact F1 -1.1 to -6.9; full mode (human only) -4 to -7 points;
  - on lncRNA-containing inputs exact F1 falls 12-15 points, because the lncRNA calls rise 2.5x.
  - Two mechanisms, both from unreviewed entries:
    - **N-terminally extended homologs:** mammalian orthologs with ~100 aa longer N-termini become the top
      hit and fire R1 (keep open) where the true start is the first ATG;
    - **lncRNA ORFs** gain hits to predicted proteins and are retained.
- **Fragments are not the problem**, so requiring an M at the subject's start does nothing.
- **Top-n consensus** (majority of 25) reduces the start-decision loss on human to -0.3 and raises the gain
  elsewhere, but needs `-k 25` hits and still leaves human below Swiss-Prot.
- **The consensus N-terminus projection is worse** than the a5b24c3 rules with either database.
- **Caveat on the exclusion:** removing human-represented clusters also removes their mammalian orthologs.
  This makes the human UniRef50 test harsher than Swiss-Prot minus human, which keeps reviewed orthologs.
  - The lncRNA effect and the extended-N-terminus effect are, however, properties of UniRef50's unreviewed
    content.
  - A human user with the full UniRef50 would mostly hit the human-represented cluster.
- **Post hoc hybrid** (Swiss-Prot first, UniRef50 fallback), isolate:
  - no loss on human (+0.1 to +0.3; bench complete -0.11);
  - most of the UniRef50 gain elsewhere: athal +0.9 / +3.7, bench fragments +0.2 to +3.0.
  - Worst complete input creinhardtii -0.57.
  - In full mode its UniRef50 hits still retain lncRNA ORFs.

### Recommendation
- **Do not replace Swiss-Prot with UniRef50 as the default or recommended database.**
- **Docs:** UniRef50 may be mentioned as an option for organisms far from well-curated model organisms
  (plants, protists, many invertebrates), with two warnings:
  - it raises lncRNA / non-coding false positives about 2-3x;
  - for mammals it worsens start decisions.
  - The practical variant is to supply Swiss-Prot hits first and UniRef50 hits only for ORFs without a
    Swiss-Prot hit (concatenated outfmt6). That keeps Swiss-Prot's start decisions but not its retention
    specificity.
- **A clean option needs a code change, not tested here:** separate hit files for retention (Swiss-Prot) and
  for the start assist (Swiss-Prot + UniRef50 fallback).
- **UniRef90 is not worth testing now.**
  - Size: 121,494,192 clusters (3.1x); FASTA 32.1 GB gz (3.66x UniRef50's 8.78 GB); XML 46.6 GB gz.
  - Expected db about 45 G letters and **about 51 GB dmnd**, plus about 3.5 GB of tables. The build streams
    ~79 GB from the FTP (about 2-2.5 h).
  - Expected search time about 3.7x UniRef50 in default mode: about **50 min per benchmark-size input** (vs
    13 min, and 0.6 min for Swiss-Prot), and about 8 h for this step's searches at 8 threads.
  - It would not raise coverage much: UniRef50 already gives hits to 98-99% of the athal / dmel candidates.
  - It would add more of the unreviewed redundancy that causes both failure modes. UniRef90 clusters are
    tighter, so even closer length-variant entries would become top hits.
  - The quota does not allow it at present.

### Caveats (step 5)
- Exclusion by the representative's taxid (and its descendants) is the closest analogue of `sprot_minus` but
  not equivalent (see "Clusters removed").
- Cluster-member flags for organisms count NCBI-merged member taxids only under their old ids. Four merged
  ids map to the 8 organisms (C. reinhardtii 2, mouse 1, dmel 1), so the own-member fractions may be slightly
  low.
- The UniRef50 searches pooled the queries per organism (double-indexed mode); single-input searches differ in
  about 0.2% of lines, all weak hits.
- The benchmark Swiss-Prot full-mode arm reuses existing runs (`pytd_v7h_homology`, searched with 3 threads).
  Their outputs do not depend on the thread count.
- Hybrid and the selected top-25 variant were defined after seeing the UniRef50 results (hybrid) or tuned on
  human d1 + d8 (top 25). The held-out sets are as in step 4.
- The machine was heavily loaded (load 28-50 on 28 cores), so the run times are indicative.
- Footprint: 16.1 GB in all (UniRef50 db 13.9 GB, tables 1.1 GB, queries / searches / hits 0.5 GB,
  taxonomy 0.17 GB, full-mode runs 0.59 GB, `_trash` 0.23 GB).
