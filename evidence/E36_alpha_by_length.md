# E36. Per-transcript alpha in the start rule from transcript and ORF length, offline evaluation (2026-10-30)

Status: exploratory, not in the paper. No TransDecoder code was changed. Recommendation below; the decision is with Brian.

## Question
The release start rule trims a 5'-open ORF to its first in-frame ATG when 0.4|u| + 17K > C(u)/mu (E13, E35).
Complete transcripts want a high alpha and fragments want a low one (E13/S7). Brian's idea: bias alpha per transcript,
alpha_i = g(L, O), using only the transcript length L and the candidate ORF length O (or related quantities: O/L,
|u|/L). Does a simple monotone form with 1-2 extra parameters, tuned on the benchmark only, improve complete
and fragment inputs at once, and does that hold on assemblies it was not tuned on?

## Claims supported
None in the manuscript. Internal decision record:
1. The relation with **length** points the same way in the benchmark and the real-assembler data. Short transcripts
   and short ORFs are more often truly truncated, so they want a low alpha. The relation with **O/L** does not
   transfer. On the benchmark fragments a 3'-open ORF (O/L near 1) is always a UTR case (a 3' partial). On the
   assemblies it is truncated in 80-94% of cases. That reversal comes from the single-cut design.
2. Three forms tuned on the benchmark pass the E35 acceptance criteria ab initio and raise every held-out simulated
   assembly. The forms use ORF length, transcript length and |u|/L. The ORF-length form,
   alpha_i = 0.35 + 0.15 log2(O/1000), was chosen by the best leave-one-organism-out objective. It picks the same
   parameters in all 8 folds.
3. The forms tuned on O/L fail on the held-out assemblies, by up to -15 points ab initio and -6 with the assist. This
   happens even where they pass leave-one-organism-out on the benchmark.

## Data (all finished runs, read only; nothing rerun)
`$P = /broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T = $P/TransDecoderPaper`. Manifest `data/alpha_by_length/manifest.tsv`
(95 runs, `pipeline/alpha_by_length/00_manifest_v1.sh`):
- The 90 E35 runs, with E35's per-candidate tables (`$P/data/kozak_floor/candidates/`):
  - benchmark ab initio: 8 organisms x complete, complete_fullutr, fragments and fragments_mixed_t1-t3 (48);
  - benchmark with the assist: complete, complete_fullutr and fragments (24);
  - human d1, d8 and d32, ab initio and with the assist;
  - Arabidopsis d1-d32, A (ab initio) and B (full homology mode with the assist).
- New:
  - human d2, d4 and d16 ab initio (`$P/work/simtrinity/depth/<d>/runs/pytd_v7`, 3a793ff). These depths have no
    assist runs; the other agent's human depth homology runs (`$P/work/start_classifier/human_depth`) were not
    finished and were not used.
  - GM12878 Trinity GG ab initio (`final_v7_abinitio`, via E34's link dir `$P/data/missed_cds/gm12878_inputs/run`).
  - GM12878 Trinity GG with the assist (`final_v7h_homology`, a5b24c3; link dir `$P/data/alpha_by_length/gm12878_inputs/v7h`).
  - Truth for both GM12878 runs: the E29 transcript-level MANE truth (E34 copy with n_codons). It covers only contigs
    carrying a MANE Select stop, so GM12878 is scored only with the E29 statistic.
- Held out (never used for tuning): all simulated assemblies and GM12878.

## Software
PyTD a5b24c3 (`$P/work/trees/PY_hstart`, unchanged) through E35's `01_candidates_v1.py` (unchanged);
`$P/env/startclf` (Python 3.11.18, pandas 3.1.0, numpy 2.4.6, matplotlib). Paper repo at 64b58a4 plus the E36 files.

## Method (`pipeline/alpha_by_length/`)
- **01_candidates_new_v1.sh** runs E35's candidate extraction on the 5 new runs. Recomputed decisions = the runs' own
  decisions (0 mismatches; also 0 against the GM12878 homology table).
- **02_labels_v1.py** labels each candidate (5'-open ORF with an in-frame ATG):
  - L, O (the open ORF, including u), |u|, C(u)/mu, K, and the homology call;
  - the oracle answer: `keep` (keeping the 5' end open gives the exact CDS: a truncated CDS), `trim` (the first ATG is
    exact: a UTR), `neither`, `nostop` or `no_truth`;
  - the per-dataset V0 constants.
- **abl_common_v1.py** holds the fast rescoring. With everything else fixed, candidate i is trimmed iff
  alpha_i x |u| + 17K > C(u)/mu. With the assist, the homology call (trim / keep) overrides where it fires, as in PyTD.
  - Only starts change, so n_pred, the stop matches and the truth-typed prediction counts are constant.
  - tp_exact and the types of predictions on negative contigs are recomputed from per-candidate flags.
- **03_check_v1.py** checks the rescoring (`check_fast_scorer.tsv`):
  - V0 = E35's V0 on all 90 runs, within E35's 6-decimal rounding (max 5e-7);
  - E35's variants V1 and V3_8.5 recomputed = E35's full-substitution scores (90 x 2, max 5e-7);
  - an alpha-form variant scored both ways (fast vs full prediction substitution) on 10 runs including negatives and
    GM12878: 0 difference;
  - new runs vs official scores: human d2/d4/d16 0.8461 / 0.8456 / 0.8246 (= `simtrinity/depth/scores_overall.tsv`);
    GM12878 E29 statistic 0.7912 / 0.9266 (ab initio) and 0.9646 / 0.9441 (assist), = E33.
- **04_descriptive_v1.py** (step 1) bins the oracle answer by L, O, O/L, |u|/L and |u|:
  - ab initio runs; benchmark groups pooled over organisms, and per organism;
  - split by whether the candidate ends in a stop. A stop-ending candidate is complete-vs-5'-partial; one without a
    stop is 3'-partial-vs-internal;
  - per bin, the best alpha on a grid.
- **05_tune_grid_v1.py** (driver `05_tune_grid_all_v1.sh`, 3 lanes, nice 10) is step 2. It scores 7,220 configurations
  on every benchmark input. Forms (alpha clipped to [0, 5]):
  - const: a0;
  - logL: a0 + b log2(L/1000);
  - logO: a0 + b log2(O/1000);
  - OL: a0 + b (O/L - 0.8);
  - uL: a0 + b |u|/L;
  - stepOL: alo if O/L < r, else ahi;
  - stepL: alo if L < l, else ahi.
  - a0, alo, ahi range over 19 values from 0 to 3. The slopes and thresholds are on fixed grids (script header).
- **06_select_v1.py** does the leave-one-organism-out selection.
  - Objective per organism, in points of exact F1 vs V0:
    - ab initio: J = mean of complete, complete_fullutr, fragments and mixed (mean of t1-t3). This keeps E13's 50/50
      complete-vs-fragmented weighting, with the fragmented half shared by fragments and mixed;
    - assist: 0.25 complete + 0.25 complete_fullutr + 0.5 fragments (there are no mixed assist runs).
  - The strict E13 objective (8 complete + 8 fragments, equal weight) is reported.
  - Two selection modes:
    - `free`: argmax J;
    - `accept`: argmax J subject to the E35 criteria on the training organisms (mean complete delta >= -0.1, no
      training input < -0.5).
  - Each arm is tuned on its own benchmark inputs. The configuration selected on all 8 organisms is then applied in
    both arms.
- **07_evaluate_v1.py / 09_summarize_v1.py** (step 3) score every selected configuration on all 95 runs. Outputs:
  - overall and per-type F1 (with negatives);
  - changed decisions by class: w2r, r2w, w2w (neither exact), masked, nostop, no_truth;
  - the E29 statistic on GM12878.
- **08_ceiling_v1.py** (step 4) computes oracle ceilings chosen on each evaluation set itself, over
  alpha in {0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 1.5, 2, 3}:
  - best constant alpha;
  - best alpha per L bin (5 bins), per O/L bin (5 bins), and per (L, O/L) bin (25 bins; the requested ceiling);
  - every rule-decided candidate decided right.
  - No transcript has two candidates, so exact counts are additive per bin.
- **10_figure_v1.py** draws `truncated_share_by_length.png/.pdf`.

The selection criterion (best LOO J among the `accept` forms) was in view before the held-out scoring was run (06 before
07), but it was not written down beforehand.

## Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; S=$T/pipeline/alpha_by_length; D=$T/data/alpha_by_length
PY="nice -n 10 $P/env/startclf/bin/python"; L=$P/data/alpha_by_length/labels
bash $S/00_manifest_v1.sh
bash $S/01_candidates_new_v1.sh                                  # ~3 min
$PY $S/02_labels_v1.py $D/manifest.tsv $L                        # ~1 min
$PY $S/03_check_v1.py $D/manifest.tsv $L $P/data/kozak_floor/eval/scores.tsv $D/check_fast_scorer.tsv
$PY $S/04_descriptive_v1.py $D/manifest.tsv $L $D
bash $S/05_tune_grid_all_v1.sh                                   # ~10 min, -> $P/data/alpha_by_length/grid/
$PY $S/06_select_v1.py $P/data/alpha_by_length/grid $D
$PY $S/07_evaluate_v1.py $D/manifest.tsv $L $D/selected.tsv $P/data/alpha_by_length/eval
$PY $S/08_ceiling_v1.py $D/manifest.tsv $L $D
$PY $S/09_summarize_v1.py $P/data/alpha_by_length/eval $D
$PY $S/10_figure_v1.py $D/oracle_by_bin.tsv $D/truncated_share_by_length
```
Tables in `data/alpha_by_length/`:
- step 1: `oracle_by_bin.tsv`, `oracle_by_bin_org.tsv`, `oracle_2d_L_OL.tsv`, `feature_medians.tsv`, and the figure;
- step 2: `loo_folds.tsv`, `loo_summary.tsv`, `selected.tsv`, `best_by_form_insample.tsv`;
- step 3: `scores.tsv`, `bench_summary.tsv`, `heldout_summary.tsv`, `delta_by_type.tsv`, `changes_summary.tsv`;
- step 4: `ceiling.tsv`, `ceiling_bins.tsv`;
- `check_fast_scorer.tsv` and `manifest.tsv`.

Not in git (~60 MB): `$P/data/alpha_by_length/{candidates,labels,grid,eval}`.

## Results

### Step 1. Descriptive: is the truncated share related to L and O? (`feature_medians.tsv`, `oracle_by_bin.tsv`, figure)
Truncated share = keep-exact / (keep-exact + trim-exact) per bin. The candidates are split by their 3' end.

| group | stop-ending: L < 500 / 500-1000 / 1000-1500 / 2000-3000 / >= 5000 | no stop (3'-open) | O/L of truncated vs UTR (medians) |
|---|---|---|---|
| benchmark complete | 0 in every bin (no truncated CDSs) | - | - / 0.88 |
| benchmark fragments | 1.00 in every bin (every stop-ending fragment is a 5' cut) | 0.00 in every bin (every 3'-open fragment is a 3' partial) | 0.76 / 1.00 |
| benchmark mixed | 0.93 / 0.84 / 0.73 / 0.62 / 0.60 | 0.78 -> 0.66 | 0.98 / 0.98 |
| sim human d1 | 0.60 / 0.26 / 0.13 / 0.09 / 0.08 | 0.89 (L < 500) -> 0.5 | 0.99 / 0.61 |
| sim human d8 | 0.70 / 0.39 / 0.26 / 0.20 / 0.18 | 0.89 -> 0.6 | 0.99 / 0.73 |
| sim human d32 | 0.79 / 0.54 / 0.43 / 0.30 / 0.37 | 0.91 -> 0.6 | 0.99 / 0.82 |
| sim athal d1 | 0.70 / 0.16 / 0.13 / 0.12 / 0.13 | 0.93 -> 0.75 | 0.99 / 0.86 |
| sim athal d32 | 0.94 / 0.53 / 0.49 / 0.41 / - | 0.96 -> 0.8 | 0.99 / 0.87 |
| GM12878 Trinity (MANE) | 0.77 / 0.41 / 0.33 / 0.22 / 0.21 | (truth has no 3'-open CDSs) | 0.64 / 0.62 |

- **Where they agree.** On every real-assembler set (simulated human and Arabidopsis Trinity, all depths; GM12878
  Trinity) and on the benchmark mixed inputs, the truncated share falls monotonically with L and with O. This holds
  within stop-ending and within 3'-open candidates.
  - The best alpha per bin rises from 0 at L or O < 500 nt to 0.5-1.0 above 2,000 nt. In shallow athal the best alpha
    above 1,000 nt is up to 3.
  - On the benchmark the complete inputs (long, never truncated) and the fragments (cut, so shorter) give the same
    direction between inputs. So the sign of a length slope is shared by the benchmark and the assemblies.
- **O/L is where the synthetic cut creates a spurious relation.**
  - On the assemblies O/L carries almost nothing within stop-ending candidates: human about 0.11-0.30 in every O/L
    bin, GM12878 0.22-0.36. Its whole signal is the jump at O/L >= 0.95, i.e. "the ORF lacks a stop". There 73-94%
    are truncated (internal fragments) on the simulated assemblies.
  - In the benchmark fragments the same 3'-open class is 100% UTR. The single cut makes a 3'-open fragment always the
    5' piece, a 3' partial, and makes no internal fragments.
  - The mixed inputs show a strong O/L gradient within stop-ending candidates (0.98 at O/L < 0.5 to 0.59 at
    0.95-1). That gradient is absent on the assemblies, because there the cut position sets O/L.
- Shallow assemblies (d1) have long UTR-carrying contigs (median L 1,300-2,400) and short truncated ones (350-550 nt).
  Deeper ones have more short fragments, so the same length relation shifts toward "keep".

### Step 2. Benchmark tuning, leave-one-organism-out (`loo_summary.tsv`, `loo_folds.tsv`, `selected.tsv`)
Ab initio arm. Mean over the 8 held-out organisms of each fold, in points of exact F1 vs V0.

| form, mode | params chosen in the 8 folds | LOO J | LOO E13 | LOO complete | LOO fragments | LOO mixed | folds with J < 0 | worst held-out input |
|---|---|---|---|---|---|---|---|---|
| const, free (alpha retuned) | 0.35 (6 folds), 0.3 (2) | -0.03 | -0.07 | -0.58 | +0.44 | +0.58 | 2 | -1.75 |
| **logO, accept** | a0 0.35, b 0.15 in all 8 | **+0.43** | +0.43 | +0.08 | +0.78 | +0.84 | 0 | -0.45 |
| logO, free | a0 0.35 (6) / 0.3 (2), b 0.15 | +0.37 | +0.37 | -0.07 | +0.81 | +0.91 | 0 | -0.49 |
| logL, accept | a0 0.35, b 0.1 (7); a0 0.3, b 0.15 (athal fold) | +0.35 | +0.23 | +0.10 | +0.37 | +0.73 | 0 | -0.55 |
| uL, accept | a0 0.5, b -1 (7); a0 0.7, b -3.5 (tthermophila fold) | +0.33 | +0.15 | +0.18 | +0.11 | +0.75 | 1 (tthermophila -0.23) | -1.44 |
| stepL, accept | l 1000, alo 0.2, ahi 0.45 (6); l 1500 (2) | +0.27 | +0.18 | -0.02 | +0.37 | +0.68 | 0 | -0.50 |
| OL, free | b 0.25-1 | -0.05 | +0.36 | -0.63 | +1.34 | +0.35 | 4 | -3.80 |
| stepOL, free | r 0.8, alo 0.15-0.25 | +0.10 | +0.43 | -0.57 | +1.42 | +0.56 | 3 | -2.43 |

The OL and stepOL accept modes select V0 or near-V0 (the criteria bind). Assist arm, tuned on the 24 assist inputs:
- stepOL (r 0.95, alo 0.35, ahi 1-3) and OL (b 1.25-4) give the best LOO J (+0.40 to +0.53, 1 fold down);
- logO, logL, uL and stepL give +0.01 to +0.11.
The stepOL and OL winners raise alpha for 3'-open ORFs, i.e. they learn the benchmark's reversed O/L relation.

Configurations selected on all 8 organisms, ab initio benchmark (`bench_summary.tsv`):

| configuration | E13 objective (0.8477 at V0) | complete | complete_fullutr | fragments | mixed | inputs < -0.5 / 48 | worst input | E35 acceptance | both up |
|---|---|---|---|---|---|---|---|---|---|
| const 0.35 (free) | 0.8477 (-0.00) | -0.39 | -0.37 | +0.39 | +0.51 | 6 | tthermophila complete_fullutr -0.83 | no | no |
| **logO: 0.35 + 0.15 log2(O/1000)** | 0.8520 (+0.43) | +0.08 | +0.01 | +0.78 | +0.84 | 0 | athal complete_fullutr -0.45 | yes | yes |
| logL: 0.35 + 0.1 log2(L/1000) | 0.8500 (+0.23) | +0.14 | +0.25 | +0.32 | +0.64 | 0 | creinhardtii fragments -0.27 | yes | yes |
| uL: 0.5 - \|u\|/L | 0.8499 (+0.22) | +0.19 | +0.25 | +0.25 | +0.80 | 0 | tthermophila fragments -0.37 | yes | yes |
| stepL: 0.2 if L < 1000 else 0.45 | 0.8497 (+0.20) | -0.08 | -0.00 | +0.49 | +0.89 | 0 | tthermophila fragments -0.36 | yes | no |
| OL free: 0.3 + 0.5 (O/L - 0.8) | +0.38 | -0.56 | -1.01 | +1.32 | +0.60 | 10 | human complete_fullutr -2.58 | no | no |
| stepOL free: 0.15 if O/L < 0.8 else 0.4 | +0.52 | -0.49 | -0.95 | +1.54 | +0.60 | 11 | human complete_fullutr -2.43 | no | no |

- Per organism, logO is not uniform:
  - complete inputs: athal -0.44, dmel -0.20, spombe -0.08, mouse -0.02; creinhardtii +0.34, tthermophila +0.85;
  - fragments: creinhardtii -0.19, tthermophila -0.05, the rest +0.6 to +1.8.
- logL:
  - complete: 2 of 8 down (athal -0.21, dmel -0.04);
  - fragments: creinhardtii -0.27 and tthermophila -0.20.
- uL: complete 1 of 8 down; fragments 2 of 8 down.
- Applied in the assist arm (same parameters), on the benchmark:

| configuration | E13 (0.8823 at V0) | complete | fragments | inputs < -0.5 / 24 | acceptance | both up |
|---|---|---|---|---|---|---|
| logO | +0.04 | +0.06 | +0.02 | 0 | yes | yes (barely) |
| logL | +0.02 | +0.10 | -0.06 | 0 | yes | no |
| uL | +0.07 | +0.16 | -0.02 | 0 | yes | no |

### Step 3. Held-out evaluation (`heldout_summary.tsv`, `scores.tsv`, `delta_by_type.tsv`, `changes_summary.tsv`)
Change in overall exact F1, points vs V0. Benchmark-selected parameters (ab initio tuning) are applied in both arms.
Human d2/d4/d16 have ab initio runs only.

| dataset (V0) | const 0.35 | logO | logL | uL | stepL | OL free | stepOL free | stepOL accept, assist-tuned | ceiling (L x O/L) |
|---|---|---|---|---|---|---|---|---|---|
| human d1 (0.8440) | -0.10 | +0.31 | +0.54 | +0.26 | +0.29 | -1.40 | -1.29 | -1.60 | +1.21 |
| human d2 (0.8461) | -0.12 | +0.45 | +0.52 | +0.39 | +0.37 | -1.29 | -1.19 | -2.02 | +1.25 |
| human d4 (0.8456) | -0.15 | +0.39 | +0.42 | +0.43 | +0.39 | -1.02 | -1.06 | -2.99 | +1.03 |
| human d8 (0.8315) | +0.22 | +1.08 | +1.09 | +1.05 | +0.94 | -0.74 | -0.68 | -5.35 | +1.91 |
| human d16 (0.8246) | +0.43 | +1.51 | +1.29 | +1.36 | +1.33 | -0.34 | -0.36 | -9.02 | +2.29 |
| human d32 (0.8063) | +0.73 | +2.65 | +2.15 | +2.29 | +2.10 | -0.17 | -0.18 | -11.20 | +3.57 |
| athal d1 A (0.8815) | -0.20 | +0.39 | +0.43 | +0.82 | +0.33 | -0.51 | -0.46 | -3.16 | +3.35 |
| athal d2 A | -0.04 | +0.66 | +0.68 | +1.03 | +0.77 | -0.40 | -0.38 | -4.82 | +2.68 |
| athal d4 A | +0.18 | +1.45 | +1.16 | +1.55 | +1.20 | -0.24 | -0.22 | -7.32 | +2.86 |
| athal d8 A | +0.60 | +2.40 | +2.03 | +2.04 | +2.05 | -0.03 | -0.10 | -11.64 | +3.79 |
| athal d16 A | +0.91 | +3.43 | +2.78 | +2.92 | +2.73 | +0.18 | +0.05 | -15.37 | +4.90 |
| athal d32 A (0.8281) | +0.89 | +3.93 | +3.07 | +3.05 | +2.89 | +0.27 | +0.11 | -15.61 | +4.95 |
| human d1 + assist (0.8641) | -0.03 | +0.07 | +0.11 | +0.05 | +0.04 | -0.22 | -0.21 | -0.17 | +0.32 |
| human d8 + assist (0.8859) | -0.05 | +0.13 | +0.17 | +0.12 | +0.14 | -0.22 | -0.17 | -0.47 | +0.48 |
| human d32 + assist (0.8963) | +0.05 | +0.22 | +0.17 | +0.18 | +0.16 | -0.02 | -0.05 | -0.54 | +0.40 |
| athal d1 B (0.9138) | -0.30 | -0.03 | +0.11 | +0.32 | -0.04 | -0.52 | -0.45 | -1.02 | +3.11 |
| athal d2 B | -0.19 | +0.05 | +0.18 | +0.39 | +0.17 | -0.43 | -0.36 | -1.55 | +2.58 |
| athal d4 B | -0.10 | +0.38 | +0.38 | +0.63 | +0.32 | -0.28 | -0.22 | -2.84 | +2.31 |
| athal d8 B | +0.16 | +1.02 | +0.88 | +0.87 | +0.89 | -0.16 | -0.18 | -4.81 | +2.02 |
| athal d16 B | +0.32 | +1.35 | +1.11 | +1.23 | +1.04 | -0.05 | -0.10 | -6.20 | +2.14 |
| athal d32 B (0.9083) | +0.30 | +1.57 | +1.24 | +1.30 | +1.21 | +0.11 | +0.06 | -5.22 | +1.96 |

(The stepOL-accept column uses the assist-tuned parameters r 0.95, alpha 0.35 / 1, in both arms. The OL accept form
tuned on the assist inputs (0.45 + 1.25 (O/L - 0.8)) gives -1.2 to -7.2 ab initio and -0.1 to -2.7 with the assist.)

**GM12878 Trinity, E29 statistic** (MANE references whose stop is predicted; percent; delta vs V0):

| configuration | ab initio: 5' partials kept open given stop (79.12) | ab initio: complete exact given stop (92.66) | assist: 5' partials (96.46) | assist: complete (94.41) |
|---|---|---|---|---|
| const 0.35 | +1.77 | -0.42 | +0.14 | -0.04 |
| logO | **+2.02** | -0.10 | +0.20 | +0.01 |
| logL | -1.94 | +0.35 | +0.14 | +0.07 |
| uL | -0.50 | +0.16 | -0.14 | +0.04 |
| stepL | -0.76 | +0.18 | +0.14 | +0.08 |
| OL free | +7.24 | -1.87 | +0.68 | -0.20 |
| stepOL free | +6.06 | -1.78 | +0.61 | -0.20 |

The GM12878 net change in MANE exact CDSs is small for every length form: logO +13 (80 w2r, 67 r2w); logL +17; uL +12.

**Per type** (simulated assemblies, f1 with negatives, logO ab initio):
- internal: +3.1 to +6.7;
- 5' partial: +0.7 to +4.4;
- 3' partial: -0.9 to -2.4;
- complete: -0.6 to +0.1.

With the assist, logO gives internal +0.2 to +3.0, 5' partial +0.1 to +1.7, 3' partial -0.1 to -2.3, complete -0.6 to
+0.1. logL and uL keep the assemblies' complete type near 0 or up (uL +0.1 to +0.5) and lose less on 3' partials
(uL -0.0 to -1.1). Their gains are again mostly internal and 5'-partial CDSs. The gain on the assemblies therefore
comes from keeping short, mostly 3'-open fragments open. Complete CDSs on long contigs gain only with logL and uL.

**Decisions changed** (logO; w2r = wrong -> right, r2w = right -> wrong):

| group | changed | w2r | r2w | net exact | neither exact (w2w) |
|---|---|---|---|---|---|
| bench complete, ab initio | 1,503 | 755 | 614 | +141 | 98 |
| bench complete_fullutr | 1,288 | 621 | 580 | +41 | 59 |
| bench fragments | 3,354 | 2,071 | 1,183 | +888 | 52 |
| bench mixed (24) | 11,175 | 7,022 | 3,925 | +3,097 | 149 |
| sim human, 6 depths | 2,663 | 1,779 | 646 | +1,133 | 30 |
| sim athal, 6 depths | 3,659 | 2,894 | 724 | +2,170 | 3 |
| GM12878 (MANE-stop contigs) | 1,204 (1,050 without truth) | 80 | 67 | +13 | 7 |
| assist: sim human (3) | 222 | 112 | 35 | +77 | 15 |
| assist: sim athal (6) | 1,873 | 1,303 | 536 | +767 | 3 |
| assist: bench complete / fragments | 1,126 / 1,705 | 567 / 780 | 445 / 837 | +122 / -57 | 88 / 51 |

Unlike a constant alpha or the Kozak floor (E35), a length-dependent alpha moves decisions both ways on complete
inputs. It trims more on long ORFs (w2r) and keeps more short ORFs open (r2w), so the complete inputs' net is small.

### Step 4. Ceiling: best alpha per bin chosen on each evaluation set itself (`ceiling.tsv`, `ceiling_bins.tsv`)
Points of exact F1 vs V0 (table above, last column; per-dataset rows in `ceiling.tsv`).

| dataset | best constant alpha (value) | L bins | O/L bins | L x O/L bins | every candidate right | logO achieves |
|---|---|---|---|---|---|---|
| human d1 / d8 / d32 | +0.30 (0.6) / +0.35 (0.3) / +2.45 (0) | +1.04 / +1.46 / +3.31 | +0.95 / +1.44 / +3.20 | +1.21 / +1.91 / +3.57 | +3.6 / +7.6 / +11.2 | +0.31 / +1.08 / +2.65 |
| athal d1 / d8 / d32 A | 0 (0.4) / +2.36 (0) / +4.56 (0) | +2.12 / +3.10 / +4.78 | +2.64 / +3.32 / +4.73 | +3.35 / +3.79 / +4.95 | +8.7 / +12.9 / +14.2 | +0.39 / +2.40 / +3.93 |
| human d1 / d8 / d32 + assist | +0.13 / 0 / +0.12 | +0.25 / +0.36 / +0.34 | | +0.32 / +0.48 / +0.40 | +0.8 / +1.4 / +1.6 | +0.07 / +0.13 / +0.22 |
| athal d1 / d8 / d32 B | +0.97 (1.5) / +0.33 / +1.45 | +2.63 / +1.38 / +1.73 | | +3.11 / +2.02 / +1.96 | +6.1 / +7.7 / +7.4 | -0.03 / +1.02 / +1.57 |
| GM12878 (F1 not meaningful) | | | | ab initio +0.10, assist +0.05 | | |

- The ceiling is chosen in-sample and per dataset. It is an upper bound for any rule that sets alpha from length bins
  without knowing the input's truncation rate.
- At depth (d16-d32), most of the ceiling is a lower constant alpha: the best constant is 0 there. That is the input
  mix, not a length relation.
- At d1-d8 a constant cannot help human (0-0.35), but length bins can (+1.0 to +1.9). logO recovers about a quarter of
  that at d1 and about half at d8.
- Shallow Arabidopsis has the largest length-specific headroom (+3.4 at d1 ab initio; +3.1 with the assist). logO
  takes only +0.39 / -0.03 of it. The athal per-bin optimum (alpha 3 above 1,000 nt) is far from the benchmark's slope.
- Every length rule stays far below the per-candidate oracle (+3.6 to +14.7 ab initio). The keep/trim error is mostly
  not a length effect (E33, E34).
- On the benchmark fragments the O/L-bin ceiling (+4.6) is much larger than the L-bin ceiling (+2.0). That is the
  spurious 3'-open relation of step 1; on the assemblies O/L adds little over L.

## Interpretation
- The length relation is real and has the same sign in the benchmark and in real-assembler output: short transcripts
  and ORFs are mostly fragments. A per-transcript alpha that falls with ORF or transcript length captures part of it.
  A logarithmic ORF-length or transcript-length slope tuned on the benchmark improves every held-out simulated
  assembly, ab initio (+0.3 to +3.9) and with the assist (-0.03 to +1.6). It keeps the benchmark within the E35
  acceptance criteria. By contrast, the UTR prior (E29/S10) and every Kozak floor (E35) raised complete inputs only by
  lowering fragments and assemblies.
- It improves "complete and fragment inputs at once" only on average. The complete-input gain is small (+0.08 logO,
  +0.14 logL, +0.19 uL) and splits by organism (athal and dmel down). Most of the gain is on fragmented inputs.
  Within the assemblies the gain is in internal and 5'-partial CDSs, with a cost on 3' partials. The effect is a
  smarter alpha for short fragments rather than a higher alpha for complete transcripts.
- O/L must not be used. Its benchmark signal is created by the single-cut fragment design and is reversed on
  assemblies. An O/L rule tuned on the assist benchmark passes leave-one-organism-out and then loses up to 15 points
  on held-out assemblies. Leave-one-organism-out does not protect against a design artefact shared by all organisms.
- GM12878 (real Trinity) neither confirms nor contradicts:
  - logO keeps 2 points more true 5' partials open at a 0.1-point cost on complete references;
  - logL does the reverse (-1.9 / +0.35);
  - with the assist, all changes are <= 0.2 points.

## Recommendation
- logO, alpha_i = max(0, 0.35 + 0.15 log2(O/1000)) with O the open candidate ORF length in nt, clearly passes:
  - it was selected by the leave-one-organism-out criterion fixed before the held-out scoring was run, with identical
    parameters in all 8 folds;
  - it meets the E35 acceptance criteria ab initio (complete +0.08, worst input -0.45) and with the assist;
  - it raises the E13 objective (+0.43 ab initio, +0.04 assist);
  - it raises every held-out simulated assembly ab initio (+0.31 to +3.93). With the assist it raises 8 of 9 (athal d1
    B -0.03, within noise);
  - on GM12878 it gives +2.0 points of 5' partials kept open.
  - In effect alpha is about 0.09 at O = 300 nt, 0.35 at 1 kb and 0.59 at 3 kb.
- Alternatives that also pass:
  - logL, 0.35 + 0.1 log2(L/1000): uses transcript length; better for complete references, worse for GM12878 5' partials;
  - uL, 0.5 - |u|/L: best balance across types and the largest shallow-athal gain, but one LOO fold down and a worst
    LOO held-out input of -1.44.
- If adopted, implement it in both PyTD and Perl, in ab initio and assist modes. Re-run the benchmark and assemblies
  with the program, because the offline evaluation substitutes decisions only, and run a Perl/PyTD parity check.
- The decision is Brian's. The complete-transcript gain is small, so do not present it as solving the
  complete-vs-fragment trade-off. It shifts alpha toward "keep" for short ORFs.

## Guard against overfitting
- Tuning used the benchmark only; the assemblies and GM12878 were scored once per selected configuration.
- 7,220 configurations were searched on the benchmark. LOO estimates out-of-organism performance, but not artefacts
  shared by all benchmark inputs (the O/L case shows this).
- The forms and the acceptance criteria came from the task. The step-1 descriptive tables (which include the
  assemblies) were seen before tuning, but they were not used to select parameters.
- The ceilings are in-sample by construction.

## Caveats
- The evaluation is an offline substitution on the runs' own selections. It is exact for the start rule (checks), but
  the program was not run with the new rule.
- Single simulation per depth; human depths share genes; differences of ~0.05 points are not interpretable.
- Human d2/d4/d16 have no assist runs. The other agent's human depth homology runs were unfinished and were not used.
- The GM12878 MANE truth covers only MANE-stop contigs, with no internal fragments and no non-MANE truncated ORFs.
  Its F1 is not meaningful; only the E29 shares are reported.
- The benchmark fragments are cut once per CDS, so they have no internal fragments and a reversed 3'-open class.
  Mixed inputs are balanced by design. Neither reproduces assembly length distributions exactly.

## Addendum (2026-10-30): simpler forms (`data/alpha_by_length/simpler_v2/`)
Question (BH): is there a simpler, more intuitive formula that works nearly as well as logO?

**Method.**
- New scripts: `abl_common_v2.py` (adds forms), `05_tune_grid_v2.py`, `05_tune_grid_all_v2.sh`, `07_evaluate_v2.py`.
  06 and 09 v1 were reused unchanged.
- Raw grid in `$P/data/alpha_by_length/grid_v2/`, evaluation in `.../eval_v2/` (not in git).
- Same objective, LOO procedure, acceptance criteria and held-out sets as steps 2-3.
- Forms tested:
  - stepO: alo if O < o else ahi;
  - rampO: a min(1, O/o);
  - linO: c O/1000 (one parameter);
  - const and logO, rescored for reference (they reproduce the v1 values).

**Results, ab initio** (points of exact F1 vs V0):

| form (selected on all 8, accept mode) | LOO J | params across 8 folds | complete | fragments | mixed | human d1/d8/d32 | athal d1/d8/d32 | GM12878 5' partials kept open |
|---|---|---|---|---|---|---|---|---|
| logO 0.35 + 0.15 log2(O/1000) | +0.43 | identical in 8 | +0.08 | +0.78 | +0.84 | +0.31 / +1.08 / +2.65 | +0.39 / +2.40 / +3.93 | +2.0 |
| rampO 0.5 min(1, O/1200) | +0.37 (1 fold -0.04) | o 1200 (6) / 1500 (2); a 0.5 (6) / 0.6 (2) | +0.15 | +0.64 | +0.70 | +0.36 / +1.02 / +2.33 | +0.43 / +2.21 / +3.39 | +0.25 |
| stepO 0.15 if O < 900 else 0.5 | +0.29 | o 750-1000; 5 distinct | -0.08 | +0.70 | +0.82 | +0.18 / +0.87 / +2.33 | +0.26 / +2.16 / +3.29 | +1.0 |
| linO 0.225 O/1000 | +0.12 (2 folds down) | 4 distinct; accept mode: none | -0.42 | +0.94 | +0.87 | -0.11 / +0.94 / +2.77 | -0.04 / +2.61 / +4.22 | +4.0 |

- linO fails the acceptance criteria: 13 of 48 inputs are below -0.5 (athal complete -1.34). Long ORFs get too
  large an alpha.
- With the assist, all forms tuned ab initio give +0.0 to +0.2 points on human and up to +1.6 on athal d32, as logO did.

**The selected logO formula rewritten.**
- 0.35 + 0.15 log2(O/1000) = 0.15 log2(O/199.6) is numerically equivalent to **alpha = 0.5 log10(O / 200 nt)**:
  - alpha is 0 for a 200-nt open ORF and rises by 0.5 for each tenfold increase in length;
  - it equals the release constant (0.35-0.4) at about 1 kb;
  - it can be stated as "+0.15 per doubling of ORF length".
- Scored on all 95 runs (817,983 candidates), the two forms differ in 140 decisions, a net +10 exact CDSs.

**Conclusion.**
- None of the simpler forms is as good as logO. In particular, a proportional ramp or a two-value step recovers
  about 70-90% of its benchmark and held-out gain, with less stable parameters across folds.
- The ramp has the better complete-input behaviour (+0.15, no organism below -0.5) but keeps far fewer true
  5' partials open on GM12878 (+0.25 vs +2.0).
- Recommended statement of the rule, if adopted: alpha = 0.5 log10(O/200 nt), floored at 0.

## Decision (2026-10-30)
Deferred. Brian: "Let's note it for a possible future addition. I'm not sure it's worth the effort given the 0.4 pt
improvement." No code change in Perl v7.0.0 or PyTD v0.2.0, and the manuscript is unchanged. If it is revisited:
- use alpha = 0.5 log10(O/200 nt), floored at 0;
- implement it in both implementations (ab initio and with the assist);
- check Perl/PyTD parity;
- rerun the benchmark, GM12878 and the simulated assemblies with the program.
