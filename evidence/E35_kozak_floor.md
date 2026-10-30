# E35. One-sided Kozak term in the start rule ("Kozak floor"), offline evaluation (2026-10-30)

Status: exploratory, NOT ADOPTED (recommendation; decision with Brian). Not in the paper. No TransDecoder code was changed.

## Question
E34 (part 2) found that 19-36% (human) and 33-42% (Arabidopsis) of the UTR extensions wrongly kept 5' open are
"Kozak vetoes": u looks like UTR (coding(u) about 0) but the first ATG has K < 0, so
UTR(u) = 0.4|u| + 17K falls below coding(u) = C(u)/mu (e.g. athal d1: 0.4 x 28.5 + 17 x (-1.39) = -10.5).
Brian approved testing a one-sided Kozak term. Does limiting the penalty for a poor context improve exact-CDS
accuracy without hurting complete inputs, ab initio and with the homology start assist?

## Claims supported
None in the manuscript. Internal decision record: no one-sided Kozak variant meets the pre-set acceptance criteria
in ab initio mode, and none is clearly better with the homology assist; the release rule (w = 17 per bit for both
signs) stays.

## Variants (fixed before evaluation; only the Kozak term changes)
Trim to the first in-frame ATG when 0.4|u| + T(K) > coding(u).

| variant | T(K) |
|---|---|
| V0 (release, control) | 17 K |
| V1 | 17 max(K, 0) |
| V2_05 / V2_1 | 17 max(K, -0.5) / 17 max(K, -1) |
| V3_8.5 / V3_4 | 17 K if K >= 0, else 8.5 K / 4 K |
| V4 | 17 max(K, 0) if coding(u) <= 0, else 17 K |

Every variant can only raise the UTR score, so it can only turn "keep open" into "trim". With homology hits the
assist's call (R2 trim; R1/R3 keep; PyTD a5b24c3 `homology_start_call`) overrides the ab initio rule where a rule
fires; conflicts and no call fall back to the (variant) ab initio rule, as in the program.

Acceptance criteria (set in the task before any result): mean exact F1 on complete inputs not lower by more than
0.1 point, and no benchmark input lower by more than 0.5 point. Also reported: the E13 objective (equal-weighted
mean exact F1 of the 8 complete and 8 fragments stranded benchmark inputs).

## Data (all finished runs, read only; nothing rerun or resumed)
`$P = /broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T = $P/TransDecoderPaper`. Manifest: `data/kozak_floor/manifest.tsv`
(90 runs; written by `pipeline/kozak_floor/00_manifest_v1.sh`, which also checks that every path exists).

| arm | datasets | runs | truth |
|---|---|---|---|
| ab initio | benchmark (E07/E28), stranded, 8 organisms x complete, complete_fullutr, fragments, fragments_mixed_t1-t3 (48) | `$P/work/bench_v3/runs/<org>/<input>/pytd_v7` (3a793ff) | `$P/data/bench_v3/inputs/<org>/<input>.truth.tsv` |
| ab initio | simulated human Trinity d1, d8, d32 (E30/E31) | `$P/work/simtrinity/runs/pytd_v7`, `.../depth/<d>/runs/pytd_v7` (3a793ff) | `$P/data/simtrinity[/depth/<d>]/labels/sim_trinity.{truth.tsv,excluded.txt,negatives.txt}` |
| ab initio | simulated Arabidopsis d1-d32, condition A (E32/E33) | `$P/work/start_classifier/athal_depth/A/depth/<d>/runs/pytd_v7` (a5b24c3) | `$P/data/simtrinity_athal/depth/<d>/labels/sim_trinity.*` |
| homology + assist | benchmark complete, complete_fullutr, fragments x 8 (24) | `.../runs/<org>/<input>/pytd_v7h_homology` (a5b24c3; E33 adoption) | as above |
| homology + assist | human d1, d8, d32 | `$P/work/start_classifier/hstart_runs/sim.human.<d>/on` (a5b24c3) | as above |
| homology + assist | Arabidopsis d1-d32, condition B (full homology mode) | `.../athal_depth/B/depth/<d>/runs/pytd_v7` (a5b24c3) | as above |

The benchmark mixed inputs have no homology runs (as in E33), so they appear in the ab initio arm only.

## Software
- PyTD a5b24c3 (`$P/work/trees/PY_hstart`, unchanged), imported for frame contrast, mu, Kozak scoring and
  `homology_start_call`. Its ab initio output equals 3a793ff (E33).
- `$P/env/startclf` (Python 3.11.18, pandas 3.1.0, numpy 2.4.6).
- Paper repo at 7a423b8 + the E35 files.

## Method (scripts `pipeline/kozak_floor/`)
1. **`01_candidates_v1.py`** (driver `01_candidates_all_v1.sh`, 4 lanes, `nice -n 10`; ~3 min): for every ORF in the
   run's `best_candidates.gff3` that does not begin with ATG and has a first in-frame ATG, recompute exactly as
   `_revise_start_sites`: |u|, coding(u) = C(u)/mu (training-ORF median per-base contrast), K (bits; the run's
   `start_refinement.kozak.tsv`), the ab initio decision; for homology runs the run's own BLAST qstart/sstart and Pfam
   env_from (`start_refinement.homology.tsv`) through `homology_start_call`. Floats written at full precision.
   Check: recomputed release decision = the run's decision (`start_revised` in `best_candidates.gff3.revised_starts.gff3`)
   and = the homology table's decisions.
2. **`02_evaluate_v1.py`**: for each run and variant, the run's `preds.tsv` with each candidate's decision substituted
   (stops never change; excluded contigs removed), scored as `score_predictions.py` / `16_score_adoption.py`: stop within
   2 nt, exact = start also within 2 nt, one credit per truth CDS; overall F1 = 2 tp / (n_truth + all scored predictions,
   incl. negatives); per type f1 = 2 tp / (n of type + predictions on transcripts of that type); per type "with negatives"
   (f1w, as `13b_depth_collect.py`) also counts predictions on negative contigs typed by their own ends. Each changed
   decision is classed by the transcript's truth CDS: w2r (exact after, not before), r2w, w2w_later (stop matched, true
   start downstream of the first ATG: trimmed to a different ATG), w2w_upstream, nostop_truth (candidate stop is not the
   truth stop), no_truth / no_truth_neg.
3. **`03_check_official_v1.py`**: V0 vs the official score tables.
4. **`04_summarize_v1.py`**, **`05_changed_profile_v1.py`**: the tables below.

## Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; S=$T/pipeline/kozak_floor; PY=$P/env/startclf/bin/python
bash $S/00_manifest_v1.sh                                   # -> data/kozak_floor/manifest.tsv (90 runs)
nohup bash $S/01_candidates_all_v1.sh &                     # -> $P/data/kozak_floor/candidates/<dataset>.tsv.gz, extract.log
nice -n 10 $PY $S/02_evaluate_v1.py $T/data/kozak_floor/manifest.tsv $P/data/kozak_floor/candidates $P/data/kozak_floor/eval
$PY $S/03_check_official_v1.py $P/data/kozak_floor/eval/scores.tsv $T/data/kozak_floor/check_v0_vs_official.tsv
$PY $S/04_summarize_v1.py $P/data/kozak_floor/eval $T/data/kozak_floor
$PY $S/05_changed_profile_v1.py $P/data/kozak_floor/eval/changed_candidates.tsv.gz $T/data/kozak_floor/changed_profile.tsv
```
Outputs in `data/kozak_floor/`: `scores.tsv` (absolute, dataset x variant), `delta_overall.tsv`, `delta_by_type.tsv`,
`group_summary.tsv`, `e13_objective.tsv`, `acceptance.tsv`, `changes_summary.tsv`, `changes_by_type.tsv`,
`changed_profile.tsv`, `kozak_veto.tsv`, `check_v0_vs_official.tsv`, `manifest.tsv`. Not in git (46 MB):
`$P/data/kozak_floor/{candidates,eval}/` (per-candidate tables, `changed_candidates.tsv.gz`, logs).

## Checks
- Recomputed release decisions equal the runs' own decisions for every candidate of all 90 runs (0 mismatches;
  `candidates/extract.log`), and equal the homology tables' ab initio and final decisions.
- V0 predictions equal each run's `preds.tsv` (90/90; `checks.tsv`).
- V0 reproduces the official scores: 328 comparisons over 90 runs, 0 differences beyond the official's rounding
  (`check_v0_vs_official.tsv`: benchmark `adoption/bench_scores.tsv` (5 dp) and `benchmark/final_v7/mixed_scores.tsv`;
  human `adoption/sim_scores.tsv` and `simtrinity/depth/scores_by_type.tsv` (with negatives); Arabidopsis A/B
  `adoption/athal_depth_scores.tsv` incl. per type with negatives). E.g. human d1 0.8440 / 0.8641 (H), athal d1
  0.8815 / 0.9138 (B), human complete 0.91993 / 0.95912 (H).
- Kozak-veto candidates (kept open by the ab initio rule only because K < 0) on the simulated assemblies agree with
  E34 (e.g. the V1 w2r count on human d1+d8+d32, 452, vs E34's 440 vetoed complete and 3'-partial UTRs).

## Results

### E13 objective and acceptance (`e13_objective.tsv`, `acceptance.tsv`)
Change in exact F1, points (100 x difference) vs V0. "complete" = the 8 stranded complete inputs (with complete_fullutr
16 inputs: within 0.04 of it). Worst = the lowest benchmark input (48 ab initio, 24 homology).

| arm | variant | mean d complete | mean d fragments | E13 objective (abs; d) | bench inputs < -0.5 | worst input | acceptable |
|---|---|---|---|---|---|---|---|
| ab initio | V0 | 0 | 0 | 0.8477 | 0 | - | (control) |
| | V1 | +1.41 | -3.59 | 0.8369; -1.09 | 32 / 48 | mouse fragments -5.08 | no |
| | V2_05 | +0.74 | -1.64 | 0.8433; -0.45 | 31 | mouse fragments -2.51 | no |
| | V2_1 | +0.34 | -0.80 | 0.8455; -0.23 | 21 | human fragments -1.36 | no |
| | V3_8.5 | +0.61 | -0.85 | 0.8465; -0.12 | 26 | mouse mixed_t3 -1.31 | no |
| | V3_4 | +1.01 | -1.75 | 0.8440; -0.37 | 32 | mouse fragments -2.52 | no |
| | V4 | +0.53 | -1.53 | 0.8427; -0.50 | 31 | mouse fragments -2.37 | no |
| homology + assist | V0 | 0 | 0 | 0.8823 | 0 | - | (control) |
| | V1 | +1.13 | -1.28 | 0.8816; -0.08 | 7 / 24 | pfalciparum fragments -2.48 | no |
| | V2_05 | +0.58 | -0.44 | 0.8830; +0.07 | 3 | pfalciparum fragments -1.01 | no |
| | V2_1 | +0.27 | -0.16 | 0.8829; +0.06 | 0 | pfalciparum fragments -0.41 | yes (benchmark criteria) |
| | V3_8.5 | +0.48 | -0.17 | 0.8839; +0.15 | 0 | pfalciparum fragments -0.43 | yes (benchmark criteria) |
| | V3_4 | +0.80 | -0.49 | 0.8839; +0.16 | 4 | pfalciparum fragments -1.07 | no |
| | V4 | +0.44 | -0.65 | 0.8813; -0.10 | 6 | creinhardtii fragments -1.16 | no |

- Every variant raises every complete and complete_fullutr input (ab initio +0.06 to +2.56; largest in spombe and
  athal, smallest in pfalciparum) and lowers every fragments and fragments_mixed input in ab initio mode
  (fragments -0.31 to -5.08, mixed -0.17 to -4.57; `group_summary.tsv`, `delta_overall.tsv`).
- Ab initio, the mean loss on fragments is 1.4-2.9 times the mean gain on complete inputs; the E13 objective falls for every
  variant (-0.12 to -1.09). The E13 trade-off is reproduced: a floor acts like a larger alpha for short u.
- With the assist, R1/R3 already keep many truncated ORFs open, so the fragment loss shrinks; V2_1 and V3_8.5 pass
  the benchmark criteria (and raise the E13 objective by +0.06 / +0.15), but see the assemblies below.

### Simulated assemblies, overall exact F1 change (points; `delta_overall.tsv`)

| dataset (V0 exact F1) | V1 | V2_05 | V2_1 | V3_8.5 | V3_4 | V4 |
|---|---|---|---|---|---|---|
| human d1 ab initio (0.8440) | -0.33 | -0.07 | -0.04 | +0.04 | -0.07 | -0.26 |
| human d8 (0.8315) | -2.10 | -1.22 | -0.72 | -0.47 | -1.14 | -0.91 |
| human d32 (0.8063) | -6.11 | -3.23 | -1.87 | -1.31 | -2.93 | -2.45 |
| athal d1 A (0.8815) | +0.20 | +0.03 | -0.09 | +0.06 | +0.21 | +0.17 |
| athal d2 A | -0.52 | -0.34 | -0.21 | -0.13 | -0.24 | -0.11 |
| athal d4 A | -1.69 | -0.80 | -0.46 | -0.49 | -0.90 | -0.58 |
| athal d8 A | -3.75 | -1.98 | -1.03 | -1.02 | -2.11 | -1.25 |
| athal d16 A | -5.78 | -3.00 | -1.67 | -1.67 | -3.24 | -2.12 |
| athal d32 A (0.8281) | -8.79 | -4.78 | -2.52 | -2.13 | -4.52 | -3.48 |
| human d1 + assist (0.8641) | -0.12 | -0.02 | -0.00 | -0.01 | -0.01 | -0.17 |
| human d8 + assist (0.8859) | -0.39 | -0.11 | -0.05 | +0.02 | -0.10 | -0.34 |
| human d32 + assist (0.8963) | -1.15 | -0.36 | -0.12 | -0.08 | -0.25 | -0.89 |
| athal d1 B (0.9138) | +0.86 | +0.41 | +0.12 | +0.31 | +0.60 | +0.41 |
| athal d2 B | +0.46 | +0.15 | +0.02 | +0.17 | +0.34 | +0.17 |
| athal d4 B | -0.23 | -0.03 | -0.06 | -0.00 | -0.04 | -0.08 |
| athal d8 B | -1.35 | -0.61 | -0.27 | -0.23 | -0.67 | -0.55 |
| athal d16 B | -2.40 | -1.11 | -0.58 | -0.57 | -1.24 | -1.10 |
| athal d32 B (0.9083) | -3.56 | -1.71 | -0.81 | -0.66 | -1.57 | -1.50 |

Stop F1 is unchanged everywhere (decisions move only starts). The sign follows the truncated share: gains only where
almost no CDS is 5'-truncated (complete inputs; athal d1-d2; human d1 marginally), losses growing with depth
(more fragments).

Per type with negatives (`delta_by_type.tsv`, column d_f1w), e.g. V3_8.5 / V2_1: complete +0.3..+1.4 / +0.2..+1.1 and
3' partial +0.2..+2.4 / +0.1..+1.1 on every ab initio assembly, but 5' partial -2.3..-4.4 / -2.1..-4.1 and internal
-1.2..-4.5 / -1.4..-4.0. With the assist (V3_8.5): human 5' partial -0.03..-0.59, internal -0.24..-0.37; athal
5' partial -1.4..-2.0, internal -1.4..-2.0. V1 costs 5' partials and internal CDSs 6-14 points ab initio.

### Decisions changed (`changes_summary.tsv`, `changes_by_type.tsv`)
Totals over each group; w2r = wrong -> right (exact gained), r2w = right -> wrong (exact lost), w2w_later = trimmed to
the first ATG while the true start is a later ATG (exact neither way).

| arm, group | variant | changed | w2r | r2w | w2w_later | net exact |
|---|---|---|---|---|---|---|
| ab initio, bench complete (8) | V1 / V2_1 / V3_8.5 / V4 | 2,416 / 639 / 999 / 981 | 1,959 / 494 / 841 / 701 | 0 / 0 / 0 / 0 | 439 / 142 / 153 / 275 | +1,959 / +494 / +841 / +701 |
| ab initio, bench complete_fullutr (8) | same | 2,239 / 594 / 913 / 939 | 1,845 / 466 / 782 / 689 | 0 | 377 / 124 / 127 / 246 | +1,845 / +466 / +782 / +689 |
| ab initio, bench fragments (8) | same | 6,653 / 1,682 / 1,998 / 2,847 | 947 / 277 / 414 / 351 | 5,439 / 1,321 / 1,492 / 2,351 | 217 / 68 / 75 / 132 | -4,492 / -1,044 / -1,078 / -2,000 |
| ab initio, bench mixed (24) | same | 19,117 / 4,399 / 6,061 / 7,823 | 2,654 / 696 / 1,180 / 923 | 15,742 / 3,454 / 4,625 / 6,439 | 649 / 228 / 235 / 435 | -13,088 / -2,758 / -3,445 / -5,516 |
| ab initio, sim human (3) | same | 2,465 / 760 / 759 / 993 | 452 / 143 / 221 / 141 | 1,777 / 542 / 472 / 729 | 116 / 44 / 36 / 76 | -1,325 / -399 / -251 / -588 |
| ab initio, sim athal (6) | same | 6,635 / 1,702 / 2,159 / 2,741 | 1,573 / 334 / 610 / 745 | 5,037 / 1,357 / 1,534 / 1,980 | 10 / 4 / 5 / 5 | -3,464 / -1,023 / -924 / -1,235 |
| assist, bench complete (8) | same | 2,041 / 542 / 826 / 886 | 1,596 / 396 / 671 / 608 | 0 | 430 / 142 / 149 / 274 | +1,596 / +396 / +671 / +608 |
| assist, bench fragments (8) | same | 3,309 / 687 / 960 / 1,573 | 770 / 208 / 336 / 301 | 2,298 / 405 / 544 / 1,137 | 213 / 66 / 70 / 129 | -1,528 / -197 / -208 / -836 |
| assist, sim human (3) | same | 651 / 134 / 143 / 399 | 104 / 26 / 44 / 31 | 372 / 53 / 53 / 270 | 110 / 40 / 31 / 71 | -268 / -27 / -9 / -239 |
| assist, sim athal (6) | same | 3,932 / 884 / 1,251 / 1,857 | 1,449 / 305 / 542 / 704 | 2,459 / 566 / 693 / 1,136 | 10 / 4 / 5 / 5 | -1,010 / -261 / -151 / -432 |

The remaining changes are candidates whose stop is not the truth stop (0-68 per group) or negative contigs on the
assemblies (3-118 per group; these change only the per-type-with-negatives denominators).

**Why the Kozak veto is mostly right** (`changed_profile.tsv`). Among V1's changed decisions with a stop-matched truth
CDS, the release keep-open was right (truly truncated CDS whose first in-frame ATG is an internal ATG) in 82-83% on the
benchmark fragments and mixed inputs, 76% on the ab initio assemblies, and 63-70% with the assist. The vetoed UTRs
(w2r) and the protected truncated CDSs (r2w) look alike on the rule's own features: median K -1.4..-1.6 vs -1.8..-2.3
bits, median |u| 24-51 vs 9-39 nt, median coding per base of u 0.02-0.22 vs 0.00-0.13, 52-81% vs 60-88% with
|u| < 60 nt. V4 (floor only when coding(u) <= 0) selects even shorter segments (median |u| 6-24 nt) where the coding
score carries no information, and is not more precise (veto right in 62-83%). E34 counted only the UTR side of the
veto; this evaluation counts both sides.

### Complete inputs: the downstream-start risk
- On complete and complete_fullutr inputs no variant can lose an exact CDS: every 5'-open candidate there is a UTR
  extension, so keeping it open is never exact. r2w = 0 on all 32 complete-type benchmark runs (both arms); every
  complete input gains.
- The risk named in the task, a true start downstream of the first ATG with a poor first-ATG context, materialises as
  w2w_later: the ORF is now trimmed to the (wrong) first ATG instead of being left 5' partial. Exact F1 does not change,
  but a 5'-partial call becomes a complete call with a wrong start. Counts (ab initio, 8 complete inputs): V1 439
  (18% of its changes), V2_05 240, V2_1 142 (22%), V3_8.5 153 (15%), V3_4 278, V4 275 (28%); complete_fullutr 124-377;
  assisted 142-430. Median K of these is -2.0 bits (lower than the fixed UTRs' -1.5) and median |u| 12-21 nt.
- In the assemblies the same class is 4-116 per group; for the truth type "complete" alone the assemblies' net change
  is positive for every variant (per-type complete F1 with negatives +0.01 to +3.05 points; ab initio +0.12 to +3.05).

## Interpretation
- A negative Kozak score at the first in-frame ATG is evidence that the ATG is internal, i.e. that the ORF is truly
  truncated. On complete transcripts that evidence is always misleading (there are no truncated ORFs), so any floor
  helps them; on fragmented inputs it is right about 3 to 5 times as often as it is wrong. Flooring or down-weighting
  it moves the rule along the same complete-vs-fragment axis as alpha (E13), and loses on the equal-weight objective
  at every setting tried.
- The E34 "Kozak veto" class is real but is the minority side of a decision whose majority side E34 did not count.
  It cannot be separated with |u|, coding(u) or K alone; distinguishing the two would need information outside u, as
  in E33 (homology or transcript context).
- With homology hits, the assist removes much of the fragment cost (R1/R3 keep the truncated ORFs that have hits), and
  V2_1 / V3_8.5 pass the benchmark criteria. They still lower 4 of 6 Arabidopsis assemblies (to -0.81 / -0.66 at d32)
  and human d32 (-0.12 / -0.08), and lower 5'-partial and internal F1 on every assisted Arabidopsis assembly
  (-0.9 to -2.0). The gains come from complete inputs and the shallow assemblies with few truncated CDSs.

## Recommendation
Do not change the release rule. No variant is acceptable ab initio (all fail the 0.5-point per-input criterion on
21-32 of 48 benchmark inputs and lower the E13 objective). With homology hits, V3_8.5 (w- = 8.5) and V2_1 (K >= -1)
are the only candidates; their gain (+0.15 / +0.06 on the E13 objective; +0.48 / +0.27 on complete inputs) is
smaller than their cost on the deeper Arabidopsis assemblies (d8-d32: -0.23 to -0.81) and on 5' partials, and they would make
the ab initio and homology rules differ. Not clearly better; not recommended.

## Guard against overfitting
- The variants and the acceptance criteria were fixed before any result; no parameter was tuned. Seven variants
  were scored once each.
- What motivated the test: E34's diagnosis on human d1/d8/d32 and athal d1/d16 (the UTR side only). Those datasets
  did not drive a choice; they are where the variants do worst at depth (human d32, athal d16) and best at d1.
- The verdict rests on the benchmark fragments and mixed inputs (32 inputs, all down ab initio for every variant) and
  on the E13 objective, which E13 used to set w = 17. The assisted-arm candidates (V2_1, V3_8.5) were identified on
  the 24 benchmark homology inputs and then failed to help on the assemblies, which they had not been chosen on.

## Caveats
- Offline substitution: the variant changes only start decisions on the runs' own selections; Predict was not rerun.
  This is exact for the ab initio rule (refinement follows selection and changes nothing else), as the 0-mismatch V0
  checks show.
- Benchmark fragments are cut once inside every CDS (more truncation than real assemblies; E13 caveat); the mixed
  inputs and simulated assemblies, with realistic truncated shares, give the same sign.
- Single simulation per depth; human depth levels share genes. Differences of ~0.05 points are not interpretable.
- GM12878 (real assembly) not evaluated; its MANE-based truth cannot see wrongly trimmed truncated ORFs (E13 caveat).
