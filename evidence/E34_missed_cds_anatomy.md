# E34. Anatomy of the missed coding regions, and of the 5' UTR-vs-coding decision on 5'-open candidates (2026-10-29)

Status: diagnostic, not in the paper. No TransDecoder code was changed.

## Question
E33 step 4 found that 80-90% of the non-exact 5'-partial and internal CDSs (after the homology start assist)
have no prediction with the true stop. Brian asked two things:
1. **Missed CDSs.** For every truth CDS with no prediction carrying its stop, where in the PyTD pipeline is it lost?
   What are the lost CDSs like, compared with the found ones, and what would the obvious relaxations
   (q 0.2, q 0.5, a lower minimum length, accepting the top candidate on empty transcripts) recover and cost?
2. **3' partials and the 5' UTR-vs-coding decision.** For truth 3' partials (and complete CDSs) whose candidate ORF is
   5'-open (no in-frame stop in the 5' UTR), how often is the open end trimmed to the true ATG, kept open, trimmed to
   another ATG, or the CDS missed? When it is kept open, is the cause a coding-looking UTR (model), a UTR too short for
   the length term, or a true start that is not the first in-frame ATG?

## Data (all read only; nothing resumed or modified)
`$P = /broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T = $P/TransDecoderPaper`. Manifests: `data/missed_cds/manifest.tsv`,
`manifest_gm12878.tsv`, `manifest_minlen30.tsv`, `manifest_gm12878_minlen30.tsv`.

| dataset | run (finished release run) | truth / labels |
|---|---|---|
| human.d1 (E30) | `$P/work/simtrinity/runs/pytd_v7` (PyTD 3a793ff) | `$P/data/simtrinity/labels/sim_trinity.*` |
| human.d8, human.d32 (E31) | `$P/work/simtrinity/depth/<d>/runs/pytd_v7` (3a793ff) | `$P/data/simtrinity/depth/<d>/labels/sim_trinity.*` |
| athal.d1.A, athal.d16.A (E32/E33) | `$P/work/start_classifier/athal_depth/A/depth/<d>/runs/pytd_v7` (a5b24c3, ab initio; E33 checked it byte-identical to E32's pytd_v7) | `$P/data/simtrinity_athal/depth/<d>/labels/sim_trinity.*` |
| athal.d1.B, athal.d16.B (E33) | `.../athal_depth/B/depth/<d>/runs/pytd_v7` (a5b24c3, full homology: DIAMOND `sprot_minus_3702` -k 1 1e-5 + Pfam GA on all LongOrfs; start assist on) | as A; hits `.../athal_depth/searches/<d>/{blastp.outfmt6,pfam.domtblout}` |
| gm12878.trinity.mane (E25/E29/E33) | `$P/work/gm12878_td/trinity/final_v7_abinitio` (3a793ff), linked as `$P/data/missed_cds/gm12878_inputs/run/` (symlinks; `preds.tsv` converted from its `transcripts.fa.transdecoder.gff3`) | E29 transcript-level MANE truth `$P/work/utr_prior/gm12878/trinity.truth.tsv` (12,172 complete + 1,915 5' partial; contigs carrying a MANE Select stop with the reference structure to their 5' end), copied with `n_codons` added to `gm12878_inputs/labels/mane.truth.tsv`; no negatives or exclusions |

- Scored set on the simulated assemblies: as `score_simtrinity.py` (excluded contigs removed; negatives = lncRNA-only contigs).
- GM12878: the truth covers only MANE-stop contigs and there are no negatives, so F1 values there are not meaningful;
  only the miss classes and the counts recovered are.

## Software
- PyTD a5b24c3 (`$P/work/trees/PY_hstart`, unchanged), imported for frame contrast, Markov scores, Kozak scoring,
  homology-start rules and `_revise_start_sites`. Ab initio output of a5b24c3 equals 3a793ff (E33).
- `$P/env/startclf` (Python 3.11.18, pandas 3.1.0, numpy 2.4.6); `$P/env/td` for the fresh PyTD runs.
- Paper repo at 38fbaaf + the E34 files.

## Method
Scripts in `pipeline/missed_cds/`.
1. **`01_trace_missed_v1.py <manifest> <dataset> <big_dir>`** reads only the run's own files: `input.fa`,
   `longest_orfs.{gff3,cds}`, `longest_orfs.cds.scores`, `decoy_qvalues.tsv`, `hexamer.scores`,
   `__checkpoints_TDpredict/top_training_orfs.cds`, `start_refinement.kozak.tsv`, `best_candidates.gff3`, `preds.tsv`,
   `predict.log` (long-ORF length), and for B the hit files the run was given.
   - It re-implements selection exactly as `predict.py` (`_select_best_orfs`, single best ORF per transcript:
     accepted = homology hit, or length >= the GC long-ORF length, or in-frame score > 0, frame contrast > 0 and
     q <= cutoff; ranked by homology count, in-frame score, length) and runs PyTD's own `_revise_start_sites` on the
     selection (training-ORF median contrast, the run's Kozak model, homology evidence for B).
   - **Check:** at q 0.05 the offline selection equals each run's `best_candidates.gff3` and the refined predictions
     equal its `preds.tsv` for every dataset (0 differences; `$P/data/missed_cds/<ds>/checks.txt`), and the scores
     reproduce the official exact F1 (human d1/d8/d32 0.8440/0.8315/0.8063; athal d1 A/B 0.8815/0.9138; d16 A/B
     0.8218/0.8989).
   - **Missed** = no prediction with the stop within 2 nt on the truth strand (`score_predictions.py` rule).
   - **Candidate** = the LongOrfs ORF whose 3' end is within 2 nt of the true 3' end (LongOrfs keeps one ORF per
     stop, from the first ATG after the previous in-frame stop or from the open frame start).
   - **Class** of each missed CDS, in this order:
     - (d) the selected prediction overlaps the CDS in its frame but ends at another stop;
     - (a) no candidate; reason `below_min_length` (the ORF LongOrfs would build for that stop is < 50 aa),
       `internal_or_3open_below_min`, or `other`;
     - (b) the candidate fails acceptance (no homology hit, shorter than the long-ORF length, and in-frame score <= 0,
       or frame contrast <= 0, or q > 0.05);
     - (c) the candidate is accepted but another ORF on the transcript ranks above it; what won: `homology_ranked`
       (more hits), `antisense`, `other_frame_overlap`, `nonoverlapping_orf`, plus whether it is longer, has a higher
       contrast, is in the same frame;
     - (e) anything else (none occurred).
   - Per CDS it records the candidate's length, in-frame score, frame contrast (total and per base), q-value,
     decoy stratum, eligibility, homology, rank, GC3 of the true CDS, source TPM, and for every missed CDS plus a 15%
     random sample of found ones (seed `E34:<dataset>`) the true CDS's own per-base frame contrast and in-frame
     Markov score (`true_contrast_pb`, `true_s0_pb`).
   - **Part 2:** for truth complete and 3' partial CDSs whose candidate is 5'-open (LongOrfs type 5prime_partial or
     internal), it recomputes the release start rule exactly as `_revise_start_sites`: first in-frame ATG (offset |u|),
     coding(u) = frame contrast of u / training median per-base contrast, Kozak K (bits) of that ATG,
     UTR(u) = 0.4|u| + 17K; whether the true start is the first in-frame ATG (within 2 nt); and the outcome from the
     final predictions (exact = 5' end at the true ATG; kept open; trimmed to another ATG; missed). For B, the run's
     `start_refinement.homology.tsv` gives the decision source.
   - **Counterfactuals** (offline, same code path): q 0.1, 0.2 (`--sensitive`), 0.5, 1.0; on transcripts with no
     accepted ORF, the top-ranked eligible ORF (in-frame score > 0, contrast > 0) or the top-ranked ORF of any score.
2. **`02_fresh_runs_v1.sh` / `02_fresh_runs_v2.sh`** (v2 takes the manifest as an argument), via **`03_fresh_lanes_v1.sh`**
   (2 processes at a time, `nice -n 10`): fresh runs in new directories `$P/work/missed_cds/<mode>/<ds>/runs/pytd_v7`:
   - `minlen30`: fresh LongOrfs `-S -G Universal -m 30` and a default Predict (human d1/d8/d32, athal d1.A/d16.A,
     GM12878 Trinity);
   - `sensitive`: the release run's LongOrfs files copied (no Predict checkpoints) and Predict `--sensitive`
     (human d1, athal d1.A). **`07_check_sensitive_v1.py`**: the program's predictions equal the offline q 0.2
     predictions exactly (human d1 46,246 / 46,246; athal d1.A 21,460 / 21,460; 0 differences;
     `data/missed_cds/check_sensitive_fresh_vs_offline.tsv`).
   - The `minlen30` runs are traced with `01` (`$P/data/missed_cds/minlen30/<ds>/`); their misses are compared CDS by CDS
     with the release run.
3. **`04_summarize_v1.py`** -> the tables below.
4. **`05_source_annotation_v1.py`**: for each truth CDS, the first source transcript's protein vs the longest
   annotated protein of its gene among the simulated coding transcripts (`kept.tsv` / `kept.fa` of E30 / E32):
   `main` (the CDS protein is that protein or a substring of it), `shared_frame` (shares a 10-aa word: same frame,
   other exons or ends), `novel` (no shared word); CDS length / longest protein; GENCODE tag (human).
5. **`06_utr_coding_check_v1.py`**: for part-2 candidates with the true start at the first ATG and the stop predicted,
   whether u, translated in the ORF frame, shares a 10-aa word with any annotated protein of the same gene (u is
   coding in another isoform: an N-terminal extension). u < 30 nt is "too short to test".

## Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper; S=$T/pipeline/missed_cds; M=$T/data/missed_cds
PY="nice -n 10 env PYTHONPATH=$P/work/trees/PY_hstart $P/env/startclf/bin/python"
for d in human.d1 human.d8 human.d32 athal.d1.A athal.d1.B athal.d16.A athal.d16.B; do
  $PY $S/01_trace_missed_v1.py $M/manifest.tsv $d $P/data/missed_cds > $P/data/missed_cds/$d.log 2>&1; done
$PY $S/01_trace_missed_v1.py $M/manifest_gm12878.tsv gm12878.trinity.mane $P/data/missed_cds
nohup bash $S/03_fresh_lanes_v1.sh &       # minlen30 x 5 sims, sensitive x 2; log $P/work/missed_cds/logs/lanes_v1.log
nohup bash $S/02_fresh_runs_v2.sh minlen30 gm12878.trinity.mane $M/manifest_gm12878.tsv &
for d in human.d1 human.d8 human.d32 athal.d1.A athal.d16.A; do $PY $S/01_trace_missed_v1.py $M/manifest_minlen30.tsv $d $P/data/missed_cds/minlen30; done
$PY $S/01_trace_missed_v1.py $M/manifest_gm12878_minlen30.tsv gm12878.trinity.mane $P/data/missed_cds/minlen30
$PY $S/07_check_sensitive_v1.py $P/data/missed_cds human.d1 athal.d1.A > $M/check_sensitive_fresh_vs_offline.tsv
$PY $S/04_summarize_v1.py $M/manifest.tsv,$M/manifest_gm12878.tsv $P/data/missed_cds $M
$PY $S/05_source_annotation_v1.py $M/manifest.tsv $P/data/missed_cds $M/b_source_annotation.tsv
$PY $S/06_utr_coding_check_v1.py $M/manifest.tsv $P/data/missed_cds $M/start5open_utr_coding_check.tsv
```
Outputs: tables in `data/missed_cds/`; per-CDS tables (`per_cds.tsv.gz`, `start5open.tsv.gz`, `source_annot.tsv.gz`,
`utr_coding_check.tsv.gz`), the offline selections and refinements (`td/`) and logs in `$P/data/missed_cds/<ds>/`.

## Results, part 1: where the missed CDSs are lost (`classification.tsv`)
Missed = no prediction with the true stop. All truth types together; per-type rows are in the table file.

| dataset | truth | exact | stop, not exact | missed | (a) | (b) | (c) | (d) | (e) | missed / non-exact |
|---|---|---|---|---|---|---|---|---|---|---|
| human d1 | 33,939 | 28,106 | 2,403 | 3,430 | 0 | 3,217 | 213 | 0 | 0 | 0.59 |
| human d8 | 20,160 | 16,224 | 1,809 | 2,127 | 0 | 2,081 | 46 | 0 | 0 | 0.54 |
| human d32 | 13,846 | 10,611 | 1,564 | 1,671 | 0 | 1,665 | 6 | 0 | 0 | 0.52 |
| athal d1 A | 21,053 | 18,165 | 1,862 | 1,026 | 0 | 994 | 32 | 0 | 0 | 0.36 |
| athal d1 B | 21,053 | 19,090 | 1,467 | 496 | 0 | 442 | 54 | 0 | 0 | 0.25 |
| athal d16 A | 19,206 | 15,331 | 2,753 | 1,122 | 0 | 1,111 | 11 | 0 | 0 | 0.29 |
| athal d16 B | 19,206 | 17,013 | 1,621 | 572 | 0 | 566 | 6 | 0 | 0 | 0.26 |
| GM12878 Trinity, MANE | 14,087 | 11,539 | 1,088 | 1,460 | 468 | 786 | 180 | 26 | 0 | 0.57 |

- Per type (human d1, missed / class b / class c): complete 2,188 / 1,990 / 198; 5' partial 412 / 404 / 8;
  3' partial 286 / 279 / 7; internal 544 / 544 / 0. The complete row reproduces E30 (2,188 complete CDSs missed at stop
  level).
- The missed share of the non-exact CDSs here (ab initio release runs) is 0.46-0.73 per type in human. E33 step 4's
  83-90% was after the homology start assist, which removes most of the start errors and leaves the selection errors.
- **(b) dominates on every simulated assembly (89-99.6% of misses; 89% only in athal d1 B, where homology has already rescued many).** Sub-reasons, human d1: q > 0.05 2,686, frame
  contrast <= 0 365, in-frame score <= 0 166. In 2,963 of the 3,217 no ORF at all is predicted on the contig; in 254
  (237 of them complete CDSs) another, failing-or-not, ORF was selected instead.
- **(a) occurs only on the real assembly**, because the simulated labeller excludes CDSs under 50 codons
  (`short_*`). GM12878: 443 MANE CDS fragments (389 of the 5' partials) are below LongOrfs' 50-aa minimum (median
  24 codons for the 5' partials); 25 "other" are contigs whose assembled CDS has an in-frame premature stop
  (19 of 25).
- **(d) occurs only on the real assembly**: 26 contigs, all with an in-frame premature stop in the assembled CDS
  (assembly or sequence error), so the predicted ORF in the right frame ends early.
- **(c)** is small: human d1 213 (complete 198), d8 46, d32 6, athal 6-54, GM12878 180 (complete 168). Human d1
  complete: the winner is a non-overlapping ORF elsewhere on the contig in 143, an overlapping ORF in another frame in
  55; it is longer in 144 and has a higher frame contrast in 135 (median winner 900 nt vs candidate 190 codons). In
  athal B, 32 of 54 lost to an ORF with more homology hits (`c_winners.tsv`). GM12878: 154 of 168 lost to a
  non-overlapping, mostly longer ORF on the same contig (single-best-ORF selection on multi-ORF contigs).
  No antisense winners: all runs are strand-specific (`-S`).

## Results, part 1: what the class (b) misses are (`b_characterisation.tsv`, `b_vs_found_by_length.tsv`, `b_by_stratum.tsv`, `b_source_annotation.tsv`)
All types together, missed (b) vs found (stop predicted):

| dataset | n (b / found) | codons, median (b / found) | < 100 codons (b / found) | contrast per base (b / found) | true CDS in-frame score <= 0 (b) | GC3 outside training 5-95% (b / found) | source TPM < 1 (b / found) |
|---|---|---|---|---|---|---|---|
| human d1 | 3,217 / 30,509 | 83 / 368 | 0.69 / 0.05 | 0.027 / 0.066 | 4.1% | 0.10 / 0.13 | 0.57 / 0.48 |
| human d8 | 2,081 / 18,033 | 80 / 272 | 0.75 / 0.11 | 0.026 / 0.070 | 5.4% | 0.16 / 0.15 | 0.11 / 0.15 |
| human d32 | 1,665 / 12,175 | 79 / 158 | 0.82 / 0.23 | 0.022 / 0.078 | 11.2% | 0.31 / 0.18 | 0.03 / 0.03 |
| athal d1 A | 994 / 20,027 | 77 / 361 | 0.78 / 0.06 | 0.019 / 0.054 | 11.4% | 0.64 / 0.30 | 0.29 / 0.10 |
| athal d1 B | 442 / 20,557 | 74 / 354 | 0.85 / 0.08 | 0.019 / 0.053 | 8.8% | 0.63 / 0.31 | 0.27 / 0.11 |
| athal d16 A | 1,111 / 18,084 | 69 / 170 | 0.89 / 0.27 | 0.022 / 0.063 | 12.5% | 0.54 / 0.33 | 0.00 / 0.00 |
| athal d16 B | 566 / 18,634 | 68 / 164 | 0.94 / 0.28 | 0.023 / 0.062 | 11.0% | 0.51 / 0.33 | 0.01 / 0.00 |
| GM12878 MANE | 786 / 12,627 | 82 / 380 | 0.69 / 0.03 | 0.033 / 0.068 | 3.2% | 0.11 / 0.13 | n/a |

Training-set GC3 (5-95%): human 0.31-0.80 (d1), athal 0.34-0.46 (d1) / 0.35-0.47 (d16); class (b) GC3 median athal
0.47 vs found 0.41 and training 0.39.
- **Length.** The misses are short: median 68-83 codons (per type up to 90), 69-94% under 100 codons (found: 5-28%), none at >= 300 codons. The fraction
  of CDSs found rises with length (human d1: 50-99 codons 0.43, 100-149 0.77, 150-199 0.93, >= 200 0.98; athal d1 A:
  0.61, 0.90, 0.97, >= 0.996).
- **Coding signal.** Within the same length bin the misses have about 40% of the per-base frame contrast of the found
  CDSs (human d1, 50-99 codons: 0.030 vs 0.075; 100-149: 0.025 vs 0.066; athal d1 A 50-99: 0.020 vs 0.050), and a lower
  per-base in-frame Markov score (0.054 vs 0.091; athal 0.038 vs 0.060). So they are short **and** weak; at equal
  length the found ones score 2-2.5 times higher. 4-12% of the misses have a true-CDS in-frame score <= 0 (vs <= 0.5%
  of the found ones): the trained model does not recognise them as coding at all.
- **q-values.** Of the (b) candidates, 13-27% are ineligible (in-frame score or contrast <= 0, no q-value), 31-45%
  have q <= 0.2 and 55-65% q <= 0.5 (all types; per type 13-42%, 27-64%, 35-78%; median q 0.07-0.62). Most sit in the short-ORF strata
  (e.g. human d1, ORF-length bins 3-7, transcripts < 500 nt (bin 0), 5'/3'-open), where the found CDSs of the same
  stratum have about twice the contrast (`b_by_stratum.tsv`).
- **GC3.** Human: no consistent difference at d1-d8 (0.10 vs 0.13 outside the training range); at d32 the missed
  internal fragments are GC3-poor (median 0.39 vs 0.55; 38% vs 23% outside). Arabidopsis: the misses are GC3-rich
  relative to the training set, 51-64% outside its 5-95% range vs 30-33% of the found ones (the A6 GC3-minority
  pattern of the benchmark).
- **Expression.** In human the misses are not less expressed than the found ones (d1: source TPM < 1 in 57% vs 48%;
  d8 11% vs 15%; d32 3% vs 3%). In Arabidopsis d1 they are enriched for TPM < 1 (29% vs 10%; d16 n/a, both ~0),
  i.e. lowly expressed short genes; within the human sets, once a CDS is assembled, being found depends on its
  sequence, not its depth.
- **What the CDSs are (source annotation).** Human complete misses are mostly short alternative-isoform CDSs: median
  45% of the gene's longest protein (found: 94%), 13% MANE Select (found: 32%), 78% sharing its frame but not equal to
  it; only 4% are in a frame not used by the gene's longest protein. Arabidopsis misses are the genes' main proteins
  (90% `main`; Araport/TAIR genes with genuinely short proteins).
- **Homology (athal B).** By construction none of the (b) candidates has a hit: any DIAMOND hit at 1e-5 or Pfam GA
  domain makes an ORF accepted. Homology retention halves the misses (d1: 1,026 -> 496; d16: 1,122 -> 572); the rest
  have no Swiss-Prot (minus Arabidopsis) hit and no Pfam domain. Their median length is 74 codons (d1).

## Results, part 1: counterfactuals (`counterfactuals.tsv`)
Recovered = missed in the release run, stop predicted in the variant (exact in brackets). Negatives = lncRNA-derived
contigs with a prediction. F1 changes in points vs the release run. Athal negatives are too few (1,372 at d1, 154 at
d16) to measure a false-positive cost; GM12878 has no negatives (its "extra predictions" are on all 429k contigs).

| dataset | variant | recovered (exact) of missed | lost | extra predictions on negatives | negatives with a prediction | d stop F1 | d exact F1 |
|---|---|---|---|---|---|---|---|
| human d1 | q 0.2 | 1,057 (888) of 3,430 | 41 | +1,334 | 8.1% -> 14.5% | -0.33 | -0.56 |
| | q 0.5 | 1,831 (1,559) | 79 | +3,705 | 25.9% | -2.40 | -2.58 |
| | top eligible on empty | 2,247 (1,944) | 0 | +9,364 | 53.1% | -8.39 | -8.05 |
| | top any on empty | 2,411 (2,058) | 0 | +12,999 | 70.6% | -11.87 | -11.33 |
| | LongOrfs -m 30 | 0 | 2 | +13 | 8.2% | -0.03 | -0.03 |
| human d8 | q 0.2 | 890 (704) of 2,127 | 11 | +375 | 12.0% -> 18.3% | +1.38 | +0.75 |
| | q 0.5 | 1,371 (1,092) | 23 | +924 | 27.7% | +1.22 | +0.41 |
| | top eligible on empty | 1,585 (1,288) | 0 | +2,151 | 48.5% | -1.09 | -1.62 |
| | LongOrfs -m 30 | 0 | 1 | +8 | 12.1% | -0.02 | -0.02 |
| human d32 | q 0.2 | 791 (575) of 1,671 | 9 | +121 | 16.6% -> 24.2% | +2.57 | +1.40 |
| | q 0.5 | 1,090 (803) | 15 | +234 | 31.2% | +3.20 | +1.72 |
| | top eligible on empty | 1,237 (926) | 0 | +493 | 47.5% | +2.84 | +1.39 |
| | LongOrfs -m 30 | 0 | 0 | +8 | 17.1% | -0.03 | -0.02 |
| athal d1 A | q 0.2 | 384 (316) of 1,026 | 1 | +57 | 6.1% -> 10.2% | +0.79 | +0.56 |
| | q 0.5 | 547 (460) | 8 | +108 | 13.9% | +1.00 | +0.74 |
| | top eligible on empty | 700 (608) | 0 | +494 | 42.1% | +0.35 | +0.19 |
| | LongOrfs -m 30 | 0 | 0 | +2 | 6.2% | -0.01 | -0.01 |
| athal d1 B | q 0.2 | 160 (130) of 496 | 0 | +56 | 6.8% -> 10.9% | +0.25 | +0.14 |
| | q 0.5 | 239 (198) | 0 | +107 | 14.6% | +0.30 | +0.17 |
| athal d16 A | q 0.2 | 505 (370) of 1,122 | 3 | +5 | 2.6% -> 5.8% | +1.33 | +0.82 |
| | q 0.5 | 633 (487) | 4 | +5 | 5.8% | +1.65 | +1.14 |
| | top eligible on empty | 831 (659) | 0 | +33 | 24.0% | +2.03 | +1.49 |
| | LongOrfs -m 30 | 0 | 0 | 0 | 2.6% | 0.00 | 0.00 |
| athal d16 B | q 0.2 | 255 (188) of 572 | 0 | +5 | 5.8% | +0.66 | +0.37 |
| | q 0.5 | 325 (251) | 0 | +5 | 5.8% | +0.84 | +0.52 |
| GM12878 MANE | q 0.2 | 232 (201) of 1,460 | 4 | (+8,815 predictions on all contigs) | n/a | n/a | n/a |
| | q 0.5 | 421 (362) | 16 | (+26,285) | n/a | n/a | n/a |
| | top eligible on empty | 552 (473) | 0 | (+105,421) | n/a | n/a | n/a |
| | LongOrfs -m 30 | 5 (5) | 0 | (+163 predictions on all contigs) | n/a | n/a | n/a |

- q 0.1 and q 1.0 and "top ORF of any score" for every dataset are in the table file.
- **q 0.2 (`--sensitive`)** recovers 31% (human d1) to 47% (human d32) of the misses (athal 32-45%). On human d1, where
  negatives are plentiful (20,815), it costs 1.3 lncRNA-contig predictions per recovered CDS and lowers exact F1 by
  0.56 points; at d8/d32, where the negatives are fewer (5,885 / 1,594), the same rule raises exact F1 by 0.75 / 1.40.
  The sign of the F1 change follows the negative:coding ratio of the input, not the rule.
- **q 0.5** recovers 48-65%; on human d1 it predicts on 26% of the negatives (-2.6 exact F1 points).
- **Accepting the top candidate on empty transcripts** (eligible / any score) recovers 62-93% but predicts on 42-71%
  of the negatives at d1 (24-44% at athal d16), approaching the longest-ORF baseline (70.6% on human d1).
- **A lower minimum length** recovers nothing on the simulated assemblies, by construction: the labeller excludes
  CDSs < 50 codons and every scored CDS has a LongOrfs candidate (class (a) = 0). Its only effect is +0-13 negative
  predictions. On GM12878 it is the one counterfactual that addresses class (a) (443 MANE fragments < 50 aa):
  with `-m 30`, class (a) falls from 468 to 284 MANE CDSs, but class (b) rises from 786 to 965: the newly reachable
  30-49-aa fragments almost all fail the coding test. Net: 5 CDSs recovered (all exact) for +163 predictions.
  Fragments this short carry too little coding signal for the q-value test.

## Results, part 2: 5'-open candidates of 3' partials and complete CDSs (`start5open_outcomes.tsv`)
Candidates = truth CDSs whose LongOrfs candidate is 5'-open (no in-frame stop between the contig's 5' end and the
true start). The right call is always "trim to the true ATG".

| dataset | type | 5'-open candidates | true start = first ATG | exact (trimmed to true ATG) | kept open | trimmed to another ATG | missed |
|---|---|---|---|---|---|---|---|
| human d1 | 3' partial | 474 | 450 | 363 (77%) | 48 | 4 | 59 |
| | complete | 7,363 | 6,960 | 5,870 (80%) | 775 | 251 | 467 |
| human d8 | 3' partial | 756 | 720 | 578 (76%) | 104 | 8 | 66 |
| | complete | 4,100 | 3,971 | 3,382 (82%) | 448 | 80 | 190 |
| human d32 | 3' partial | 687 | 653 | 518 (75%) | 102 | 2 | 65 |
| | complete | 1,718 | 1,677 | 1,435 (84%) | 204 | 21 | 58 |
| athal d1 A | 3' partial | 201 | 199 | 111 (55%) | 63 | 2 | 25 |
| | complete | 4,996 | 4,980 | 3,859 (77%) | 1,007 | 10 | 120 |
| athal d1 B | 3' partial | 201 | 199 | 127 (63%) | 59 | 2 | 13 |
| | complete | 4,996 | 4,980 | 3,937 (79%) | 989 | 11 | 59 |
| athal d16 A | 3' partial | 525 | 517 | 372 (71%) | 114 | 2 | 37 |
| | complete | 2,181 | 2,177 | 1,794 (82%) | 347 | 2 | 38 |
| athal d16 B | 3' partial | 525 | 517 | 383 (73%) | 114 | 2 | 26 |
| | complete | 2,181 | 2,177 | 1,841 (84%) | 326 | 3 | 11 |
| GM12878 MANE | complete | 3,242 | 3,034 | 2,660 (82%) | 293 | 129 | 160 |

("true start = first ATG" excludes a few candidates with no in-frame ATG: human 10-18 3' partials.)

**Why the kept-open ones were kept** (`start5open_kept_open.tsv`; ab initio rule unless stated):

| dataset, type | kept open | true start not the first ATG | homology keep (R1/R3) | u scores coding: coding(u) >= 0.4\|u\| | Kozak veto: coding(u) < 0.4\|u\| but K < 0 | \|u\| < 60 nt |
|---|---|---|---|---|---|---|
| human d1 complete | 775 | 127 | - | 494 | 154 | 275 |
| human d1 3' partial | 48 | 6 | - | 33 | 9 | 15 |
| human d8 complete | 448 | 42 | - | 274 | 132 | 194 |
| human d8 3' partial | 104 | 5 | - | 61 | 35 | 55 |
| human d32 complete | 204 | 17 | - | 114 | 73 | 119 |
| human d32 3' partial | 102 | 5 | - | 58 | 37 | 63 |
| athal d1 A complete | 1,007 | 5 | - | 646 | 356 | 661 |
| athal d1 A 3' partial | 63 | 0 | - | 37 | 26 | 45 |
| athal d1 B complete | 989 | 5 | 94 | 561 | 329 | 651 |
| athal d16 A complete | 347 | 1 | - | 199 | 147 | 247 |
| athal d16 A 3' partial | 114 | 2 | - | 73 | 38 | 90 |
| athal d16 B 3' partial | 114 | 2 | 11 | 67 | 33 | 90 |
| GM12878 MANE complete | 293 | 51 | - | 204 | 38 | 105 |

Scores of the kept-open ones (medians):
- "u scores coding": |u| 87-107 nt (human d1), 33-51 (athal); coding(u) 0.71-0.96 per base of u, i.e. about
  typical coding signal (e.g. human d1 complete: coding 72.0, UTR score 50.8, K +0.20).
- "Kozak veto": coding(u) about 0 per base in athal (-0.06 to +0.07), 0.12-0.29 in human, with K -1.0 to -1.9 bits;
  |u| 15-36 nt in athal, 33-69 in human. E.g. athal d1 complete: |u| 28.5, coding 0.6, K -1.39, UTR score
  0.4 x 28.5 + 17 x (-1.39) = -10.5. These UTRs look like UTR; the segment is too short for 0.4|u| to outweigh a
  negative Kozak score (one Kozak bit = 42 nt of length term).
- "True start not the first ATG": the first-ATG rule cannot be right. Human d1 complete: 127 kept open plus 251 trimmed
  to an earlier ATG, 378 of 7,363 candidates (5.1%); human d8/d32 3.0% / 2.2%; athal <= 0.3%; 3' partials 1-2%;
  GM12878 complete 180 of 3,242 (5.6%).

**Kept-open rate by |u|** (trim-needed candidates: true start = first ATG, stop predicted; `start5open_by_u.tsv`):

| dataset, complete | 3-29 nt | 30-59 | 60-149 | 150-299 | >= 300 | all |
|---|---|---|---|---|---|---|
| human d1 | 0.15 | 0.09 | 0.10 | 0.07 | 0.14 | 0.10 (648 / 6,518) |
| human d32 | 0.15 | 0.16 | 0.09 | 0.05 | 0.27 | 0.12 |
| athal d1 A | 0.25 | 0.19 | 0.18 | 0.22 | (3) | 0.21 (1,002 / 4,861) |
| athal d16 A | 0.19 | 0.16 | 0.14 | 0.12 | (2) | 0.16 |

Fraction of these real UTR segments with coding(u) >= 0.4|u|: human 0.29-0.32 at 3-29 nt, 0.05-0.09 at 150-299 nt;
athal 0.29 -> 0.13-0.24. Per-base coding of the trim-needed segments has median -0.10 (human d1) and 0.00-0.01 (athal);
the median Kozak score is 1.7 bits (human d1) vs 0.9 (athal d1).

**Are the coding-looking UTRs coding elsewhere?** (`start5open_utr_coding_check.tsv`) Of the kept-open UTRs that
score as coding, 15-18% (human) and 0-4% (athal) share a 10-aa word with another annotated protein of the same gene
(an N-terminal extension in another isoform), vs 2-9% of the correctly trimmed UTRs. The rest have no annotated
coding counterpart.

## Interpretation
- **Missed CDSs are a coding-test problem, not a selection or LongOrfs problem.** On the simulated assemblies,
  89-99.6% of missed CDSs (no prediction with the true stop) have a LongOrfs candidate that fails acceptance (class b);
  0.4-11% lose to another ORF on the contig (class c, mostly a longer ORF elsewhere on the same contig), and none are
  below the LongOrfs minimum or carry a same-frame wrong stop. On the real GM12878 assembly, 30% of the MANE misses are
  fragments under 50 aa (class a) and 2% have an in-frame premature stop (class d).
- **The class (b) misses are short and weakly coding.** Median 68-83 codons (69-94% under 100 vs 5-28% of found CDSs),
  and at the same length they carry ~40% of the per-base frame contrast of found CDSs; 4-12% score <= 0 in frame.
  In human they are mostly short alternative-isoform CDSs (median 45% of the gene's longest protein; 13% MANE); in
  Arabidopsis they are short main proteins whose GC3 lies outside the training set's range (51-64% vs 30-33%).
  Expression does not explain them in human (Arabidopsis d1 misses are enriched for TPM < 1). With homology (athal B) half are rescued; the rest have no Swiss-Prot (minus
  Arabidopsis) or Pfam evidence.
- **Relaxing the test trades them for noncoding calls at about 1:1 or worse.** q 0.2 recovers 31-47% of the misses
  and costs 1.3 lncRNA-contig predictions per recovered CDS on human d1 (exact F1 -0.6), but gains +0.4..+1.4 points
  where negatives are scarce (human d8/d32, Arabidopsis). Accepting the top candidate everywhere approaches the
  longest-ORF baseline's 50-70% negative rate. A lower minimum length recovers nothing on the simulations (the labeller
  removes < 50-codon CDSs) and only 5 of 1,460 MANE misses on GM12878, whose 30-49-aa fragments then fail the coding test.
- **The 5' UTR-vs-coding decision on 5'-open candidates fails mainly because UTR segments score as coding, not because
  of the first-ATG rule.** For true 3' partials and complete CDSs with an open 5' UTR, the end is kept open in 10-16%
  (human) and 15-36% (Arabidopsis) of the trim-needed cases. About 55-70% of those keeps are UTRs whose frame contrast
  is >= 0.4 of typical coding per base (median 0.7-0.8, near typical coding), and only 15-18% (human) / 0-4% (athal) of
  them are coding in another annotated isoform: a model-discrimination problem, where the UTR's coding-score
  distribution overlaps that of coding sequence (worse at short |u|: 29-32% of real UTRs under 30 nt pass 0.4 per base
  vs 5-9% at 150-299 nt). 13-42% (human 19-36%, Arabidopsis 33-42%, GM12878 13%) are short UTRs (median 15-36 nt in Arabidopsis) that look like UTR but whose length
  term cannot outweigh a negative Kozak score (a weak-Kozak, short-|u| problem, larger in Arabidopsis, whose Kozak
  scores are weaker). A true start that is not the first ATG accounts for 2-5% of the human candidates and < 0.5% in
  Arabidopsis.

## Checks
- Offline selection and refinement reproduce every release run exactly (selected ORF IDs and final coordinates;
  `checks.txt`), and the official exact F1 to 4 decimals.
- The fresh `--sensitive` Predict runs equal the offline q 0.2 predictions (0 differences, 2 datasets).
- Each `minlen30` run (5 simulated + GM12878) is reproduced by its own offline selection (0 differences).
- Human d1 complete misses = 2,188, as in E30.
- GM12878: missed 5' partials 727 of 1,915 (ab initio stops predicted 1,188, as in E33).

## Caveats
- Single simulation and run per level; differences of ~0.01 F1 are not interpretable.
- The counterfactual F1 values on the simulations depend on the number of negatives (lncRNA-derived contigs at real
  expression), which falls with depth and is small in Arabidopsis; the human d1 numbers are the informative cost.
- "Contrast per base" of the true CDS uses all missed CDSs and a 15% random sample of found ones.
- `true start = first ATG` uses a 2-nt tolerance, as the scorer.
- `05`/`06` compare to annotated proteins of the simulated transcript set (GENCODE v50 kept protein-coding
  transcripts; TAIR12 mRNAs); "novel" means no identical 10-aa word, not proof of a different frame.
- GM12878 truth is MANE-stop contigs only (E29); its F1 values are not meaningful.
