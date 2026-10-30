# E29. Start placement on 5'-open ORFs: error anatomy, real UTR/truncation data, UTR-length prior (2026-10-25)

## Question
On stranded fragments, PyTD release (3a793ff) finds the stop of 5'-partial references in 0.871 of
cases (median over 8 organisms) but leaves the ORF open at the fragment start (exact) in only 0.735;
Prodigal and CodAn reach 0.91-0.94. (1) Where do the wrong 5' ends lie? (2) Brian's hypothesis: when
an open end is trimmed to an ATG very close to the sequence start, the implied 5' UTR is very short,
so a prior on UTR length could inform keep-open vs trim. Is any effect an artefact of the synthetic
inputs (N(0.5,0.1) cuts, no UTR on 5' partials, UTRs clipped to 250 nt)? (3) Does a UTR-length prior
help on the benchmark and on GM12878?

## Data
`$P = /broad/hptmp/bhaas/TRANSDECODER_PAPER`.
- Release runs `$P/work/bench_v3/runs/<org>/<variant>/pytd_v7/` (stranded complete, complete_fullutr,
  fragments); GM12878 release runs `$P/work/gm12878_td/{lraa,stringtie,trinity}/final_v7_abinitio/`.
- Balanced mixed inputs `fragments_mixed_t{1,2,3}` (E28; 25% each complete / 5' partial / 3' partial /
  internal), run here with lambda = 0 (identical to the release; `preds.tsv` byte-identical to the
  coordinator's `pytd_v7` runs where they exist: dmel, creinhardtii, athal t1).
- Reference annotations `$P/data/annot/<org>/` (as listed in `benchmark_v2/truth_sources.tsv`); MANE
  v1.5 GTF for GM12878; assemblies and structures as in `genome_annotation/scripts/gm12878_final_v7.sh`.

## Method (scripts in `pipeline/derive/`, tables in `data/start_refinement/utr_prior/`)
- `03_utr_features.py`: for every selected ORF before refinement (`best_candidates.gff3`) re-derives
  the release decision (coding score, Kozak bits, alpha*seg + w*K, trim) and
  d = 5' UTR implied by trimming to the first in-frame ATG. Re-derived decisions agree with the runs'
  `revised_starts.gff3` in all 51 runs (0 mismatches).
- `03_utr_anatomy.py` -> `release.anatomy.tsv`, `release.decision_by_d.tsv` (benchmark),
  `release_gm12878.*` (GM12878 with transcript-level truth, below).
- `03_utr_annotation.py` + `03_utr_prior_table.py` -> `annotation_utr5_summary.tsv`,
  `trunc_d_summary.tsv`, `utr_prior_lr.tsv`: all eligible protein-coding transcripts (truth-builder
  checks, not one per gene, no added or clipped UTRs; human = all MANE; mouse = all GENCODE M39
  protein_coding without *_NF tags, 166,288 transcripts). "Open" UTR = no in-frame stop upstream of
  the start (the only UTRs for which a 5'-open ORF arises on a complete transcript).
- `03_utr_gm12878_truth.py` + `03_utr_gm12878_summary.py` -> `gm12878_five_prime_truncation.tsv`,
  `missing_frac_hist.tsv`: each assembled transcript (LRAA/StringTie GTF; Trinity primary minimap2
  alignment) that carries a MANE Select stop codon with the reference splicing from the stop to the
  transcript's 5' end is classed complete (contains the start; UTR length, open or not) or
  5' partial (begins inside the CDS; fraction of CDS missing). Also writes a transcript-level truth
  table, used for the GM12878 anatomy and offline simulation.
- Prototype (PyTD branch `exp/utr-prior`, worktree `$P/work/trees/PY_utr`, commit c9d7875, local
  only, not merged/pushed): `--start-utr-prior-weight lambda` (default 0 = off) adds
  lambda * w * log2[P(d|UTR)/P(d|truncated)] to the UTR score (w = 17, i.e. the prior's bits weigh
  like Kozak bits). The table is fixed and organism-independent, from the annotations only
  (`utr_prior_lr.tsv`: P(d|UTR) = open 5' UTRs whose start is the first in-frame ATG; P(d|trunc) =
  distance from a uniform cut inside a CDS to the first in-frame ATG; per-organism histograms averaged
  with equal weight). lambda = 1 was fixed a priori (Bayesian weighting), not tuned on the benchmark.
- `03_utr_runs.sh` (fresh runs, new dirs `$P/work/utr_prior/runs/`, `$P/work/utr_prior/gm12878_runs/`,
  eval `$P/work/utr_prior/gm12878_eval/`): lambda = 1 on stranded complete, complete_fullutr,
  fragments (8 orgs); lambda = 0 and 1 on the 24 mixed inputs and on GM12878 LRAA and StringTie.
  `03_utr_eval.py` -> `utr_prior_eval.{bench,bench_types,gm12878}.tsv`.
- `03_utr_prior_sim.py` -> `offline_sim.release*.tsv`: the same rule applied offline to the release
  decisions (stops are unchanged by refinement, so exact matches follow from the decisions), for
  lambda 0.5/1/2 and a "short" variant (penalty only for d < 20, Brian's hypothesis as stated). The
  offline lambda = 1 counts equal the fresh runs to within 1 exact match per input (ties at the 2nd
  decimal of the rounded scores).

## Results
### 1. Error anatomy (release, stranded)
- 5'-partial fragments, medians over organisms: stop found 0.871, exact 0.735; exact given stop
  0.835 (0.771 athal - 0.908 creinhardtii). **Every** stop-correct non-exact 5'-partial is a trim to
  an internal in-frame ATG (9,397 ORFs over 8 organisms; no other error type). Median d (fragment
  start to chosen ATG) 46-63 nt per organism (median 55.5); 18-25% have d < 20, 42-52% d < 50; the
  trim removes a median 8.6% of the reference CDS (5.9-11.4%).
- Start-present references (complete, complete_fullutr, 3' partials): left 5'-open through the true
  ATG in 0.7-8.5% of stop-correct predictions (median 5.5% complete, 4.6% complete_fullutr, 5.9% 3'
  partials); true UTR length of those median 24-128 nt per organism (pooled median 52-56); trims to an
  upstream in-frame ATG 0.5%; ATG-started ORFs with another start 1.8-2.0% (mostly creinhardtii, where
  the annotated start is often a later ATG; E09).
- Decision vs d (`release.decision_by_d.tsv`; complete + fragments, all organisms): truncated ORFs
  are wrongly trimmed at 30% (d 0-19), 26% (20-49), 20% (50-99), 14% (100-199), 10% (200-499), 7%
  (>= 500). True UTR extensions whose first ATG is the start are wrongly kept open at 36% (d 0-9),
  23% (10-19), 18% (20-49), 15% (50-99), 12% (100-199). Both error types peak at short d: a short
  segment carries little coding signal and its ATG has a partly missing Kozak window, so short d is
  where the rule has least information, not where one answer dominates. Same shape in every organism.
- GM12878 (transcript-level truth from MANE, release runs): 5' partials with the right stop were left
  open in 73.6% (LRAA), 81.2% (StringTie), 79.1% (Trinity); wrong trims median d 60-73 nt, 17-23% at
  d < 20; truncated ORFs at d < 20 were trimmed 37-45%. Complete references were left 5'-open in
  2.6-3.9% of stop-correct predictions.

### 2. Is it the test data?
- Annotated 5' UTRs (all eligible transcripts): median 97-497 nt; < 10 nt 0.2-1.6%, < 20 nt
  0.4-9.0%, < 50 nt 3.6-29.7%. But the UTRs that create an open ORF on a complete transcript (no
  in-frame stop; 4-47% of UTRs) are much shorter: median 27-118 nt; < 10 nt 0.9-8.1%; < 20 nt
  2.0-32.0%; < 50 nt 14.7-81.1% (spombe and pfalciparum shortest). Short open UTRs are therefore
  not rare in real annotations.
- Distance from a cut inside a CDS to the first in-frame ATG (annotation CDSs): < 20 nt in 9-12% of
  cuts, median 87-127 nt; uniform vs N(0.5,0.1) cut positions give the same d distribution (within
  0.01), so the cut model does not drive the d distribution of truncated ORFs.
- Pooled likelihood ratio log2 P(d|UTR)/P(d|trunc) (`utr_prior_lr.tsv`): -0.67 bits (d 3-9),
  +0.3..+1.1 (10-99), ~0 (100-149), -0.9..-6 (>= 150). The very-short-UTR effect exists but is weak
  (2-4 fold per organism at d < 10, opposite sign at 10-50 nt); the strong information is at long d,
  where the coding score already decides.
- Benchmark UTR side is realistic: open-UTR rate among complete inputs 31.9%, open-UTR median 75 nt,
  < 20 nt 9.5%; GM12878 assembled complete transcripts: open 26.8-38.3%, median 64-101 nt, < 20 nt
  5.1-14.1%.
- Real 5' truncation (GM12878, transcripts carrying a MANE stop with reference structure): begin
  inside the CDS 9.1% (LRAA, 1,365/14,949), 7.5% (StringTie, 1,017/13,521), 13.6% (Trinity,
  1,915/14,087). Fraction of CDS missing: median 0.50 / 0.55 / 0.69, IQR 0.24-0.73 / 0.17-0.83 /
  0.38-0.88, broad (StringTie bimodal; 11% / 19% / 7% miss < 10%). Benchmark fragments: all 5' partials
  miss 0.50-0.95 (IQR 0.53-0.62); mixed inputs 0.05-0.80 (median 0.47).
- Share of truly truncated ORFs among 5'-open ends (pi): GM12878 0.21 (LRAA), 0.20 (StringTie), 0.37
  (Trinity); benchmark fragments 0.76, mixed inputs 0.76, complete inputs 0.
  -> Per sequence, the 5'-partial error is real (the d distribution and the d-dependence are not
  artefacts of the cut model, and occur on GM12878 too). What the benchmark exaggerates is its
  weight: fragments and mixed inputs have ~4x (vs Trinity ~2x) the share of truncated open ends seen in
  real assemblies, so the 0.735 vs 0.91-0.94 gap measures behaviour that matters less in real output.
  (The comparators' high 5'-partial exact rates come with lower start accuracy on 3' partials, e.g.
  athal 3' partials exact 0.639 Prodigal, 0.448 CodAn vs 0.742 PyTD; `partial_types.tsv`.)

### 3. UTR-length prior (lambda = 1, fresh runs; `utr_prior_eval.*`)
- Stops identical in every run (stop F1 unchanged); all changes are start codons.
- complete: exact F1 +0.02 (creinhardtii) to +2.60 (spombe) points, all 8 up (median +0.7);
  complete_fullutr: +0.19 to +2.47, all up. Complete-transcript exact F1 did not fall.
- fragments: -0.60 to -1.30 points, all 8 down (median -0.74).
- mixed inputs (24 = 8 orgs x 3 trials), whole input exact F1: -0.28 to -1.43 points, all down. Per
  type (sensitivity, mean of 3 trials; range over the 24): complete +0.80 (-0.14..+2.51), 3' partial
  +0.72 (-0.09..+2.05), 5' partial -2.26 (-4.09..-0.96), internal -2.39 (-3.78..-0.62) points.
- GM12878: MANE exact given assembled LRAA 0.9291 -> 0.9336 (+28 CDSs), StringTie 0.9149 -> 0.9179
  (+27); GENCODE complete exact +47 / +51; distinct 5'-partial predictions 2,502 -> 2,370 /
  2,646 -> 2,624.
- Offline: lambda = 1 behaves like a slightly cheaper alpha increase (complete +1,098 exact vs
  fragments -1,040, against +1,069 / -1,388 for alpha +0.1) but combined with a lower alpha it is worse
  than alpha alone (lambda 1 + alpha -0.1: -242 / -169), so it adds little beyond the class-mix
  trade-off already set by alpha. The "short" variant (penalty only at d < 20): fragments +156 net
  exact (189 gained, 33 lost), mixed +573 net over 24 inputs (about +0.1 point each), complete -72 and
  complete_fullutr -75 (about -0.06 points), GM12878 -4 / -8 / +4 (LRAA / StringTie / Trinity).

## Conclusion
The 5'-partial deficit is a real property of the keep-open vs trim rule, not an artefact of the cut
model, but it is not mainly a short-UTR problem: wrong trims sit at a median 55 nt from the fragment
start, and true UTR extensions are wrongly kept open at the same short distances. UTR length carries
little evidence in the 10-150 nt range where most decisions fall (|log2 LR| <= 1.1), because UTRs that
leave an ORF open are themselves short. An annotation-derived UTR-length prior (lambda = 1) mostly
shifts the balance toward trimming: it helps complete transcripts and the real GM12878 assemblies
(+0.3-0.45 points MANE exact), and costs 0.6-1.4 points on fragments and mixed inputs; a penalty only
for very short implied UTRs changes ~0.1 point. The benchmark over-weights truncation (pi = 0.76 vs
0.20-0.37 in real assemblies), which is why the release's alpha = 0.4 (chosen with complete and
fragments weighted equally) looks conservative on fragments.

**Recommendation:** do not adopt the prior as a default; keep it unmerged. For the paper, report the
5'-partial gap with the GM12878 truncation rates and pi as context, and keep the existing guidance that
`--start_contrast_alpha` can be lowered for fragmentary assemblies and raised for full-length ones.

## Caveats
- GM12878 truth covers transcripts concordant with MANE from stop to 5' end only (MANE, one CDS per
  gene); discordant isoforms are excluded.
- Benchmark complete inputs clip UTRs to 250 nt, which creates extra long open UTRs (utr_first at
  d >= 200: 4,917 in complete vs 2,718 in complete_fullutr).
- Offline and fresh counts differ by <= 1 exact match per input (score ties).


## Use in the manuscript (2026-10-25)
Results, Start-codon refinement, paragraphs "Truncated ORFs and short UTRs are hard to tell apart"
(error anatomy vs distance, cut-model and annotation checks, UTR-length prior not adopted, GM12878
truncated share 0.20-0.37 vs 0.76) and Discussion Limitations (76% vs 20-37%).

## Notebook section (2026-10-30)
`notebooks/03_start_refinement.Rmd`, section "Truncated ORFs and short UTRs are hard to tell apart"
(anchors `#truncated-orfs-and-short-utrs-are-hard-to-tell-apart` and
`#how-many-open-ends-are-truly-truncated`), recomputes every number of the Results paragraphs and the
Discussion sentence from `data/start_refinement/utr_prior/` (31 checks, all agree). Notes:
- "about 5%" = median over organisms of the per-organism share (5.2%); per input type 5.5 / 4.6 / 5.9%.
- Distance error rates use the complete and fragments inputs only (complete_fullutr and mixed not
  pooled); with complete_fullutr pooled in, the d < 10 open rate would be 35%, not 36%.
- The 76% (fragments 0.7558, mixed trial 1 0.7594) is the transcript-level π of
  `gm12878_five_prime_truncation.tsv`, the same definition as for GM12878. A prediction-level
  cross-check (correct-stop open ends before refinement) gives 75.5% / 75.4%.
- LRAA π is 0.2086 (21%), StringTie 0.2035 (20%): the Results text "20% (LRAA, StringTie)" was changed
  to "20--21% (StringTie, LRAA)".
