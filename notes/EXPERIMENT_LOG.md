# Experiment log: what was tried, what worked, what did not

Created 2026-10-28. Internal; not part of the manuscript.

## Purpose and how to use it
This is a quick map of the method changes, analyses and benchmark design choices tried for the
TransDecoder (Perl) / PyTransDecoder (PyTD) paper. For each one it gives the decision and the lesson.
**The evidence records (`evidence/E*.md`) remain the detailed provenance:** data, commands, software
versions and full tables. This log only summarises them and points to them. If a number here
disagrees with its E-record, the E-record wins; fix the log.

- Read the summary table first, then the entry for the topic you care about.
- **Outcomes:** ADOPTED (in the release code or the paper's design), NOT ADOPTED (tried and rejected),
  PARTIAL (adopted in part, or adopted during development and later replaced), OPEN (no decision yet,
  or the facts needed are not recorded), ABANDONED (stopped before a result).
- **Paths:** `$P` = `/broad/hptmp/bhaas/TRANSDECODER_PAPER`, `$T` = `$P/TransDecoderPaper`.
  "notes/NN" = `$T/notes/NN_*.md`.
- **Code versions:** the release is PyTD `main` 3a793ff and Perl `devel` 41fc5ee (E24, E25). "Development
  defaults" (`pytd_dev_v2`, PyTD 48c870c / Perl de3f57e) means the 2026-10-20 rule with the fixed
  short-ORF test, before the q-value test.
- **Adding entries:** when an experiment ends, add or update its entry here as well as its E-record.

## Summary

| ID | Topic | Outcome | One-line lesson |
|---|---|---|---|
| A1 | Perl/PyTD parity fixes | ADOPTED | Compare intermediates (hexamer table, 6-frame scores, selections), not just final ORFs |
| A2 | Regenerated GC long-ORF table | ADOPTED | The 2017 table had no code; low-GC rows were 3-10% low from a small sample |
| A3 | GC from G+C counts; code-aware rule; tie-break | ADOPTED | 2q(C) is wrong on stranded data; small estimate changes move thresholds a full step |
| A4 | Keep the long-ORF (length) fallback | ADOPTED | 19% of the CDSs it rescues have no homology hit; it also admits some false calls |
| A5 | Drop the S1 > 0 requirement | NOT ADOPTED | Gains almost no true CDSs; floods random sequence with candidates |
| A6 | GC3-matched / alternative training sets | OPEN | GC3 minority classes explain many failures; no decision on using it |
| A7 | Contrast test on all ORFs as `--precise` | NOT ADOPTED | Post-filters on the chosen ORF did not beat a stricter test; replaced by q settings |
| A8 | Smooth short-ORF thresholds (taper etc.) | NOT ADOPTED | The cliff matters; the shape of the replacement matters less than its strength |
| A9 | Target-decoy q-value acceptance | ADOPTED | Stratify by transcript length too; it costs some real translated small ORFs |
| A10 | Long-ORF accept under the q-value test | ADOPTED (`gc`) | A decoy-derived length admits long random ORFs in TGA-only codes |
| A11 | Named q settings (`--precise`, `--sensitive`) | ADOPTED | Remaining lncRNA calls are long, strong ORFs no cutoff removes |
| A12 | Presets for small inputs | ADOPTED | Self-training fails below ~1,000 sequences; presets must include fragments and be clade-matched |
| A13 | Frame contrast against all five frames | NOT ADOPTED | Reverse frames hurt AT-rich, TGA-only fragments (T. thermophila -7.7 points) |
| S1 | Fix of the PWM start-refinement threshold | PARTIAL | Correct fix lowered accuracy: the old behaviour acted as a permissive threshold |
| S2 | Stabilise PWM training (k-fold, fixed window) | NOT ADOPTED | PWM starts are seed-dependent; the PWM mostly learned the coding/UTR boundary |
| S3 | Coding gate on start moves | NOT ADOPTED | Blocks correct moves where 5' UTRs look coding (athal 76%) |
| S4 | Trim open 5' ends by frame contrast | ADOPTED | Most of the exact-start gain comes from trimming into UTRs, not from start context |
| S5 | Choose among in-frame ATGs by PWM | NOT ADOPTED | Later ATGs by context help only C. reinhardtii |
| S6 | Kozak term: multiplicative vs additive | ADOPTED (additive) | Additive helps fragments more; complete inputs barely change |
| S7 | alpha x w grid, leave-one-organism-out | ADOPTED | alpha trades complete vs fragment accuracy; default is a 50/50 compromise |
| S8 | Kozak in bits (w = 17/bit) | ADOPTED | Rounding the weight changes at most 0.13% of predictions |
| S9 | First-ATG-of-longest-ORF baseline | ADOPTED (reported) | Beats TD on complete exact F1 (0.956 vs 0.912 at release); fails on fragments and random sequence |
| S10 | UTR-length prior for keep-open vs trim | NOT ADOPTED | Short open UTRs are common in real annotation; the prior just shifts alpha's trade-off |
| S11 | Keep-open vs trim for 5'-open ORFs (E33) | ADOPTED (PyTD only, homology mode; classifiers not adopted) | Self-trained classifiers fail (benchmark artefact "ORF lacks a stop"; context does not transfer). Homology hits (subject N-terminus at the ATG -> trim; alignment or Pfam domain upstream -> keep) give +2 to +9 exact F1 points on simulated assemblies; in PyTD a5b24c3. Top-n hits, looser e-values (step 4) and UniRef50 (step 5) did not improve on top-1 Swiss-Prot at 1e-5 |
| S12 | Later annotated ATGs (C. reinhardtii) | OPEN | 19.6% of C. reinhardtii starts are later ATGs; no tested rule recovers them |
| S13 | One-sided Kozak term ("Kozak floor", E35) | NOT ADOPTED | A negative Kozak score at the first ATG is right ~3/4 of the time on fragmented inputs (an internal ATG of a truncated CDS); flooring it helps complete inputs and loses more on fragments and assemblies |
| S14 | Per-transcript alpha from transcript / ORF length (E36) | OPEN (recommended to Brian: logO) | Length relates to truncation the same way on benchmark and assemblies; O/L does not (single-cut artefact, -15 points held out). alpha = 0.35 + 0.15 log2(O/1kb) passes E35 acceptance and raises every held-out assembly ab initio (+0.3 to +3.9) |
| O1 | Minimum ORF length 100 -> 50 codons | ADOPTED | Needs a coding test below 100 codons, else lncRNA calls triple |
| O2 | Fixed short-ORF guard (C = 200 below 100 codons) | PARTIAL | Kept short-protein gain at -m 100 noise; replaced at release by A9 |
| O3 | Post-hoc minimum length on GM12878 | NOT ADOPTED | 50 vs 100 aa barely changes precision; the step is at 100-150 aa |
| O4 | Transcript-level lncRNA filter | NOT ADOPTED | Classifiers reject lncRNAs better but cannot place CDSs in fragments; recommend combining tools |
| H1 | Homology retention on the benchmark | ADOPTED (optional mode) | Gains mostly on fragments; at release it doubles lncRNA calls |
| H2 | Clade-excluded homology series (GM12878) | ADOPTED (design) | Exclude the organism's own clade; Pfam is not distance-controlled |
| H3 | Exclude TE/Alu-derived homology hits | NOT ADOPTED | One Alu-containing protein (UTY_PANTR) drove thousands of retained ORFs |
| H4 | DIAMOND `--more-sensitive`, top-n hits, looser e-value (E33 step 4) | NOT ADOPTED | Best variant +0.1-0.2 exact F1 points held out, 2-5x slower, and lowers some complete inputs |
| H5 | UniRef50 instead of Swiss-Prot (E33 step 5) | NOT ADOPTED | Helps poorly covered clades (athal d16 +5 points) but human lncRNA calls 3,268 -> 8,220 and unreviewed orthologs keep UTR extensions open; 21x slower |
| B1 | bench_v2 (RefSeq) -> bench_v3 (current annotations) | ADOPTED | Do not quote bench_v2 numbers |
| B2 | GC-matched random and GENCODE lncRNA negatives | ADOPTED | Stop F1 on complete coding transcripts alone does not separate methods |
| B3 | Balanced mixed inputs (complete/5'/3'/internal) | ADOPTED | Single-type fragment inputs favour one start behaviour |
| B4 | Matching and scoring rules | ADOPTED | Our own `fp` column double-counted until fixed 2026-10-21 |
| B5 | Baselines (longest ORF, all ORFs) | PARTIAL | All-ORFs baseline scored but not reported |
| B6 | ANGEL as a comparator | ABANDONED | 17-22 h per A. thaliana run |
| B7 | CodAn | PARTIAL | Needs line-wrapped input; too slow for the simulated assemblies |
| B8 | Non-standard genetic code handling | ADOPTED | Five comparators cannot use the ciliate code: mark n/a, not "lost" |
| B9 | Old-vs-new critique tests; bootstrap CIs | OPEN | CIs are narrow (+-0.002-0.006); with 8 organisms use sign counts |
| B10 | Release-history benchmark | ADOPTED (supplement) | r2012 coordinates need clipping before scoring |
| T1 | Simulated GM12878 Trinity assembly | ADOPTED | Real-assembly truth is mostly unresolvable (91%); simulate instead |
| T2 | Labeller end tolerance (10 nt) | ADOPTED | Decided after a strict pass; report both passes |
| T3 | Read-depth ablation | ADOPTED | Depth controls fragmentation; d4-d8 match the real Trinity assembly |
| T4 | Simulated Arabidopsis Trinity assembly (E32) | ADOPTED (in paper: Results, Methods, Suppl. Fig. S7; manuscript e3b134d) | Few expressed noncoding transcripts (0.4% of reads): a frame/end-placement test. Ab initio TD leads to d8, longest ORF passes it at d16-d32; homology + start assist leads at every depth (0.899-0.915) |
| T5 | Trinity de novo on GM12878 | ABANDONED | Projected to take weeks on this host |
| T6 | Anatomy of missed CDSs; keep/trim errors on 5'-open 3' partials (E34) | OPEN (diagnostic) | Misses are short, weakly coding CDSs failing the q test; relaxing it costs ~1 lncRNA call per CDS; kept-open UTRs score coding, not a first-ATG issue |
| G1 | Three GM12878 assemblies vs MANE/GENCODE | ADOPTED | Sensitivity is limited by assembly, not CDS calling |
| G2 | Ribo-seq support (Cenik 2015) | ADOPTED | A+T1 data: weak periodicity, so "covered" is more informative than the frame test |
| G3 | MS peptide support (Sinitcyn 2023) | ADOPTED | Peptides on "novel" ORFs mostly come from paralog/pseudogene copies |
| G4 | Expressed-gene sensitivity | ADOPTED | Gene-level TPM suits MANE, not all GENCODE isoforms |
| G5 | Protein-match scoring on GM12878 | ADOPTED | Lets rules be compared on real data without a CDS truth set |
| G6 | Release reruns on GM12878 | ADOPTED | q-value test removes more catalogued Ribo-seq small ORFs than the length rules |

Outcome counts: see the end of this file.

---

## 1. Algorithm and coding test

### A1. Perl/PyTD parity fixes (2026-10-17 onward)
- **Question:** do the Perl and Python implementations give the same output?
- **Tried:** regression on the 5 TransDecoder sample sets, comparing LongOrfs output, hexamer model,
  6-frame scores, pre-refinement selections and final ORFs; then fixes in PyTD (background base
  frequencies parsed wrongly, training tie order, minimum length for 3'-partials, ORF numbering,
  gene/CDS lines treated as candidates, genetic-code names) and in both (genome-propagation gene IDs,
  usage text).
- **Outcome: ADOPTED.**
- **Numbers:** before the fixes (v6.0.0 vs PyTD 43f4293), StringTie sample 22 Perl-only / 30 PyTD-only
  ORFs; after, 0 / 0 with only PWM start differences (12/53/97/4/0 per set). At development defaults,
  84/84 benchmark inputs identical; at release, 84/84 identical and decoys, q-values and predictions
  byte-identical on 7 parity cases. PyTD runs 2.71x faster (median, release).
- **Lessons:**
  - Check intermediates with a parser, not by bytes: float printing and FASTA case differ (Perl keeps
    soft-masked lowercase).
  - A random-number generator must be portable to get identical output: the release uses
    xoshiro128** seeded by FNV-1a in both languages.
  - The "10x faster" and "Perl 679 vs Python 682" claims in notes/02 predate the PWM port; do not use them.
- **Pointers:** E08, E24, E25; notes/03, notes/04; `$T/regression/`; `data/regression/`.

### A2. Regenerated GC-calibrated long-ORF table (2026-10-17)
- **Question:** the 2017 table (TransDecoder 33247f5) of random-ORF length quantiles had no surviving
  simulation code. Can it be reproduced and documented?
- **Tried:** seeded simulation of LongOrfs ORFs in i.i.d. random sequence, 5% GC steps, >= 85,000 ORFs
  per level; low-GC levels pooled from chunks.
- **Outcome: ADOPTED** in both implementations (PyTD f5b12e7, Perl 694e204).
- **Numbers:** q0.999 from 519 nt (25% GC) to 2,436 nt (80% GC); within ~3% of 2017 at >= 50% GC;
  2017 was 2.8-10.3% lower at 25-45% GC; bootstrap 95% CI half-widths 1-2%.
- **Lessons:**
  - In the 2017 table the q0.999 and q0.9999 values were nearly equal at low GC, a sign of a small sample.
  - The quantity is a per-ORF 0.999 quantile, not an FDR; the old usage text said "1% FDR".
  - Gaps: the 50-80% rows come from runs whose command lines were not saved, and the final table was
    merged by hand.
- **Pointers:** E10; notes/01 §6, notes/03, notes/04; `$T/analyses/long_orf_threshold/`;
  `data/long_orf_threshold/final_thresholds.tsv`.

### A3. GC estimate, code-aware length rule, deterministic ties (2026-10-18 to 10-20)
- **Tried:** GC = (G+C)/(A+C+G+T) from counts instead of 2 x freq(C); length retention only for codes
  whose stops are TAA/TAG/TGA; no length retention above 80% GC; deterministic ORF tie-break.
- **Outcome: ADOPTED** (PyTD efc1450, 48c870c; Perl f7e75a1, de3f57e).
- **Numbers:** athal complete GC 39.4% -> 42.17% (threshold 615 -> 663 nt); C. reinhardtii 64.4% ->
  68.29% (1,083 -> 1,347 nt). Perl/PyTD differing prediction lines on 7 inputs fell to 0 after the GC
  fix and the tie-break.
- **Lesson:** 2 x freq(C) assumes complement symmetry, which fails with `-S`. Because the table is
  looked up in 5% steps, a small change in the estimate can move the threshold a whole step.
- **Pointers:** E10; `benchmark_v2/launch_v3/launch15.sh`, `launch16.sh`.

### A4. Keep the length (long-ORF) fallback (2026-10-21; release update 2026-10-25)
- **Question:** true CDSs fail the coding test; does the length rule earn its place?
- **Tried:** classified every truth CDS on bench_v3 complete inputs (8 organisms) by which part of the
  test it fails, the cause, and whether the length rule or homology rescues it.
- **Outcome: ADOPTED** (rule kept; analysis in Results, Fig. 5).
- **Numbers:** 4.9% of 135,429 matched CDSs fail; 67% of CDSs < 300 nt fail (53% at release). The
  length rule rescues 447; 85 (19%) of those have no homology hit. Predictions kept by length alone
  are 0-1.4% of all; on human lncRNA input 49 true vs 159 on lncRNAs (release).
- **Lessons:**
  - Failure causes are GC3/codon-usage minorities (fly chr4/Y, low-GC3 Chlamydomonas, Tetrahymena
    ribosomal proteins), low-complexity proteins (KRTAP, LCE), and unconserved, possibly
    mis-annotated models.
  - The rule does not apply to T. thermophila (TGA-only code).
- **Pointers:** E12 (and its Release update); `analyses/long_orf_failures/`;
  `pipeline/derive/05_*_v7.*`; `data/long_orf_failures/v7/`; `notebooks/05_coding_failures.Rmd`.

### A5. Drop the S1 > 0 requirement (2026-10-21, exploratory)
- **Outcome: NOT ADOPTED.**
- **Numbers:** development rule: +6 to +537 true CDSs as candidates per input, while random transcripts
  with a candidate rose 154 -> 3,548 (human) and 467 -> 19,684 (T. thermophila). Release rule: true
  CDSs -0.66 to +0.65%, random sequences with a candidate 1.4- to 7.9-fold (median 2.5).
- **Lesson:** under the q-value test the cost is smaller, but there is still essentially no gain.
- **Pointers:** E12 (`alt_rules_negatives.tsv`).

### A6. GC3-matched and alternative training sets (2026-10-21)
- **Tried:** training the Markov model on all annotated CDSs, 500 random CDSs, 500 longest, or 500
  GC3-quintile-matched CDSs.
- **Outcome: OPEN.** The analysis is in Results, but no record says whether a GC3-stratified model was
  considered for the code.
- **Numbers (development):** failure rate (ORFs >= 300 nt) 2.38% with TD's model vs 1.95% with GC3-matched.
  GC3-matched models rescue 25% of score failures (69% in D. melanogaster); at release, 27% / 63%.
- **Lesson:** the training-set choice matters little except for GC3 minority classes.
- **Pointers:** E12 (`oracle_models.tsv`).

### A7. Contrast test on all ORFs as a `--precise` mode (2026-10-23)
- **Question:** Brian asked for a precision/recall curve for extra score thresholding to support a
  `--precise` option.
- **Tried:** 16 settings (short-ORF length limit x contrast factor, with or without length auto-accept)
  over 20 inputs, rerun from stored checkpoints; post-filters on the selected ORF (per-base contrast,
  total contrast, both); a minimum-length design at -m 50/98/150 on GENCODE lncRNAs.
- **Outcome: NOT ADOPTED.** Brian then asked for no `--precise` flag of this kind (E19). The final
  `--precise` is a q-value setting (A11).
- **Numbers:** `sall_c200` (contrast test at every length): lncRNA calls 18.9% -> 6.0% at -2 points stop
  sensitivity (complete) and -6 (fragments). GM12878: non-matching proteins -70% / -72% for -3.2% /
  -3.6% true positives.
- **Lessons:**
  - Raising the threshold for short ORFs alone saturates: lncRNA false positives are mostly 100-300-codon
    ORFs that pass the length gate.
  - Dropped GM12878 ORFs have about half the Ribo-seq frame support of kept ones (0.33-0.36 vs
    0.62-0.66), so some are real.
  - LNCipedia could not be downloaded (404 on 2026-10-23), so GENCODE lncRNAs were used instead.
- **Pointers:** E18; `benchmark_v2/scripts/stringency_sweep_pytd.sh`; `data/benchmark/stringency/`.

### A8. Smooth short-ORF threshold shapes (2026-10-23/24)
- **Question:** remove the 100-codon cliff in the fixed short-ORF test.
- **Tried:** taper (C x mu at 100 codons falling linearly to 0 at the auto-accept length), inverse
  (C x 100/L) and fixed (C at every length), 172 fresh runs; later an "oracle curve" from
  annotation-derived prior log-odds (64 runs).
- **Outcome: NOT ADOPTED.** Superseded by the q-value test (A9); the experimental rules were removed
  at release (PyTD be23241).
- **Numbers:** cliff: stop sensitivity 0.60 at 75-99 codons vs 0.91 at 100-124, with lncRNA calls piled
  up at 100-124 codons (4,971 vs 877 at 75-99). `taper_c150`: lncRNA calls 18.9% -> 10.7%, at -0.4 points
  stop sensitivity (complete) and -1 (fragments). The oracle curve was no better than the linear taper.
- **Lessons:**
  - The strength C matters more than the shape; all shapes lie on nearly the same trade-off curve.
  - The taper needs a fallback length (300 codons) where there is no auto-accept length, e.g. TGA-only codes.
- **Pointers:** E19, E20 (Earlier steps); PyTD branch `exp/short-orf-taper` (local);
  `benchmark_v2/scripts/short_orf_rule_run.sh`; `data/benchmark/rules/`.

### A9. Target-decoy q-value acceptance (2026-10-24; release 2026-10-24)
- **Question:** a single, length-free acceptance rule.
- **Tried:**
  - Decoy ORFs come from dinucleotide-shuffled copies of each transcript, found and scored like targets.
  - Each candidate gets a q-value for its frame contrast, stratified by ORF length x transcript length
    x complete/partial type.
  - Accept if S1 > 0, C > 0 and q <= cutoff. Prototype as a post-filter, then 126 fresh runs.
- **Outcome: ADOPTED** as the release rule (PyTD 3a793ff, Perl 41fc5ee).
- **Numbers** (release medians vs development fixed short-ORF test):

  | input | stop sensitivity | negatives called |
  |---|---|---|
  | complete + random | 0.964 vs 0.957 | 0.29% vs 0.48% |
  | fragments + random | 0.872 vs 0.815 | 0.23% vs 0.11% |
  | complete + lncRNA | 0.963 vs 0.976 | 5.0% vs 18.9% |
  | fragments + lncRNA | 0.858 vs 0.858 | 5.9% vs 18.7% |

  All 84 inputs: fragments exact F1 0.812 vs 0.784 (median).
- **Lessons:**
  - Stratifying by transcript length was the key ingredient (Brian's suggestion): with ORF-length strata
    alone the prototype only matched the taper. Fragments cap the possible ORF length.
  - Calibration differs by negative type: far below nominal on random sequence, about nominal or above
    on lncRNAs (see Inconsistencies).
  - Costs:
    - GM12878 catalogued Ribo-seq ORF stop matches fell from 178/217/455 (fixed test) to 101/95/186
      (LRAA/StringTie/Trinity).
    - MANE exact given assembled fell 0.5-1.2 points, for 2.4- to 4.4-fold fewer unannotated CDSs.
  - Call lncRNA hits "predictions on annotated lncRNAs": only 1-2% end at a catalogued Ribo-seq ORF
    stop, but lncRNAs can encode peptides.
  - The decoy step makes PyTD runs 1.75x slower (median).
- **Pointers:** E20, E22, E24, E25; `pytransdecoder/core/decoys.py`; `data/benchmark/rules/`,
  `data/benchmark/final_v7/`; `notebooks/04_short_orfs_qvalues.Rmd`.

### A10. Long-ORF auto-accept under the q-value test (2026-10-24)
- **Tried:** none, the GC table (`gc`), or a decoy-derived length (0.999 quantile of decoy ORF lengths).
- **Outcome: ADOPTED (`gc`).**
- **Numbers:** without auto-accept the decoy rule loses ~0.4% of long true CDSs (stop sensitivity at
  >= 400 codons 0.9950 vs 0.9982). `gc` restores it (0.9984). `decoy` gives 0.9991, but in
  T. thermophila (TGA-only) its 1,183 nt length raises random-sequence calls 43 -> 333.
- **Pointers:** E20.

### A11. Named q-value settings (2026-10-24)
- **Outcome: ADOPTED:** default q 0.05, `--precise` 0.01, `--sensitive` 0.2, `--max-qvalue`.
- **Numbers (release):** `--precise` lncRNA calls 4.25% / 4.56% at -1.3 / -6.4 points stop sensitivity
  vs default. `--sensitive` fragments + random 0.895 stop sensitivity; q 0.3 and 0.5 add almost nothing
  (0.898) while lncRNA calls climb to 15-20% at q 0.5, so q 0.5 is not recommended.
- **Lesson:** the remaining lncRNA calls at q 0.01 are long (human: 613 of 1,322 >= 200 codons) and
  hardly change with q; no score cutoff removes them without losing true CDSs.
- **Pointers:** E22, E25.

### A12. Presets for small or partial inputs (2026-10-24)
- **Tried:** held-out test (MD5 split by ID) of training on the input vs same-species,
  same-species "mix" (complete + fragments) and other-species presets, at 50-5,000 sequences; 912 fresh runs.
- **Outcome: ADOPTED.** Seven bundled presets (mammal, insect, plant, fungus, apicomplexan, ciliate,
  green_alga), built from complete transcripts plus fragments; a warning below 1,000 transcripts.
- **Numbers:**
  - Training on the input with 50-500 sequences (half random) calls ORFs on 8-82% of random sequences;
    at 1,000 it works again (0.8-0.9% complete).
  - Presets keep random calls at 0-0.4%.
  - Complete-only presets lose on fragments (0.806 vs 0.880); "mix" closes the gap (0.873).
  - Human -> Arabidopsis/Drosophila presets lose on fragments (0.72 vs 0.87).
- **Lesson:** the bundled presets are built from the benchmark annotations, so never score them on
  benchmark inputs; accuracy claims come from the held-out test.
- **Pointers:** E21, E24; `pytransdecoder build-preset`; `data/benchmark/presets/`.

### A13. Frame contrast against all five other frames (2026-10-25)
- **Question (Brian):** v5.0.0 dropped the reverse-strand frames from the contrast on a small benchmark.
  Does restoring them help under the q-value test?
- **Tried:** `--contrast-frames all` (C = S1 - max(S2..S6) for targets and decoys); 156 fresh runs plus
  GM12878.
- **Outcome: NOT ADOPTED** (forward frames kept; Methods now cites this result).
- **Numbers:**
  - At q 0.05, worse on 75 of 84 inputs (exact F1 median -0.33 points).
  - At matched negative rates the two are equal, except stranded fragments, where the release rule is
    ~1 point ahead.
  - T. thermophila fragments -7.7 points: 18% of candidates removed before the q test.
  - Run time 1.41x. GM12878: -306 to -2,226 matching proteins.
- **Lesson:** in an AT-rich genome with a single stop codon (TGA), the reverse complement of a coding
  fragment often scores above its own frame.
- **Pointers:** E26; PyTD branch `exp/contrast-frames` 83a4806 (local, not merged); `data/benchmark/frames/`.

## 2. Start codons

### S1. Fix of the PWM start-refinement threshold (2026-10-17)
- **Tried:** the v6.0.0 refinement learned its threshold on length-normalised PWM scores but applied it
  to raw sums, on a window shifted 1 nt. Both implementations were fixed (Perl 8cf8d5f, PyTD 0525f5a).
- **Outcome: PARTIAL.** The fix was released, then the PWM refinement itself was replaced by S4.
- **Numbers:** the fix lowered complete-transcript exact accuracy; revised starts, complete, fix vs
  v6.0.0: human 3,540 vs 4,207; T. thermophila 615 vs 4,593. Reverting it restored v6.0.0-level F1.
- **Lesson:** a "correct" fix can lower accuracy when the old behaviour happened to act as a permissive
  threshold, which suits complete transcripts and hurts truncated ones.
- **Pointers:** E09; notes/01 §5-6, notes/03, notes/04; `$P/work/analysis/fork_findings.md` (not in git).

### S2. Stabilising PWM training (2026-10-17/18)
- **Tried:** seeds 1-5; k-fold CV with a deterministic noise filter and several threshold objectives;
  a fixed 20/10 window; "first passing ATG".
- **Outcome: NOT ADOPTED** (superseded by S4).
- **Numbers:** starts differing from seed 1: PASA 29-61, StringTie 47-144 (stops identical). The best-AUC
  window swings, e.g. PASA (15,2) -> (20,2) -> (3,7). k-fold gave athal complete 23,444-23,530 exact vs
  legacy 22,384-22,815, but athal fragments fell to 15,110-15,311 (none 17,055).
- **Lessons:**
  - The PWM learned mainly the coding/non-coding transition, not the start context.
  - Run on bench_v2 (RefSeq) truth, now superseded.
  - The Trinity seed range in notes/03 (14-26) has no run directory and cannot be verified.
- **Pointers:** E09; notes/03, notes/04; `analyses/start_refinement_stability/`; `$P/work/start_grid/`.

### S3. Coding gate on start moves (2026-10-17)
- **Tried:** move the start only if the skipped segment scores <= 0 under the Markov model.
- **Outcome: NOT ADOPTED** (recorded as a negative result).
- **Numbers:** it removed most truncated-ORF errors but blocked 76% of correct moves in athal (16% in
  C. reinhardtii), because athal 5' UTR segments often score > 0.
- **Lesson:** a hard gate is too blunt; S4's graded threshold replaced it. The 76%/16% numbers exist
  only in notes/04 (no saved table).
- **Pointers:** notes/04; E09 Gaps.

### S4. Trim open 5' ends to the first ATG by frame contrast (2026-10-19)
- **Tried:** trim when the segment's coding score is below alpha x |u| x mu (alpha = 0.4); later with a
  Kozak term (S6).
- **Outcome: ADOPTED** (release rule: alpha 0.4, Kozak weight 17 per bit).
- **Numbers:** exact TP summed over 8 organisms (complete): no refinement 83,307; fixed PWM 90,163;
  v6.0.0 96,957; frame contrast 114,400; final 116,384. Median complete exact F1 0.687 -> 0.910 (release
  0.912).
- **Lessons:**
  - Most of the gain is from trimming open 5' ends into UTRs; the Kozak term adds about 1% on complete
    transcripts.
  - The `pytd_dev_v2` label also includes -m 50 and other changes, so fragment gains are not
    attributable to the start rule alone.
- **Pointers:** E09, E13; `figures/fig3_start/`; `notebooks/03_start_refinement.Rmd`.

### S5. Choose among in-frame ATGs by PWM (2026-10-19)
- **Outcome: NOT ADOPTED.**
- **Numbers:** `pytd_dev_pwm` vs `pytd_dev`, complete: lower in 6 organisms, higher in C. reinhardtii
  (12,661 vs 12,592), equal within 1 in P. falciparum.
- **Pointers:** E09.

### S6. Kozak term: multiplicative vs additive (2026-10-19/20)
- **Tried:** multiplicative (alpha + 0.15K) x L (PyTD 0886981) vs additive alpha x L + 25K per nat
  (PyTD f5e6a39).
- **Outcome: ADOPTED (additive).**
- **Numbers:** additive minus multiplicative exact F1: fragments +0.005 to +0.023; complete -0.006 to
  +0.015 (mean +0.001).
- **Lessons:**
  - Self-trained Kozak models agree with annotation-trained ones (r 0.84-0.98) and recover the -3 A/G
    preference in all 8 organisms.
  - The reasons for the additive form and weight 25 were first recorded only in commit messages;
    E13 later supplied the grid.
- **Pointers:** E09, E13; `figures/kozak/`.

### S7. alpha x w grid and leave-one-organism-out (2026-10-21)
- **Tried:** 78 (alpha, w) cells x 16 stranded inputs (complete and fragments weighted equally);
  leave-one-organism-out; GM12878 sub-grid.
- **Outcome: ADOPTED:** (0.4, 25 per nat) confirmed.
- **Numbers:**
  - Mean exact F1 0.8322 at (0.4, 25); eight cells within 0.5 points (alpha 0.2-0.5, w 10-50).
  - Leave-one-organism-out picks (0.4, 25) in 7 of 8 folds (T. thermophila: 0.3).
  - w 25 vs 0: +0.6 (complete), +2.7 (fragments) points.
  - Complete inputs alone prefer alpha 3; fragments prefer alpha 0-0.2.
- **Lesson:** alpha sets the complete-vs-fragment balance, so the "best" value depends on the input mix.
  The 50/50 objective itself was not cross-validated.
- **Pointers:** E13; `analyses/start_param_sensitivity/`; `data/start_param_sensitivity/`.

### S8. Kozak score in bits, w = 17 per bit (2026-10-22)
- **Outcome: ADOPTED** (PyTD 44e417a, Perl ab38ca6).
- **Numbers:** w = 17 vs the exact equivalent 17.33: 0-26 predictions differ per input (<= 0.13%),
  exact F1 within 0.0003; regression sample sets: 4 of 4,277 ORFs change.
- **Pointers:** E13 addendum, E08.

### S9. First ATG of the longest ORF as a start rule (2026-10-23)
- **Tried:** baseline `longest_orf_atg`: trim each 5'-open longest ORF to its first ATG.
- **Outcome: ADOPTED as a reported baseline.** On 2026-10-28 Brian approved one sentence in Results
  (benchmark). It was first OPEN: E07 said "not yet in the manuscript".
- **Numbers:** complete exact F1 median 0.942 (final defaults 0.910; higher in all 8 organisms; bootstrap
  diff +0.032, CI 0.028-0.037); fragments exact F1 0.433 (final 0.789).
- **Lesson:** on complete transcripts this simple rule wins; TransDecoder's case rests on fragments and
  noncoding sequence. That is why the manuscript presents the noncoding test first.
- **Release recompute (2026-10-28):** same definition, on the `pytd_v7` LongOrfs output (50-aa minimum).
  Medians over 8 organisms, release vs baseline:
  - complete exact F1 0.912 vs 0.956 (difference +0.045, 95% CI 0.042-0.050; higher in all 8 organisms);
  - fragments exact F1 0.822 vs 0.484 (lower in all 8);
  - complete + random stop precision 0.994 vs 0.550, with 80.8% of random sequences called (release 0.28%);
  - mixed inputs exact F1 0.843 vs 0.497.
- **Pointers:** E07 Caveats and Addendum (2026-10-28), E17; `benchmark_v2/scripts/longest_orf_with_start.py`;
  `pipeline/derive/02_first_atg_baseline.py`; `data/benchmark/final_v7/first_atg_baseline.tsv`.

### S10. UTR-length prior for keep-open vs trim (2026-10-25)
- **Question (Brian):** a very short implied 5' UTR may signal truncation; does a prior on UTR length help?
  Is the 5'-partial deficit a test-data artefact?
- **Tried:** error anatomy; annotation UTR statistics; a fixed annotation-derived likelihood-ratio
  prior (lambda = 1, set a priori); fresh runs plus offline simulation.
- **Outcome: NOT ADOPTED** (kept unmerged, PyTD branch `exp/utr-prior` c9d7875). The manuscript says so.
- **Numbers:**
  - Every stop-correct, non-exact 5' partial is a trim to an internal ATG, at a median 55 nt from the
    fragment start.
  - The prior changes complete exact F1 by +0.02 to +2.60 points, fragments by -0.60 to -1.30, and
    lowers all 24 mixed inputs.
  - GM12878 MANE exact given assembled +0.3 to +0.45 points.
  - Short-UTR-only variant: about 0.1 point.
- **Lessons:**
  - The deficit is real, not a cut-model artefact. But open 5' UTRs in real annotation are short
    (median 27-118 nt), so UTR length carries little information (|log2 LR| <= 1.1 over 10-150 nt).
  - The benchmark over-weights truncation: 76% of 5'-open ends in fragments and mixed inputs are truly
    truncated, against 20-37% in GM12878 assemblies.
- **Pointers:** E29; `pipeline/derive/03_utr_*`; `data/start_refinement/utr_prior/`.

### S11. Keep-open vs trim classifier for 5'-open ORFs (2026-10-28)
- **Question:** can rules or classifiers using |u|, C(u)/mu, Kozak, transcript length, ORF length,
  3' UTR length etc. decide better than the release rule? Do decoys give a UTR null?
- **Tried:**
  - Step 1: oracle-label ceiling with logistic regression and boosted trees under several train/test
    schemes, on 48 benchmark inputs and 7 simulated assemblies.
  - Step 2: self-trained variants without labels: synthetic cuts of confident ORFs, with |u| x
    transcript-length matching; Saerens EM prior adjustment; an unsupervised mixture over context.
  - No TransDecoder code was changed.
- **Outcome:** classifiers (steps 1-2) NOT ADOPTED and kept out of the manuscript; the homology rules
  (step 3) ADOPTED in PyTransDecoder only, in homology mode (Brian, 2026-10-28; PyTD a5b24c3, v0.2.0).
- **Numbers:**
  - With u-intrinsic features only, fitted models do not beat the rule outside their training domain
    (bench -> sim AUC 0.92-0.94 vs 0.93).
  - In-domain context models with true labels reach accuracy 0.93-0.94 vs 0.87 (human d1), +533 to +1,141
    exact ORFs per assembly: a ceiling, not achievable without labels.
  - Best self-trained overall (`b_rule_ctx`): +0.4 to +2.7 points on human assemblies, but 5'-partial
    F1 -14 to -34 points.
  - The only per-type-neutral variant (`a_u_match_noem`): +0.6 to +1.6 on assemblies, but -3.0/-3.5
    on S. pombe complete.
- **Lessons:**
  - **Benchmark artefact:** in `fragments` inputs every 5'-open UTR-extension candidate is a 3' partial
    (no stop) and every truncated one has a stop, so "ORF lacks a stop" separates the classes
    perfectly (context AUC 0.98-1.00). Fragment-input gains from context features mean nothing for
    real data.
  - **No transfer:** benchmark-trained boosted trees recognise the input type from context (long
    3' UTRs on complete inputs) and fail on the assemblies (hgb_full on d32 accuracy 0.583 vs rule 0.822).
  - Synthetic 5' cuts do not reproduce how truncation co-varies with 3' completeness and transcript
    length in real assemblies.
  - Prior adaptation through miscalibrated synthetic models is unreliable.
  - Decoys approximate the UTR null at |u| < 200 nt but understate its tail.
  - Untested next idea: the u-only matched model combined with the mixture's pi estimate.
- **Step 3, homology override (2026-10-28; ADOPTED in PyTD a5b24c3, see below):**
  - Rules on the longest-ORF DIAMOND (Swiss-Prot minus own organism, -k 1) / Pfam hits:
    - trim when the subject's N-terminus (sstart <= 3) aligns within 5 aa of the first ATG (R2);
    - keep open when the alignment starts >= 5 aa upstream of the ATG with sstart > 1 (R1);
    - keep open when a Pfam envelope starts upstream (R3);
    - otherwise the release rule decides.
  - Exact F1 vs release:
    - R2 + R1 + R3: human d1 / d8 / d32 +2.0 / +5.4 / +9.0, athal +2.0; every CDS type improves;
      5' partials +11 to +14.
    - R2 + R1: +1.9 / +5.3 / +8.7 / +1.5.
    - R2 alone: +1.3 / +1.9 / +1.7 / +0.5, and never negative on the 16 complete benchmark inputs
      (+0.02 to +3.4; precision >= 0.99).
    - The keep rules cost up to 0.3 (R1) or 0.6 (R3) points on complete inputs.
  - Lesson: homology does what self-training could not.
  - Next steps:
    - validate on GM12878 with full Swiss-Prot;
    - test a stricter keep rule (k = 20).
  - Pointers: E33 step 3; `data/start_classifier/step3/`; searches in `$P/work/start_classifier/homology/`.
- **Adoption (2026-10-28):**
  - PyTD main a5b24c3 adds the step-3 policy H to start refinement:
    - BLAST rules only with `--retain-blastp-hits`, the Pfam rule only with `--retain-pfam-hits`;
    - on by default whenever hits are supplied; opt-out `--no-homology-start-assist`;
    - per-ORF table `start_refinement.homology.tsv`.
  - Checks:
    - ab initio output byte-identical to 3a793ff (5 sample sets + 4 benchmark inputs);
    - in-program decisions = offline H on 150,356 benchmark and 42,043 simulated ORFs (0 mismatches).
  - Fresh runs with reused searches. Exact F1 vs homology for selection only:
    - simulated human d1 / d8 / d32 / athal +2.0 / +5.4 / +9.0 / +2.0, every type up;
    - benchmark fragments +5.9 (median, 34 inputs, none down);
    - benchmark complete inputs median ~0, range -0.58..+2.92;
    - GM12878 Trinity: MANE exact 7,063 -> 7,190; 5' partials kept open given the stop 0.78 -> 0.96.
  - Stops and selection are unchanged. Perl is unchanged, so Perl = PyTD now holds for ab initio runs only.
  - Pointers: E33 "Adoption in PyTD"; `data/start_classifier/adoption/`.
- **Arabidopsis depth series (E32 levels d1-d32, 2026-10-28):**
  - Full homology mode (DIAMOND vs Swiss-Prot without Arabidopsis, Pfam, on all LongOrfs peptides) vs ab initio:
    - exact F1 +3.2 (d1) to +8.0 (d32) points: 0.914 / 0.915 / 0.909 / 0.901 / 0.899 / 0.908;
    - of this, retention alone (assist off) gives +0.8 to +1.0, and the start assist gives +2.3 to +7.2.
  - Best overall method at every level.
  - Ab initio reproduces E32's TransDecoder row exactly.
  - Pointers: E33 "Arabidopsis depth series"; `data/start_classifier/adoption/athal_depth_scores.tsv`.
- **Step 4, hit depth and e-value (2026-10-28; NOT ADOPTED, a5b24c3 settings kept):**
  - Question (Brian): is the top DIAMOND hit at 1e-5 (first HSP) too stringent, and what does the assist still get wrong?
  - Error anatomy of the assist-on runs (human d1 / d8 / d32, athal d1, and athal depth-series B d1):
    - 85-90% of the non-exact 5' partials and 83-89% of internal CDSs in human have no prediction with the true stop (nearly always none on the transcript), so starts cannot fix them;
    - the remaining decision errors are mostly a hit that fires no rule, because the alignment starts inside the subject (sstart > 3);
    - wrong R2 trims: 4-12 per assembly;
    - 3' partials are candidates only when the 5' UTR has no in-frame stop (28-43% of them); their errors are keep calls;
    - complete CDSs left open: R1 / R3 on homologs with longer N-termini (human), non-firing hits (athal 784 of 974).
  - Deep searches of the 5'-open candidates (DIAMOND `-k 25 --max-hsps 0 --evalue 10`, default and `--more-sensitive`; hmmsearch `--domE 1` with `-Z 1 --domZ 1`):
    - 4 assemblies and 16 benchmark inputs (complete + fragments, 8 organisms).
  - 3,030 offline rule variants per dataset: e-value 1e-10..10 x mode; top 1 vs top 3-25 with majority / unanimous / bitscore / any-upstream-unless-trim / first-firing; best vs all HSPs; Pfam GA vs domE thresholds.
    - Chosen on human d1 + d8 only.
  - Numbers (exact F1 points vs a5b24c3):
    - best on tuning (more-sensitive, e <= 10, all HSPs, top-3 bitscore-weighted): +0.09 tuning, +0.15 to +0.17 d32, +0.10 to +0.18 athal, +0.17 to +0.18 bench fragments (median), but complete inputs 4 of 8 down (tthermophila -0.24 with GA, -1.17 with domE);
    - top 3-5 hits on the default search: <= +0.04; top 10-25: negative;
    - e-value 1e-5 -> 10 in default mode: no change; 1e-10: -0.01 to -0.03;
    - all HSPs: no change;
    - Pfam looser than GA: down to -8.9 on complete inputs;
    - oracle ceiling +1.6 to +6.6.
  - Deeper hits give uncontested right-call support to only 4-22% of the wrong decisions, and contradict 6-20% of the right ones.
  - Lessons:
    - the remaining ceiling is in selection (missed CDSs) and in homologs that do not reach the subject N-terminus, not in hit depth;
    - more hits add more conflict than signal.
  - Pointers: E33 "Step 4"; `data/start_classifier/step4/`; scripts `pipeline/start_classifier/step4_*`; searches in `$P/data/start_classifier/step4/searches/`.
- **Step 5, UniRef50 instead of Swiss-Prot (2026-10-29; NOT ADOPTED, Swiss-Prot stays recommended):**
  - Question (Brian): would UniRef50 (or UniRef90) work better than Swiss-Prot for homology mode, especially for
    the start assist?
  - Setup:
    - UniRef50 2026_03 DIAMOND db with the representative's taxid. The FASTA `TaxID=` is the cluster's LCA, so
      the representative's taxid was taken from the XML; NCBI-merged taxids had to be resolved, or
      `--taxon-exclude` silently drops their clusters.
    - Organism taxon and its descendants excluded.
    - Searches: `-k 1` and `-k 25` at 1e-5 of the 5'-open candidates of 21 inputs; `-k 1` on all LongOrfs of 21
      full-mode inputs, including `.wRand` and `.wLnc`.
    - Offline rule variants (top-n, consensus N-terminus, subject-begins-with-M), chosen on human d1 + d8;
      fresh a5b24c3 Predict runs in full mode.
  - Numbers:
    - candidates with a hit: athal / dmel 49-59% -> 98-99%, tthermophila 18-27% -> 72-81%, human 93-97% -> 98-100%;
    - start decisions with the a5b24c3 rules (exact F1 vs Swiss-Prot):
      - human / mouse -1.1 to -6.9;
      - athal +0.7 / +3.8 (d1 / d16);
      - benchmark fragments of the sparse organisms +0.8 to +3.2;
    - full mode:
      - human -4 to -7, and -12 to -15 with lncRNA negatives (lncRNA calls 3,268 -> 8,220);
      - athal, dmel, tthermophila +0.7 to +5.1;
      - random negatives unchanged.
  - Mechanisms:
    - unreviewed mammalian orthologs with ~100 aa longer N-termini become the top hit and keep UTR extensions
      open (wrong keeps 394 -> 1,659 in human d1);
    - lncRNA ORFs gain hits to predicted proteins.
    - Not fragments: an M filter changes nothing.
  - Variants:
    - top-25 majority cuts the human loss to -0.3 and gives +2.9 / +4.8 on athal;
    - the consensus N-terminus projection is worse.
    - Post hoc hybrid (Swiss-Prot hit first, UniRef50 only when Swiss-Prot has none): isolate +0.1..+0.3 on
      human and most of the UniRef50 gain elsewhere, but in full mode it still adds the lncRNA false positives.
  - Cost:
    - DIAMOND 21x slower per input (13.4 min vs 37.5 s, 8 threads); db 13.9 GB.
    - UniRef90 estimated at ~51 GB and ~50 min per input; not worth testing now.
  - Lessons:
    - more sequences give coverage where Swiss-Prot is sparse, but unreviewed N-termini and predicted proteins
      cost start accuracy and specificity where it is not;
    - excluding UniRef clusters by taxon also removes the organism's orthologs (a harsher test than
      sprot_minus);
    - retention and the start assist would need separate hit files to use UniRef50 safely.
  - Pointers: E33 "Step 5"; `data/start_classifier/step5/`; scripts `pipeline/start_classifier/step5_*`; db and
    searches in `$P/data/start_classifier/step5/`.
- **Human depth series (2026-10-30; in the manuscript with Arabidopsis):** fresh A ab initio / B homology + start assist /
  C homology without assist on E30/E31 d1-d32 (sprot_minus_9606 + Pfam on all LongOrfs, PyTD a5b24c3). A = E31 release
  rows (byte-identical predictions). Exact F1 B 0.882 (d1) -> 0.949 (d32) vs A 0.844 -> 0.806 (+3.8 to +14.3 points;
  assist alone +2.5 to +11.1); best at every level. Cost: lncRNA contigs with a prediction 1,688 -> 3,319 (d1), 264 -> 361
  (d32). Pointers: E33 "Human depth series"; `data/start_classifier/adoption/human_depth_scores.tsv`; scripts 20-22.
- **Pointers:** E33; `pipeline/start_classifier/`; `data/start_classifier/{step1,step2,step3,step4,step5}/`; large
  intermediates in `$P/data/start_classifier/`; env `$P/env/startclf`.

### S12. Later annotated ATGs; other open start issues
- **Outcome: OPEN.**
- In C. reinhardtii the true start is a later in-frame ATG in 19.6% of transcripts (0.7-5.0%
  elsewhere); no tested rule recovers most of them (E09). E33 found ~30% of C. reinhardtii UTR
  candidates have a later true ATG.
- The GM12878 `start_downstream_in_tx` label (1,931 LRAA, 1,127 StringTie) is documented as "should not
  happen" and has not been investigated (E09 Caveats).
- **Pointers:** E09, E33.

### S13. One-sided Kozak term ("Kozak floor") (2026-10-30)
- **Question:** E34 found that 19-42% of the UTR extensions wrongly kept open are "Kozak vetoes" (u looks like UTR,
  but K < 0 at the first ATG pushes 0.4|u| + 17K below coding(u)). Does limiting the negative Kozak term help?
- **Tried:** offline, on 90 finished runs (benchmark 48 ab initio + 24 homology inputs; human d1/d8/d32 and athal
  d1-d32, ab initio and with the homology start assist), with the runs' own selections and scores:
  V1 17 max(K, 0); V2 floor at -0.5 / -1; V3 negative K weighted 8.5 / 4 per bit; V4 floor only when coding(u) <= 0.
  Acceptance set beforehand: mean complete-input change >= -0.1 point, no benchmark input below -0.5; E13 objective
  reported. V0 reproduced every official score (328 comparisons, 0 differences).
- **Outcome: NOT ADOPTED** (recommendation; no code changed).
- **Numbers (exact F1 points vs release):**
  - Ab initio: complete inputs +0.34 (V2_1) to +1.41 (V1) on average, every one up; fragments -0.80 to -3.59, every
    fragments and mixed input down (21-32 of 48 inputs below -0.5); E13 objective -0.12 (V3_8.5) to -1.09 (V1).
  - Assemblies ab initio: human d32 -1.3 to -6.1, athal d32 -2.1 to -8.8; only athal d1 (+0.03..+0.21, all but V2_1) and human d1
    with V3_8.5 (+0.04) gain. 5'-partial and internal F1 fall 1-14 points.
  - With the assist: V2_1 and V3_8.5 pass the benchmark criteria (E13 objective +0.06 / +0.15) but lower athal
    d8-d32 by 0.23-0.81 and 5'-partial / internal F1 on every athal assembly (-0.9 to -2.0).
  - Decisions changed: on complete inputs only gains (r2w = 0; 15-28% of the changes trim to the first ATG while the
    true start is a later ATG, neutral for exact F1); on fragments and assemblies the keep-open was right in 76-83%
    (ab initio) / 63-70% (assist) of the changed decisions.
- **Lessons:**
  - E34 counted only the UTR side of the veto. A poor first-ATG context mostly means an internal ATG of a truncated
    CDS; the vetoed UTRs and the protected truncated CDSs have similar |u|, coding(u) and K, so no u-intrinsic floor
    separates them.
  - Any floor moves the rule along the same complete-vs-fragment axis as alpha (E13/S7).
- **Pointers:** E35; `data/kozak_floor/`; scripts `pipeline/kozak_floor/`; per-candidate tables in `$P/data/kozak_floor/`.

### S14. Per-transcript alpha from transcript and ORF length (2026-10-30)
- **Question (Brian):** set alpha_i = g(L, O) per transcript (L = transcript length, O = candidate ORF length, or O/L,
  |u|/L) so that complete transcripts get a high alpha and fragments a low one.
- **Tried:** offline, on 95 finished runs. The E35 harness and candidate tables were reused, plus human d2/d4/d16 and
  GM12878 Trinity ab initio and with the assist.
  - Step 1: descriptive tables of the oracle keep/trim answer by L, O, O/L and |u|/L.
  - Step 2: 7 forms (const, logL, logO, OL, uL, stepOL, stepL; 7,220 configurations), tuned on the benchmark only
    (complete, complete_fullutr, fragments and mixed; E13-type 50/50 objective), with leave-one-organism-out selection,
    free or under the E35 acceptance criteria.
  - Step 3: held-out evaluation on the simulated human d1-d32 and athal d1-d32 assemblies and the GM12878 E29
    statistic, ab initio and with the assist.
  - Step 4: oracle ceilings per (L, O/L) bin.
  - Exact fast rescoring, checked against E35 (90 runs) and full substitution (0 difference).
- **Outcome: OPEN** (logO recommended; decision with Brian; no code changed).
- **Numbers (exact F1 points vs release):**
  - logO = 0.35 + 0.15 log2(O/1000): same parameters in all 8 LOO folds.
    - Benchmark ab initio: E13 objective +0.43, complete +0.08, fragments +0.78, mixed +0.84; 0 of 48 inputs < -0.5.
    - Held out ab initio: human +0.31..+2.65, athal +0.39..+3.93.
    - With the assist: +0.07..+1.57 (athal d1 B -0.03).
    - GM12878: 5' partials kept open +2.0, complete exact given stop -0.1.
  - logL (0.35 + 0.1 log2(L/1000)) and uL (0.5 - |u|/L) also pass; stepL passes the criteria but not "both up".
  - O/L forms:
    - ab initio benchmark-tuned: -1.4 to +0.3 on the assemblies;
    - tuned with the assist: they pass LOO (+0.40-0.53), then lose 1-15 points on the assemblies.
  - Ceiling (alpha per L x O/L bin, in-sample): human +1.2..+3.6, athal +2.7..+5.0 ab initio; at depth mostly a lower
    constant alpha. The per-candidate oracle is +3.6..+14.7.
- **Lessons:**
  - Short transcripts and ORFs are mostly fragments on every data set, so a length slope transfers.
  - O/L >= 0.95 means "no stop". That class is all UTR in the single-cut benchmark fragments but mostly truncated
    in assemblies. LOO by organism cannot catch a design artefact shared by all organisms.
  - The gain is mostly on internal and 5'-partial CDSs; the complete-input gain is small and organism-dependent.
- **Pointers:** E36; `data/alpha_by_length/`; scripts `pipeline/alpha_by_length/`; per-candidate tables in `$P/data/alpha_by_length/`.

## 3. Short ORFs and minimum length
- **Addendum (2026-10-30): simpler forms** (BH asked for a more intuitive formula; E36 addendum).
  - A two-value step (0.15 if O < 900 nt else 0.5) and a proportional ramp (0.5 min(1, O/1200)) recover about
    70-90% of logO's gain, with less stable parameters across folds.
  - A one-parameter linear form (0.225 O/kb) fails acceptance.
  - logO is equivalent to alpha = 0.5 log10(O/200 nt): 140 of 817,983 decisions differ.

### O1. Minimum ORF length 100 -> 50 codons (2026-10-20)
- **Outcome: ADOPTED** (LongOrfs default 50 at PyTD 429a89a), always with a coding test for short
  ORFs (O2, then A9).
- **Numbers:** human complete + lncRNA: 50-99 aa sensitivity 0.31 -> 0.86, but lncRNAs with a call
  6,608 -> 21,416 of 33,668 and precision 0.74 -> 0.47; ~35 lncRNA calls per correct call at 50-99 aa.
- **Lesson:** random sequence is a much easier negative than lncRNA (0.013-0.159 calls per random
  sequence at -m 50).
- **Pointers:** E04; notes/05; `data/benchmark/short_orf_summary.tsv`, `by_length.tsv`.

### O2. Fixed short-ORF guard: contrast >= 200 x mu below 100 codons (2026-10-20)
- **Tried:** guard factor 150 / 200 / 250 at -m 50.
- **Outcome: PARTIAL.** In the development defaults (PyTD 7e9ffa2 / 48c870c), removed at release in
  favour of A9.
- **Numbers:** human complete: short sensitivity 0.31 -> 0.56, lncRNA calls 6,608 -> 6,943, precision
  0.737 -> 0.729. GM12878 catalogued Ribo-seq ORF matches 458 -> 178 (LRAA) vs -m 50 without the guard.
- **Lessons:**
  - The guard was tuned coarsely: only 3 factors, the 100-codon limit was never varied, and there was
    no cross-validation.
  - Composition-based tests remove most translated noncanonical sORFs: precise annotation of canonical
    short proteins and capture of the "dark proteome" conflict.
- **Pointers:** E04; notes/05.

### O3. Post-hoc minimum length on GM12878 predictions (2026-10-21)
- **Outcome: NOT ADOPTED** (analysis only).
- **Numbers:** 50 vs 100 aa barely changes the GENCODE-exact share (LRAA 0.504 vs 0.506); at 150 aa
  unannotated predictions fall to 1-4% but MANE exact given assembled falls to ~0.85.
- **Lesson:** a post-hoc filter differs from a different `-m`, because training and per-transcript
  selection see the full candidate set.
- **Pointers:** E04; notes/06; `genome_annotation/scripts/minaa_sweep.sh`.

### O4. A transcript-level lncRNA filter inside TransDecoder (2026-10-20)
- **Outcome: NOT ADOPTED** (Brian, 2026-10-20: Discussion only; recommend pairing TransDecoder with a
  coding-potential classifier).
- **Numbers (development defaults):** lncRNAs called: TransDecoder 20.6%, CPAT 10.5%, CPC2 7.9%,
  RNAsamba 9.1%; but on fragments TransDecoder exact F1 0.68 vs <= 0.47.
- **Lesson:** RNAsamba calls 6.2% of lncRNA genes present in GENCODE v28 vs 10.9% of later ones,
  consistent with training leakage. Era-split negatives are worth keeping.
- **Pointers:** notes/05; E04 (`lnc_by_era.human.tsv`).

## 4. Homology

### H1. Homology retention on the benchmark (2026-10-18 to 10-26)
- **Tried:** DIAMOND vs Swiss-Prot without the organism's own taxid plus Pfam `--cut_ga`; old code
  (27bfb5e), development defaults (84 inputs) and release (84 inputs).
- **Outcome: ADOPTED** as the existing optional mode (no change to retention).
- **Numbers (release):** fragments: stop sensitivity +0.42 to +9.20 points, exact F1 +0.11 to +4.65;
  complete-type: exact F1 -0.35 to +2.09. Human lncRNAs with a prediction 1,639 -> 3,268.
- **Lessons:**
  - Homology-supported ORFs bypass the q-value test, so once that test removes chance ORFs, homology
    adds a large relative share of lncRNA calls.
  - Which lncRNA ORFs gain hits was not examined (OPEN).
  - Five T. thermophila runs were re-predicted at a pinned commit after the shared tree moved during
    the run (see I-section).
  - The per-taxid Swiss-Prot databases have no committed builder.
- **Pointers:** E07, E25; `data/benchmark/homology_vs_abinitio.tsv`; `homology_vs_abinitio.dev_v2.tsv`.

### H2. Clade-excluded homology series on GM12878 (2026-10-21)
- **Tried:** Swiss-Prot minus human / primates / mammals / vertebrates / Metazoa, BLAST only, Pfam
  only, and both; 33 configurations.
- **Outcome: ADOPTED** (evaluation design; in Results).
- **Numbers:** MANE exact given assembled, minus human: +1.8 / +1.1 / +1.4 points (LRAA / StringTie /
  Trinity); the gain shrinks with distance (Trinity 0.9196 -> 0.9095 at minus Metazoa vs 0.9055 ab initio).
- **Lesson:** Pfam alone recovers most of the gain, but Pfam HMMs include human sequences, so it is not
  distance-controlled.
- **Pointers:** E06; `genome_annotation/scripts/homology_series*.sh`; `data/gm12878/gm12878_homology_series.tsv`.

### H3. Screen out transposable-element homology hits (2026-10-21)
- **Outcome: NOT ADOPTED** (Brian, 2026-10-21: keep TE-homology-retained ORFs as candidate neo-ORFs;
  users can screen them by best hit or repeat overlap).
- **Numbers:** in Trinity, 4,692 of 4,981 CDSs found only with minus-human (not minus-primates) are
  unannotated; 4,465 have chimp UTY (Q6B4Z3, Alu-encoded C-terminus) as best hit; 95% are < 100 aa.
- **Lesson:** a single curated protein with a TE-derived segment can dominate homology retention; check
  best-hit concentration whenever homology adds many unannotated ORFs.
- **Pointers:** E06; `data/gm12878/gm12878_homology_primate_retention.txt`.

### H4. DIAMOND sensitivity, hit depth and e-value
- **Outcome: NOT ADOPTED (2026-10-28).** Tested for the start assist in E33 step 4 (see S11): `--more-sensitive`,
  top 3-25 hits with voting schemes, all HSPs, e-values 1e-10 to 10. The best variant gained +0.1-0.2 exact F1
  points held out, was 2-5x slower, and lowered some complete benchmark inputs. Loosening the e-value in default
  mode changed nothing. Pfam below GA thresholds hurt (-1.6 median on complete inputs). Not tested for
  retention alone (E06 series used default mode).

### H5. UniRef50 instead of Swiss-Prot
- **Outcome: NOT ADOPTED (2026-10-29).** E33 step 5 (see S11). UniRef50 raises 5'-open candidates with a hit
  from 49-59% to 98-99% in athal and helps clades Swiss-Prot covers poorly (athal d16 full mode 0.899 -> 0.950),
  but in human and mouse it is worse: unreviewed mammalian orthologs with longer N-termini keep UTR extensions
  open, and human lncRNAs with a prediction rise 3,268 -> 8,220 (1,639 ab initio). DIAMOND is 21x slower
  (13 min vs 37 s per benchmark input). Swiss-Prot stays the recommended database; UniRef90 not tested. The
  UniRef50 database was deleted in the 2026-10-29 disk cleanup (rebuild: `pipeline/start_classifier/step5_*`).

## 5. Benchmark design

### B1. Truth sets: bench_v2 (RefSeq) -> bench_v3 (current annotations) (2026-10-17/18)
- **Tried:** RefSeq `rna.gbff` truth (bench_v2, first runs on UGER) replaced by each organism's current
  reference annotation: MANE v1.5, GENCODE M39, FlyBase r6.69, TAIR12, PomBase, TGD v6, RefSeq for
  P. falciparum and C. reinhardtii. Eight organisms span 22-68% GC and the ciliate code; UTRs are
  trimmed to 250 nt, with a full-UTR set kept.
- **Outcome: ADOPTED.**
- **Lessons:**
  - bench_v2 numbers (e.g. notes/04 "TD 0.96-0.985 vs Prodigal 0.92-0.94") must not be quoted.
  - The strategy grid and the early start-fix comparison used bench_v2 truth.
  - C. reinhardtii v6.1 needs a JGI login, so v5.5 was used.
- **Pointers:** E07; `benchmark_v2/`; notes/02 (2018 design), notes/04.

### B2. Negatives: GC-matched random and GENCODE lncRNAs (2026-10-18 / 10-20)
- **Outcome: ADOPTED.**
- **Numbers:** on stranded complete transcripts the longest-ORF baseline has stop F1 0.976 (median, vs
  0.978 final defaults); with random negatives its precision falls to 0.803 (0.499 T. thermophila)
  vs 0.993.
- **Lessons:**
  - Without negatives, stop-level accuracy does not separate methods, which is why Results presents
    the noncoding test first.
  - lncRNA negatives are not perfectly clean; random sequence is the clean negative.
- **Pointers:** E04, E07.

### B3. Balanced mixed inputs (2026-10-25)
- **Question (Brian):** self-training methods could be skewed by an input of one partial type; the
  `fragments` inputs have only 5' and 3' partials cut near the CDS middle.
- **Tried:** 25% each of complete, 5' partial, 3' partial and internal, 3 seeded trials x 8 organisms.
  Earlier single-type `fragments_internal` / `fragments_5p` inputs were abandoned and moved to `_trash`.
- **Outcome: ADOPTED** (Table 1, Table S9).
- **Numbers:** TransDecoder stop F1 0.970, exact F1 0.843 (v6.0.0: 0.916 / 0.768); CodAn best on
  5' partials (exact F1 0.971).
- **Lessons:**
  - Inputs are named `fragments_*` on purpose, so that `run_tool.sh` picks partial models for CodAn
    and RNAsamba.
  - `score_all.py` walks every run directory: filter out `fragments_mixed_t*` before a full rescore.
- **Pointers:** E28; `pipeline/benchmark/make_mixed_inputs.py`, `score_mixed.py`; `data/benchmark/final_v7/`.

### B4. Matching and scoring (2026-10-18 onward)
- **Outcome: ADOPTED:** stop within 2 nt; exact = start within 2 nt as well; one credit per truth CDS
  (from the 2018 evaluation).
- **Lessons:**
  - `score_predictions.py` counted each non-duplicate FP twice in the `fp` column until 2026-10-21;
    sensitivity, precision and F1 were unaffected. Do not quote `fp` from older tables.
  - "Predicted length" bins include the stop codon.
- **Pointers:** E07 Caveats, E04 Caveats.

### B5. Baselines (2026-10-18 / 10-23)
- **Outcome: PARTIAL.** The longest-ORF baseline is reported. The all-ORFs (>= 100 codons) baseline was
  scored but not reported (Brian, 2026-10-23). The first-ATG baseline is S9 (reported; release recompute 2026-10-28).
- **Lesson:** on mixed inputs and simulated assemblies the longest-ORF baseline comes from `pytd_v7`
  LongOrfs, not the older `pytd` run.
- **Pointers:** E07, E22, E28.

### B6. ANGEL as a comparator (2026-10-19)
- **Outcome: ABANDONED** (stopped with Brian's agreement).
- **Numbers:** 3 completed A. thaliana runs took 17.4-21.9 h each (runtime 62,503 s vs minutes for
  others); 3 more were partial.
- **Lesson:** ANGEL's fallback reimplements TransDecoder's longest-ORF + hexamer scheme; it is
  deprecated since 2021.
- **Pointers:** E07 Gaps, E27; `methods_benchmark.tex`.

### B7. CodAn (2026-10-18 to 10-27)
- **Outcome: PARTIAL.** CodAn is in the main benchmark and the mixed inputs, but was not run on the
  simulated assemblies (too slow).
- **Numbers:** athal complete runtime 35,114 s; T. thermophila runs took days.
- **Lessons:**
  - CodAn needs line-wrapped input (all bench_v2 CodAn results were invalid).
  - Its env needed a `libnsl.so.1` symlink.
  - 5 of our 8 organisms are CodAn training species.
  - P. falciparum (22% GC) is below its lowest GC band: stop F1 0.710.
- **Pointers:** E07, E16, E27, E28; notes/04; memory `transdecoder-paper-project.md`.

### B8. Non-standard genetic code handling (2026-10-18 / 10-25)
- **Outcome: ADOPTED:** GeneMarkS-T, CodAn, CPC2, RNAsamba and Borf are n/a for T. thermophila, because
  none has a user-settable code. GeneMarkS-T was run with code 1 (it supports 1/4/11 only).
- **Pointers:** E07, E27.

### B9. Earlier critiques old vs new; bootstrap CIs (2026-10-23)
- **Outcome: OPEN** (planned supplement; tables marked interim pending CodAn).
- **Numbers:**
  - v5.7.1 vs new, complete: non-ATG predictions 15.7% -> 5.6%.
  - Exact 75.8% -> 89.7% of truth.
  - lncRNA calls were 19.6% -> 20.6% at development defaults; A9 later cut them to ~5%.
  - Bootstrap CIs are about +-0.002-0.006.
- **Lesson:** with 1,000 transcript resamples, CIs reflect transcript sampling only. Across 8 organisms
  a per-organism sign count is the more honest summary.
- **Pointers:** E17; `benchmark_v2/scripts/critiques_old_vs_new.py`, `bootstrap_ci.py`; E16 (literature).

### B10. Release-history benchmark (2026-10-24)
- **Outcome: ADOPTED** (Supplementary Note 1).
- **Numbers:** exact start (complete + random) 0.666 (r2012) -> 0.758 (v5.7.1) -> 0.897 (v7 start rule);
  lncRNA calls 13-19% in every release before the q-value test, 5.0/5.8% after.
- **Lessons:**
  - r2012 reports 3'-partial ORFs 3 nt past the transcript end; clip before scoring (otherwise fragment
    sensitivity 0.37 instead of 0.76). Present this as normalisation.
  - Text guidance (Brian): describe what each release changed in method; do not discuss defects of
    earlier releases.
- **Pointers:** E23; `data/benchmark/history/`; `notebooks/S1_release_history.Rmd`.

## 6. Simulated assemblies

### T1. Simulated GM12878 Trinity assembly (2026-10-25/26)
- **Question:** benchmark inputs are annotated transcripts cut at random. How do methods do on real
  assembler output with known CDSs?
- **Tried:**
  - Polyester reads from GENCODE coding (complete CDS) and lncRNA transcripts at salmon-estimated
    GM12878 expression (52.75M pairs).
  - STAR, then genome-guided Trinity, then minimap2.
  - Each contig labelled by projecting the source CDS (complete / 5' / 3' / internal / noncoding /
    excluded).
- **Outcome: ADOPTED** (Results, Methods, Supp. Fig. S5).
- **Numbers:**
  - 84,972 contigs; 33,939 scored CDSs (83% complete); 20,815 negatives; 32.6% unresolved.
  - TransDecoder stop F1 0.916, exact F1 0.844 (best); Prodigal 0.630 / 0.453, predicting on 52% of
    negatives.
- **Lessons:**
  - Labelling the real GM12878 Trinity assembly against GENCODE leaves 90.9% of contigs unresolved, so
    real assemblies cannot serve as truth.
  - The simulation is more complete than the real assembly (7.0% vs 13.6% of contigs begin inside the CDS).
  - Predictions on excluded contigs count neither way; UTR-only contigs are not counted as false positives.
- **Pointers:** E30; notes/07; `pipeline/simtrinity/`; `data/simtrinity/`; `notebooks/02_benchmark.Rmd`.

### T2. Labeller end tolerance (2026-10-26)
- **Tried:** a strict first pass, then up to 10 nt terminal overhang and in-frame extension through
  <= 10 clipped bases.
- **Outcome: ADOPTED.**
- **Numbers:** unresolved 28,893 -> 27,666; +68 complete, +103 / +189 / +350 partials, +406 negatives.
- **Lesson:** rules added after inspecting a first pass must be reported with both passes' counts.
- **Pointers:** E30.

### T3. Read-depth ablation (2026-10-26/27)
- **Tried:** nested subsamples of the STAR BAM to 1/2 ... 1/32 (seed 20261026); assemble, label, run and
  score each level.
- **Outcome: ADOPTED** (Supp. Figs. S5-S6). The run stopped at d32 by Brian's decision before the
  share-stability rule was met. d64 was killed unscored.
- **Numbers:**
  - Complete share 83% (d1) -> 28% (d32); the E29 truncation share brackets the real assembly (13.6%)
    between d4 and d8.
  - TransDecoder stop F1 0.916-0.929 at all levels; exact F1 0.844 -> 0.806.
  - CPAT exact F1 0.80 -> 0.43.
- **Lessons:**
  - Per-type precision needs a rule for typing predictions on negatives (by the ORF's own ends).
  - With negatives counted, Prodigal passes TransDecoder on 5' partials from d16 and GeneMarkS-T from
    d8; Prodigal and the longest-ORF baseline lead on internals at every level.
  - There is one subsample per level, so differences of about 0.01 are not interpretable.
- **Pointers:** E31; `pipeline/simtrinity/11_*`, `12_*`, `13_*`; `data/simtrinity/depth/`.

### T4. Simulated Arabidopsis Trinity assembly, "E32" (started 2026-10-27/28)
- **Question:** after E31, Brian wanted to consider repeating the ablation on Arabidopsis for comparison
  (memory `transdecoder-paper-project.md`).
- **What exists:**
  - Scripts `pipeline/simtrinity_athal/` (TAIR12 tables; salmon on ENA ERR4140056 Col-0 leaf reads;
    STAR index; Polyester; STAR; Trinity; subsampling; labelling and runs). Their headers cite "E32".
  - Lane logs: d1 finished 2026-10-28 04:06; d2 and d4 started 04:00 (`$P/work/simtrinity_athal/depth/logs/`).
  - E33 used the finished d1 `pytd_v7` run (release exact F1 0.882; truncated share 0.33).
- **Outcome: OPEN.** There is no `evidence/E32*.md` and no row in the evidence index;
  `$T/data/simtrinity_athal/` is empty. Results, status of d2/d4 and any decision are unknown.
- **Pointers:** `pipeline/simtrinity_athal/`; `$P/data/simtrinity_athal/`; `$P/work/simtrinity_athal/`; E33 Data.

### T5. Trinity de novo on GM12878 (2026-10-20)
- **Outcome: ABANDONED:** projected to take weeks on this machine. Partial output is in
  `$P/data/gm12878/trinity_denovo_killed/` and is not used.
- **Pointers:** E11 Caveats.

### T6. Anatomy of missed CDSs and of the 5' UTR-vs-coding decision (E34, 2026-10-29)
- **Question:** where are the truth CDSs with no prediction carrying their stop lost in the PyTD pipeline, what
  would relaxations recover and cost, and why are open 5' ends of true 3' partial / complete CDSs kept open?
- **Tried:** every truth CDS of simulated human d1/d8/d32, athal d1/d16 (ab initio A, homology B) and the GM12878
  Trinity MANE contigs traced through the release runs' own intermediates (selection and refinement re-run offline
  with PyTD's code: 0 differences from the runs); offline q 0.1/0.2/0.5/1.0 and top-candidate fills; fresh runs with
  LongOrfs `-m 30` and Predict `--sensitive` (equal to the offline q 0.2). No code changed.
- **Outcome: OPEN (diagnostic; no method change proposed).**
- **Numbers:**
  - Misses on the simulations are 89-99.6% "candidate fails the coding test" (class b), 0.4-11% "lost to another ORF
    on the contig"; none below the LongOrfs minimum (the labeller drops < 50-codon CDSs). GM12878: 30% of MANE misses
    are fragments < 50 aa, 54% class b, 12% class c, 2% contigs with a premature in-frame stop.
  - Class (b): median 68-83 codons (69-94% < 100 codons vs 5-28% of found CDSs), ~40% of the per-base frame contrast
    of found CDSs of the same length; human ones are mostly short alternative-isoform CDSs (median 45% of the gene's
    longest protein), Arabidopsis ones short main proteins with GC3 outside the training range (51-64% vs 30-33%).
  - q 0.2 recovers 31-47% of misses; on human d1 it costs 1.3 lncRNA-contig calls per recovered CDS (exact F1 -0.6),
    elsewhere F1 rises because negatives are few. Top candidate everywhere: 62-93% recovered, 42-71% of negatives
    called. `-m 30`: 0 recovered on the simulations, 5 of 1,460 on GM12878.
  - 5'-open candidates of 3' partials / complete CDSs are kept open in 10-16% (human) / 15-36% (athal) of trim-needed
    cases: 56-70% of those keeps are UTRs scoring >= 0.4 of typical coding per base (only 15-18% / 0-4% coding in
    another isoform), 13-42% short UTRs (< ~40 nt) vetoed by a negative Kozak score, and 2-5% (human) / < 0.5% (athal)
    of candidates have a true start that is not the first ATG.
- **Lessons:**
  - The remaining sensitivity gap is short, weakly coding CDSs; any acceptance relaxation pays in noncoding calls at
    about 1:1 on a realistic lncRNA load. Which way F1 moves on a simulated assembly depends on its negative count.
  - The keep/trim errors come from the UTR coding-score distribution overlapping coding (model / short-segment
    noise) and from the Kozak term on short UTRs, not from the first-ATG rule.
- **Pointers:** E34; `data/missed_cds/`; `pipeline/missed_cds/`; `$P/data/missed_cds/`.

## 7. GM12878 application

### G1. Three assemblies vs MANE/GENCODE (2026-10-20/21; release 2026-10-24)
- **Outcome: ADOPTED.**
- **Numbers:** exact CDS given an assembled transcript 0.935 / 0.927 / 0.906 (LRAA / StringTie /
  Trinity; MANE, development). Release: 0.929 / 0.915 / 0.893, with unannotated CDSs 341 / 1,687 /
  2,528 (development 885 / 4,127 / 11,235).
- **Lessons:**
  - "MANE precision" is low by construction; GENCODE exact/compatible shares are the fair precision
    measures.
  - 4.8% of Trinity ORFs were not projected to the genome (soft-clipped or split alignments).
- **Pointers:** E03, E11, E25; notes/06; `genome_annotation/scripts/eval_genome_cds.py`.

### G2. Ribo-seq support (Cenik 2015, GSE65912) (2026-10-21)
- **Outcome: ADOPTED.**
- **Numbers:** exact matches 83-89% covered (>= 10 P-sites); unannotated 0.6-4.3%; at TPM 5-20, 89% vs
  3.4%, so low abundance does not explain it.
- **Lessons:**
  - RNase A+T1 digestion gives weak periodicity (frame-0 ~0.55), so the frame test is conservative.
  - Rep1 contributes little (1.7M of 11.5M footprints).
  - Unique alignments only, so repeat-derived ORFs are invisible.
- **Pointers:** E01; `analyses/riboseq/`.

### G3. MS peptide support (Sinitcyn 2023) (2026-10-21)
- **Outcome: ADOPTED.**
- **Numbers:** exact >= 100 aa 92-94% with a peptide (release 94-95%). Unannotated with a peptide: 9-24%
  (development; first written as 10-25%) and 30-49% at release (smaller denominator). One credible
  non-reference case: C3orf86P (now predicted only with `--sensitive`).
- **Lessons:**
  - The published search did not contain our ORFs, so this confirms known proteins and cannot discover
    new ones.
  - Peptides on "novel" ORFs mostly come from paralog/pseudogene copies or variant alleles.
- **Pointers:** E02, E25; `analyses/proteomics/`.

### G4. Expressed-gene sensitivity (2026-10-21)
- **Outcome: ADOPTED.**
- **Numbers:** StringTie exact MANE CDS for 74% (gene TPM >= 1) to 81% (TPM >= 10) of expressed MANE
  genes; exact given assembled stays 0.91-0.94.
- **Lesson:** gene-level TPM suits MANE (one isoform per gene) but not all GENCODE isoforms. The
  expression command was not scripted at run time.
- **Pointers:** E05.

### G5. Protein-match scoring on GM12878 (2026-10-23 onward)
- **Tried:** score predicted proteins as identical to, or sharing an exact 30-aa window with, a GENCODE
  v50 protein; "non-matching" otherwise. Join to Ribo-seq/MS support of dropped ORFs.
- **Outcome: ADOPTED** for comparing rules (E18-E22, E25, E26).
- **Numbers:** release default LRAA 32,751 matching / 277 non-matching; StringTie 35,294 / 1,185 (fixed
  short-ORF test 33,587 / 1,005 and 36,547 / 4,119).
- **Pointers:** `benchmark_v2/scripts/score_protein_matches.py`; E25.

### G6. Release reruns on GM12878 (2026-10-24)
- **Outcome: ADOPTED** (all GM12878 numbers in the paper use the release runs).
- **Numbers:** catalogued Ribo-seq ORF stop matches 101 / 95 / 186 (release) vs 178 / 217 / 455 (fixed test);
  precise 79 / 81, sensitive 144 / 165.
- **Lesson:** stricter acceptance trades unannotated noise for real translated small ORFs. Report the cost.
- **Pointers:** E25; `pipeline/derive/06_*`, `07_*`, `08_*`.

## 8. Infrastructure and process lessons
Each item cites its source. Items marked "undocumented" were listed by the coordinator, but no record
or note describes the incident; details are unknown.

**Compute and hosts**
- **No UGER.** The first ~500 benchmark and start-grid jobs ran on UGER before 2026-10-17; UGER is
  billed, so all runs since are local (`run_all.sh`, `xargs -P`). Low-GC threshold chunks ran on UGER;
  `run_low_gc_chunks.sh` now runs them locally. (notes/04; E10 Gaps; memory `no-uger-run-locally.md`)
- **Shared host.** methods01 has 28 cores; `nproc` prints 8 because OMP_NUM_THREADS=8; load is ~40
  from other users. Runtimes are wall clock under load: use ratios, not absolute times. (memory
  `no-uger-run-locally.md`; E07, E08 Caveats)
- **Use methods02 for memory.** The simulated-library STAR run (peak 37 GB) ran on methods02, because
  methods01 had 26-31 GB free. methods02 sees `/broad/hptmp`; launch with `nohup setsid`. (E30;
  memory `methods02-for-memory.md`)
- **STAR BAM sort RAM.** `--limitBAMsortRAM 8000000000` failed at sorting, because simulated reads are
  written grouped by source transcript and the sort needed ~30 GB; use 40 GB. A memory-gated launcher
  (`run_star_when_memory_v2.sh`, ulimit -v) waited for >= 42 GB. (E30; `$P/work/simtrinity/`)
- **Thread caps for sklearn.** HistGradientBoosting used all OpenMP threads by default (43 min CPU in
  6.6 min wall); later runs set OMP_NUM_THREADS=3. (E33)

**Running scripts and job control**
- **Don't edit running bash scripts on hptmp.** Bash reads scripts incrementally. An in-place edit of
  `short_orf_rule_run.sh` broke running jobs (2026-10-24), and even a write-then-`mv` replacement of
  `run_tool.sh` made multi-day CodAn jobs die with "unexpected EOF" before writing `runtime_sec.txt`.
  Write a new versioned script (`*_v2.sh`) instead. (memory `dont-edit-running-scripts.md`; versioned
  launchers such as `run_star_when_memory_v2.sh`, `01_features_all_v1.sh`)
- **Don't move a shared code tree under running jobs.** launch17 ran from `work/trees/PY_dev`, which
  moved to 44e417a mid-run, so 5 T. thermophila Predict steps had to be re-run from a worktree pinned at
  48c870c. Release runs use pinned worktrees (`PY_v7`, `TD_v7`). (E07 Caveats, E24)
- **Run labels track trees, not commits; runs can predate their commits.** The development final-default
  runs finished minutes before 48c870c/de3f57e were committed; only 2 inputs were re-verified.
  `pytd_version.txt` said be23241 while the tree had c78845d's content. Record the commit in each run
  directory and commit before launching. (E07, E09, E10 Gaps, E24)
- **Fresh runs, never resume, for rule comparisons** (Brian). Incomplete directories are moved aside,
  not deleted. (E19, E20; `short_orf_rule_run.sh`)
- **`pgrep -f` in waiter loops (self-matching).** Observed 2026-10-27 (main session cleanup). Two
  background waiters had loops of the form `until grep -q ALL_DONE $L || ! pgrep -f bench_v3/launch.sh;
  do sleep 600; done` (and one on `b6lq997nx_never`). The pattern also matched the waiter's own
  `bash -c` command line, so `pgrep -f` always succeeded and the loops ran for 2.5–9 days after their
  launchers had ended. They were stopped by task/PID. The first simulated-Trinity waiter
  (`$P/work/simtrinity/run_trinity_after_star.sh`, `pgrep -f run_star_when_memory_v2.sh`) had the same
  hazard; `run_trinity_after_star_v2.sh` waits on log lines only. Rule: wait on a PID (`kill -0 PID`) or
  on a marker line in a log, never on `pgrep -f` of a string that appears in your own command.
- **Duplicate waiter shells.** Observed 2026-10-27 during E31: the depth-ablation agent started a new
  background `until grep -q LEVEL_DONE ... laneE.out; do sleep 60; done` waiter on each check (~every
  10 min, when its foreground tool call timed out), leaving ~30 identical loops over ~5 h. Harmless
  for results (they only read a log) but they cluttered the process table; all but the newest were
  killed by PID. Rule: keep one waiter per condition; reuse it or check the log directly.
- **Kill by exact PID.** The d64 Trinity run and two stalled `flock` processes were killed by PID. (E31)

**Filesystem and tools**
- **NFS `flock` stalls in Trinity ParaFly.** At d32 one partition stalled ~1.5 h at the end of Butterfly:
  the ParaFly xargs drop-in was waiting on an NFS `flock` on `butterfly_commands.completed`. Killing the
  two stalled `flock` processes let it finish; 18,352 / 18,352 partitions completed. (E31)
- **Trinity conda environment.**
  - Bioconda ParaFly ran commands serially whatever `-CPU` said; it was replaced by a bash/xargs
    drop-in (`ParaFly_xargs.sh`; bioconda-recipes issue #70548).
  - salmon must be pinned to 1.10.3 (2.x breaks `salmon_runner.pl`); one env was broken by a salmon
    install that downgraded Trinity.
  - Unset PERL5LIB for Trinity.
  - Next time, use the Trinity container via apptainer.
  - Sources: E11; memory `use-containers-for-complex-tools.md`.
- **Polyester patches.**
  - Fragment lengths < 50 nt are redrawn (normal draws <= 0 crash `subseq`).
  - Position 101 reuses the position-100 `illumina5` profile.
  - With `strand_specific = TRUE`, R1 = Polyester mate 2 for an RF library; strand was verified at
    transcript and genome level.
  - Source: E30.
- **CodAn input** must be line-wrapped (80 nt); its env needed a `libnsl.so.1` symlink. (E07; notes/04;
  memory `transdecoder-paper-project.md`)
- **GeneMarkS-T** needs a licence key in the home directory; neither key nor software goes in the repo. (E07)
- **External downloads can disappear:** LNCipedia returned 404 on 2026-10-23. (E18)

**Numbers, notebooks and figures**
- **Double rounding.** Medians computed from 4-decimal rates differed by 0.001 from medians of counts;
  recompute rates from counts before taking medians. (E23 Correction; `notebooks/S1_release_history.Rmd`)
- **Rounded model files.** The Kozak model written to 3 decimals made 0-4 predictions per input differ
  between nominally equivalent weights. A Kozak-score rounding near-tie was the last Perl/PyTD
  difference (eaea367). (E13 addendum; E09 Checks)
- **Knitting re-renders figures.** `scripts/knit_notebooks.sh` writes figure PDFs/PNGs as well as HTML,
  so any knit changes figure files in git (commit 90465f8 "Re-knit figure PDFs" changed
  `fig2_benchmark.pdf` and `figS4_rankings.pdf` bytes). Knit only the notebook you changed, and check
  `check_report()` (0 mismatches) afterwards. (notebooks/README.md; E28, E31 Checks)
- **Rescoring scope.** `score_all.py` walks every run directory, so superseded or separately scored
  inputs (`fragments_mixed_t*`, `fragments_internal`) leak into the 84-input summary unless filtered or
  moved to `$P/_trash`. (E28 Housekeeping)
- **Records vs notes.** Where early notes disagree with committed files, the files and E-records take
  precedence (notes/05, notes/06 "Corrections"; E07 Gaps, E10 Caveats).
- **Large outputs stay out of the repo** (`$P/work`, `$P/data`); E33's 177 MB of intermediates were moved
  out after analysis. (evidence/README.md; E33 Storage note)

## Inconsistencies noticed while building this log (2026-10-28)
1. **Index statuses are stale.**
   - evidence/README.md lists E18, E19, E20 and E21 as "planned".
   - E20's rule and E21's presets were adopted in the release (E24, E25); E18/E19 were superseded.
   - E07 is listed as "complete except CodAn scores (71/80 runs)" and E15 as "Fig. 2 ... to redraw
     after CodAn", but `results_benchmark.tex` Table 1 has CodAn values and `discussion.tex` has a
     todo "CodAn is now in Table".
   - E07 itself does not record when all 80 main-benchmark CodAn runs finished.
   **Resolved 2026-10-28:** the E07, E15, E18-E21 index rows and records now carry dated status updates. CodAn is 79/80 in the main benchmark; the missing *T. thermophila* complete.unstranded.wRand run stopped without output and is not applicable. Fig. 2 includes CodAn.
2. **E28 contradicts itself.** Its Status says all CodAn runs are done and scored (2026-10-27), but its
   Caveats still say "CodAn results are pending".
   **Resolved 2026-10-28:** the E28 caveat is annotated as resolved.
3. **E25 says the loss of catalogued Ribo-seq ORFs** (101/95/186 vs 178/217/455) "must be stated in
   Results/Discussion". No such statement was found in `manuscript/sections/*.tex`: "catalog" appears
   only for lncRNA predictions in `results_short_orfs.tex` and `discussion.tex`.
   *Open (manuscript): investigated 2026-10-28; proposed fix sent to the coordinator.*
4. **Calibration on lncRNAs.** E20 calls the lncRNA calibration "roughly nominal" (development runs).
   `results_short_orfs.tex` reports observed FDR above nominal at release (0.080 / 0.102 at q 0.05) and
   carries a todo about the wording.
   *Open (manuscript): investigated 2026-10-28. 19-29% of release lncRNA calls are long enough to be accepted on length alone (>= L\* = 726 nt), so 1 - precision is not a q-value calibration.* **Update 2026-10-28:** the Results now say this; the tally is reproducible with `pipeline/derive/04_lnc_calls_by_autoaccept_length.py` -> `final_v7/lnc_calls_by_autoaccept_length.tsv` (E25).
5. **Stale caveats about outdated Methods text.** E08, E09 and E10 Caveats say the Methods still describe
   the PWM procedure, "differ only in start codons" and GC = 2q(C). None of these phrases is in the
   current `methods_algorithm.tex`; the caveats appear out of date.
   **Resolved 2026-10-28:** E08, E09 and E10 (and E07's methods_benchmark caveat) are annotated as out of date, with status updates.
6. **E12 Release update** says Fig. 5d,e still use the development homology runs; E15 says panels d,e
   come from the release homology runs (2026-10-26). E12 appears out of date.
   **Resolved 2026-10-28:** E12 was wrong. `org_meta.tsv` gives `homology_source = pytd_v7_homology` for all 8 organisms, and Fig. 5d,e use the release comparison. E12 now has a status update.
7. **Mixed-input trial agreement.** `results_benchmark.tex` says no median exact F1 differed by more than
   0.025 between trials. E28 Results gives at most 0.008 (before CodAn), and E28 Status gives at most
   0.016 for CodAn. The text is a looser bound, not a contradiction, but the source of "0.025" is unclear.
   *Investigated 2026-10-28: 0.025 is the largest per-type `f1_exact_trial_range` in `final_v7/mixed_medians.tsv` (0.0248, TransDecoder v6.0.0 internal fragments). Over whole inputs the maximum is 0.0084. Manuscript wording open.*
8. **E32 has no record.** The Arabidopsis simulation is called "E32" in script headers and in E33, but
   there is no `evidence/E32*.md`, no index row, and `data/simtrinity_athal/` is empty.
9. **Regression at release commits.** `results_algorithm.tex` has a todo: the sample-set comparison is
   from PyTD 44e417a / Perl ab38ca6, not rerun at the release commits (E08).
   *Open: no release-commit regression run exists (`data/regression/` ends at kozak_bits_44e417a_ab38ca6).*
10. **notes/06 values superseded by records.** notes/06 still says 10-25% of unannotated predictions share
    a peptide; E02 corrects this to 9-24%. Its "26-38%" figure was corrected in notes/06 itself.
   **Resolved 2026-10-28:** notes/06 now carries a dated correction citing E02.
11. **First-ATG baseline.** E07 says the `longest_orf_atg` result (exact F1 0.942 vs 0.910) is "not yet in
    the manuscript". It is still absent (no "0.942" in the sections). Whether to report it is undecided.
   **Resolved 2026-10-28:** Brian approved one sentence. The baseline was recomputed at the release (E07 Addendum; `first_atg_baseline.tsv`): complete exact F1 0.956 vs 0.912, fragments 0.484 vs 0.822.

## Outcome counts
Entries: 58 (A1-A13, S1-S14, O1-O4, H1-H5, B1-B10, T1-T6, G1-G6).

| Outcome | Count | Entries |
|---|---|---|
| ADOPTED | 33 | A1-A4, A9-A12, S4, S6-S9, S11, O1, H1, H2, B1-B4, B8, B10, T1-T4, G1-G6 |
| NOT ADOPTED | 14 | A5, A7, A8, A13, S2, S3, S5, S10, S13, O3, O4, H3, H4, H5 |
| PARTIAL | 4 | S1, O2, B5, B7 |
| OPEN | 5 | A6, S12, S14, B9, T6 |
| ABANDONED | 2 | B6, T5 |

Rejected variants inside ADOPTED entries: the decoy-derived long-ORF length (A10), the multiplicative
Kozak form (S6), the self-trained keep/trim classifiers (S11, steps 1-2), and top-n / looser-e-value / UniRef50 evidence for the start assist (S11 steps 4-5; also H4, H5). Open sub-questions inside other entries: which lncRNA ORFs gain homology hits (H1);
stating the catalogued Ribo-seq cost in the manuscript (G6).
