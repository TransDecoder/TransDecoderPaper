# E12. True CDSs that fail the coding criteria

## Claims supported
Cited by Results, "True coding regions that fail the coding test, and two fallbacks"
(results_coding_failures.tex, Fig. 5) and Methods, "GC-calibrated long-ORF rule". It supports:

- **Methods, "GC-calibrated long-ORF rule"** (`methods_algorithm.tex`). The rule "is intended to recover genes with atypical composition". This record is the analysis that E10 Gaps lists as missing.
- **Results/Discussion, if added:**
  1. 4.9% of annotated CDSs (bench_v3 complete, 8 organisms) fail rule 3. In order, the failure modes are:
     - the short-ORF guard;
     - S1 ≤ 0;
     - a shifted frame scoring at least S1.

     Failures fall with length: 67% of CDSs under 300 nt fail, 0.3–0.6% of CDSs of 1 kb or more.
  2. The long-ORF rule rescues 447 of the 6,627 failing CDSs (397 end up predicted). It does not apply to *T. thermophila*.
  3. The failures have three main causes:
     - a GC3/codon-usage minority class (fly chromosome 4/Y, low-GC3 *Chlamydomonas*, *Tetrahymena* ribosomal proteins). A GC3-matched model rescues most of these.
     - low-complexity or repetitive proteins (KRTAP, LCE, FLG/HRNR, mucins, ZNF arrays). These are 78% of failures in human and mouse.
     - unconserved, possibly mis-annotated models, mostly *Tetrahymena* "hypothetical protein".
  4. In GM12878, most assembled MANE CDSs that are not predicted exactly have a different start. 17% (StringTie) and 30% (LRAA) of them fail the coding test, mostly through the short-ORF guard.
  5. Predictions kept only by the long-ORF rule are 0–1.4% of predictions. Among them are false calls on lncRNA and random negatives.

## Data
- **Benchmark inputs:**
  - `$P/data/bench_v3/inputs/<org>/{complete,complete_fullutr,complete.wLnc,complete.wRand}.fa` and their `.truth.tsv` files;
  - 8 organisms (record E07).
- **Annotation, for gene names, descriptions and chromosomes:** `$P/data/annot/<org>/...`, listed in `$T/benchmark_v2/truth_sources.tsv`.
- **Final-default runs (inputs to this analysis):** `$P/work/bench_v3/runs/<org>/<input>/pytd_dev_v2/`. Files used:
  - `longest_orfs.{gff3,cds,cds.scores}`;
  - `hexamer.scores`;
  - `base_freqs.dat`;
  - `__checkpoints_TDpredict/top_training_orfs.cds`;
  - `predict.log`;
  - `preds.tsv`;
  - `input.fa.transdecoder.gff3`.
- **Homology annotation of CDSs:** `$P/work/bench_v3/runs/<org>/<input>/pytd_homology/{pfam.domtblout,blastp.outfmt6,input.fa.transdecoder_dir/longest_orfs.gff3}`.
  - These are old-code runs: LongOrfs −m 100, DIAMOND against Swiss-Prot minus the organism's taxid, and hmmsearch with `--cut_ga`.
  - CDSs are matched to the homology-run ORFs by transcript and stop coordinate.
- **GM12878 runs:** `$P/work/gm12878_td/{stringtie,lraa,trinity}/final_abinitio/`.
- **GM12878 evaluation:**
  - `$P/work/gm12878_td/eval/status/<asm>_final_abinitio.{ref_status,pred_status}.tsv`;
  - assembly structures `$P/data/gm12878/star_stringtie/stringtie.gtf` and `$P/data/gm12878/longread_lraa/lraa_denovo/LRAA.gtf`;
  - MANE `$P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`.

## Software
- **Python:** 3.11 (`$P/env/td`), with numpy 2.4.6, scipy 1.17.1 and pandas 3.0.6.
- **PyTD:** `$P/work/trees/PY_dev` at 48c870c. The analysis imports it read-only: `Translator`, `genome.read_structures`, and `eval_genome_cds.py` from `$T/genome_annotation/scripts`.
- **This repo:** at 5f3c269 when the analysis was run. The scripts were committed later
  (analyses/long_orf_failures/, last change d90c38a):
  - `analyses/long_orf_failures/long_orf_failures.py`;
  - `analyses/long_orf_failures/markov.py`, a vectorised copy of PyTD's Markov training and scoring;
  - `analyses/long_orf_failures/features.py`;
  - `analyses/long_orf_failures/summarize.py`.

## Commands
```bash
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; T=$P/TransDecoderPaper
cd $T/analyses/long_orf_failures
$P/env/td/bin/python long_orf_failures.py all --jobs 4
# = bench (complete), bench --input complete_fullutr, gm12878, altrules, summarize
```

The run in this record was done step by step, on 2026-10-21, on the shared host with at most 4 processes:
1. `bench`;
2. `bench --input complete_fullutr`;
3. `gm12878`;
4. `altrules`;
5. `summarize`.

## Outputs
All outputs are in `$T/analyses/long_orf_failures/results/`.

- **Per-CDS tables:** `per_cds/<org>.per_cds.tsv.gz` and `per_cds/<org>.complete_fullutr.per_cds.tsv.gz`, one row per truth CDS. Each row has:
  - the matched ORF and its scores S1–S6, the contrast, and the rule flags;
  - the category and outcome;
  - GC, GC3 and codon log-likelihood ratio vs the training set;
  - amino-acid composition and low-complexity measures;
  - Pfam and Swiss-Prot hits;
  - pass/fail under alternative Markov models.

  Per-organism parameters are in `per_cds/<org>.meta.tsv`.
- **Summary tables:**
  - `summary_by_org.tsv`, `summary_by_length.tsv`;
  - `oracle_models.tsv` (alternative training sets);
  - `features_fail_vs_pass.tsv` (failing vs length-matched passing CDSs);
  - `fail_rate_by_feature.tsv` (quintiles);
  - `family_keyword_enrichment.tsv`, `fail_by_chromosome.tsv`, `ribosomal_proteins.tsv`;
  - `explanation_flags.tsv` and `failing_cds_attribution.tsv.gz`;
  - `utr_effect.tsv` (trimmed vs full UTRs);
  - `examples_failing.tsv`, `org_meta.tsv`.
- **Long-ORF rule and alternative rules:**
  - `length_rule_only_predictions.tsv`: predictions kept by the length rule alone;
  - `alt_rules_negatives.tsv`: exploratory, rules without S1 > 0, on random and lncRNA negatives.
- **GM12878:**
  - `gm12878_mane_assembled.tsv.gz`;
  - `gm12878_summary.tsv`;
  - `gm12878_fail_examples.tsv`;
  - `gm12878_fail_vs_pass_features.tsv`.
- **Write-up:** `FINDINGS.md`.

## Results
**Rule 3 failures on bench_v3 complete** (`summary_by_org.tsv`):

| org | matched CDSs | fail rule 3 | short guard | S1 ≤ 0 | shifted ≥ S1 | rescued by L\* (predicted) | lost | fail (ORF ≥ 300 nt) | L\* |
|---|---|---|---|---|---|---|---|---|---|
| athal | 26,772 | 1,260 (4.71%) | 820 | 232 | 208 | 60 (54) | 1,200 | 0.92% | 663 |
| creinhardtii | 17,621 | 1,025 (5.82%) | 322 | 652 | 51 | 8 (7) | 1,017 | 2.09% | 1,347 |
| dmel | 13,784 | 797 (5.78%) | 311 | 438 | 48 | 183 (181) | 614 | 2.80% | 804 |
| human | 19,144 | 406 (2.12%) | 195 | 35 | 176 | 42 (41) | 364 | 0.86% | 804 |
| mouse | 21,361 | 637 (2.98%) | 224 | 54 | 359 | 150 (110) | 487 | 1.71% | 804 |
| pfalciparum | 5,277 | 155 (2.94%) | 144 | 7 | 4 | 1 (1) | 154 | 0.12% | 519 |
| spombe | 5,081 | 165 (3.25%) | 97 | 62 | 6 | 3 (3) | 162 | 0.92% | 615 |
| tthermophila | 26,389 | 2,182 (8.27%) | 234 | 1,611 | 337 | 0 | 2,182 | 6.18% | off |
| all | 135,429 | 6,627 (4.89%) | 2,347 | 3,091 | 1,189 | 447 (397) | 6,180 | 2.38% | |

**By CDS length, pooled** (fail % / still lost %):

| CDS length (nt) | Fail | Still lost |
|---|---|---|
| < 300 | 67.0% | 67.0% |
| 300–599 | 9.35% | 9.33% |
| 600–999 | 2.32% | 1.87% |
| 1,000–1,999 | 0.60% | 0.16% |
| 2,000–4,999 | 0.30% | 0.02% |
| ≥ 5 kb | 0.41% | 0.02% |

The lost CDSs of 1 kb or more are all in *T. thermophila* or *C. reinhardtii*.

**Training set and alternative models** (ORFs ≥ 300 nt, pooled; `oracle_models.tsv`):
- TD training ORFs that are annotated CDSs: 499–500 of 500 in every organism.
- Failure rates:

  | Model | Failure rate |
  |---|---|
  | TD's model | 2.38% |
  | all annotated CDSs (2-fold) | 2.08% |
  | 500 random CDSs | 2.49% |
  | 500 longest CDSs | 2.47% |
  | 500 GC3-quintile-matched CDSs | 1.95% |

- Score failures rescued by the GC3-matched model: 25% overall; 69% dmel, 52% mouse, 36% human, 35% athal, 20% creinhardtii, 8% tthermophila, 7% spombe.

**Explanations** (ORFs ≥ 300 nt; failing vs length-matched passing; `explanation_flags.tsv`, `features_fail_vs_pass.tsv`, `fail_by_chromosome.tsv`, `family_keyword_enrichment.tsv`, `ribosomal_proteins.tsv`):
- **dmel:** failing median GC3 is 0.385 vs 0.655 for controls. GC3 quintile 1 fails at 13.5%, the other quintiles at 0.04–0.3%. Chromosome 4: 40/80 fail (39 rescued by the GC3 model). Chromosome Y: 8/22 fail.
- **human/mouse:** 78% of failures are low-complexity or repetitive, vs 13% of controls. Shifted frame ≥ S1 accounts for 134/162 (human) and 318/358 (mouse). Families:
  - KRTAP 52/72 (human) and 78/92 (mouse);
  - LCE 12/14 and 20/20;
  - human chromosome Y 10/42.
- **Human failing vs passing composition:** Cys 14.3% vs 2.5%; SEG-like low-complexity fraction 0.42 vs 0.11.
- **tthermophila:** 97% of failures have no Swiss-Prot or Pfam hit, vs 58% of controls. 26% have a stronger out-of-frame ORF covering ≥ 50% of the CDS, vs 0.6%. Ribosomal proteins: 21/60 fail S1 > 0 (median GC3 0.47 vs 0.23 genome-wide), and 17 of those pass with the GC3-matched model.
- **Exclusive attribution, all organisms (n = 3,114):**
  - low complexity: 32%;
  - GC3: 33%;
  - stronger out-of-frame ORF: 6%;
  - no homology: 26%;
  - unexplained: 2.4%.

**UTRs** (`utr_effect.tsv`): rule 3 failure is unchanged with full UTRs (all organisms within about 1 point). The number of passing CDSs outcompeted by another ORF rises, for example dmel 12→157 and tthermophila 79→199.

**Predictions kept by the length rule alone** (`length_rule_only_predictions.tsv`):
- Complete inputs: 0–1.37% of predictions, 398 of 402 being the truth CDS.
- Negatives, true vs false calls:

  | Input | True | False |
  |---|---|---|
  | human lncRNA | 47 | 99 |
  | human random | 40 | 82 |
  | mouse lncRNA | 103 | 48 |
  | creinhardtii random | 8 | 58 |

- GM12878 StringTie: 148 such calls, 116 distinct CDSs. Of these, 11 match MANE exactly, 35 are compatible with GENCODE, and 46 do not overlap GENCODE.

**Alternative rules, exploratory** (`alt_rules_negatives.tsv`): dropping S1 > 0 adds +6 to +537 true CDSs as candidates per input (0.06–2.0%). Random transcripts carrying a candidate rise from 154 to 3,548 (human) and from 467 to 19,684 (tthermophila).

**GM12878 MANE, assembled but not exact** (`gm12878_summary.tsv`):

| | StringTie | LRAA |
|---|---|---|
| Assembled but not exact | 651 | 412 |
| Right stop, different start | 441 | 256 |
| Passes, outcompeted | 75 | 18 |
| Fails rule 3 | 113 (17.4%) | 125 (30.3%) |
| – short guard | 80 | 88 |
| – S1 ≤ 0 | 17 | 19 |
| – shifted frame | 16 | 18 |
| No ORF | 22 | 13 |
| All assembled MANE failing rule 3 | 1.43% | 2.07% |
| Failures that also fail in bench human | 90/113 | 81/125 |

## Checks
- **`markov.py` against PyTD.** Training on `top_training_orfs.cds` reproduces `hexamer.scores`, max abs diff 9e-16. 5,000 ORFs per organism rescored against `longest_orfs.cds.scores` gave 0 mismatches at 0.01 (human, tthermophila). μ was recomputed and matches `predict.log` (assert to 1e-3, all organisms).
- **Rule application against TD's own output.** I compared my classification with TD's own output across all 8 organisms:
  - Every "pass_coding" CDS is either predicted (128,660) or outcompeted (142).
  - 397 of the 447 rescued CDSs are predicted.
  - 6,179 of the 6,180 "fail_coding_lost" CDSs have no prediction with their stop. The exception is *T. thermophila* g29832.t1: the ±2 nt stop tolerance matched an overlapping ORF in a different frame (stop 433 vs 431).
  - The length-rule-only prediction count in human complete (45) is consistent with the 41 rescued-and-predicted truth CDSs plus 4 non-truth ORFs.
- **GM12878 counts.** The MANE assembled and exact counts equal `MANE_assembled`/`MANE_exact` in the eval summaries (StringTie 8,919/8,268; LRAA 6,294/5,882).
- **Truth-set organelle check.** Only 7 *S. pombe* mitochondrial CDSs are in the truth sets, and none fails the score test.

## Caveats
- **Reviewed and cited** (status updated 2026-10-22; earlier drafts of this record said uncommitted).
- **The truth sets are annotations.** Unconserved *Tetrahymena*, *S. pombe* and *Chlamydomonas* models, and mouse Gm genes, may not be real proteins. Their failure may be correct behaviour.
- **Pfam and Swiss-Prot come from the old-code `pytd_homology` LongOrfs (−m 100).** CDSs under 100 aa have no homology data. "No Swiss-Prot hit" means no curated homolog outside the organism; it does not mean species-specific.
- **Low-complexity measures are proxies:**
  - a SEG-like 12-aa window under 2.2 bits;
  - top-3 residue fraction;
  - 3-mer recurrence.

  These are not SEG, a disorder predictor or repeat annotation.
- **Alternative-model training sets are split two-fold by gene, not by family.** Paralogs can leak across the split, which may flatter the truth-trained models.
- **Alternative-rule and length-only negative counts are candidate counts** before single-best selection (upper bounds). They are exploratory, and the rules were not tuned.
- **Input runs predate their commits.** The `pytd_dev_v2` runs predate the 48c870c/de3f57e commits (E10 Gaps). Equivalence is verified only for spombe and creinhardtii complete.
- **Interpretations not tested here:**
  - HGT: no gene set tested;
  - reduced codon-usage selection on fly chromosome 4/Y;
  - frame errors in *Tetrahymena* gene models.

## Addendum: length fallback vs homology fallback (2026-10-21)
```bash
$P/env/td/bin/python $T/analyses/long_orf_failures/length_vs_homology_fallback.py   # -> results/length_vs_homology_fallback.tsv
```
Of the 447 true CDSs that fail the coding test and are retained by the long-ORF rule:
- 362 (81%) also have a homology hit (Swiss-Prot minus the organism's own entries, or Pfam; from the `pytd_homology` runs).
- **85 (19%) have no homology hit**, so only the length rule recovers them:

| Organism | Length-only |
|---|---|
| fly | 40 |
| mouse | 21 |
| *A. thaliana* | 14 |
| *C. reinhardtii* | 5 |
| human | 2 |
| *S. pombe* | 2 |
| *P. falciparum* | 1 |

Swiss-Prot keeps close relatives of each organism, so for a truly novel organism this is a lower bound. Manuscript use: the Results subsection "True coding regions that fail the coding test motivate a length-based fallback".
