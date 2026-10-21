# GM12878: TransDecoder on three assemblies, final defaults (2026-10-21)

**Setup**
- PyTD at PY_dev 48c870c, ab initio: LongOrfs -m 50 + short-ORF guard 200, additive Kozak start rule.
- Each run: Predict, then genome-map, then `eval_genome_cds.py` (MANE v1.5, GENCODE v50), then Ribo-seq stop matching.
- Scripts: `work/gm12878_td/final_defaults.sh` and `overnight.sh`.
- Tables: `work/gm12878_td/eval/{lraa,stringtie,trinity}_final_abinitio.summary.tsv`.

**Assemblies**
- LRAA: long reads, minimap2 -> LRAA.
- StringTie: STAR -> StringTie.
- Trinity: genome-guided on the same STAR BAM (ENCSR000AED rep1, RF).
  - Trinity 2.15.2 + salmon 1.10.3, with the parallel ParaFly replacement.
  - 77,695 partitions, none failed; 428,960 transcripts (319,068 Trinity "genes").
  - Transcripts mapped back with minimap2 splice:hq; 99.8% mapped.

| | LRAA | StringTie | Trinity GG |
|---|---|---|---|
| transcripts (mapped) | 39,530 | 83,342 | 428,223 |
| distinct predicted CDS | 23,922 | 34,404 | 72,323 |
| complete / partial | 18,807 / 5,115 | 27,034 / 7,370 | 52,123 / 20,200 |
| MANE exact (sens) | 5,882 (0.306) | 8,268 (0.430) | 6,951 (0.362) |
| MANE compatible (sens) | 0.395 | 0.569 | 0.547 |
| MANE exact, given assembled | 0.935 | 0.927 | 0.906 |
| GENCODE complete exact (sens) | 0.192 | 0.234 | 0.209 |
| GENCODE complete compatible (sens) | 0.269 | 0.352 | 0.446 |
| GENCODE exact, given assembled | 0.889 | 0.886 | 0.876 |
| predicted complete: GENCODE-exact frac | 0.639 | 0.511 | 0.265 |
| all predictions: no GENCODE CDS overlap | 3.7% | 12.0% | 15.5% |
| Ribo-seq stop matches (exact) | 178 (33) | 217 (61) | 455 (130) |

**Reading**
- Given a correctly assembled transcript, TransDecoder recovers the exact annotated CDS 91-93% of the time (MANE) and 88-89% of the time (GENCODE), whatever the assembler. The ceiling is set by assembly, not by CDS prediction.
- Trinity GG has the highest GENCODE compatible sensitivity (0.45). But many of its transcripts are fragments, with 20k partial ORFs, so its exact-match fraction is lower.
- 4,624 of 97,088 Trinity ORFs (4.8%) were not projected to the genome, because the CDS lay partly outside the aligned region (soft-clipped or split alignments).

The homology-retention run on Trinity (DIAMOND + Pfam) is still running.

## Minimum ORF length vs sensitivity / precision (post-hoc filter, 2026-10-21)
- Method: `eval_genome_cds.py --min-aa L` on the final-default predictions. Script: `genome_annotation/scripts/minaa_sweep.sh`. Table: `genome_annotation/results/gm12878_minaa_sweep.tsv`.

**50 vs 100 aa barely changes precision.**

| | GENCODE-exact share of predictions, 50 -> 100 | no-GENCODE-overlap, 50 -> 100 | MANE exact given assembled, 50 -> 100 |
|---|---|---|---|
| LRAA | 0.504 -> 0.506 | 3.6% -> 3.2% | 0.934 -> 0.916 |
| StringTie | 0.404 -> 0.427 | 11.6% -> 10.3% | 0.927 -> 0.912 |
| Trinity | 0.194 -> 0.207 | 14.5% -> 12.3% | 0.905 -> 0.889 |

**The step is between 100 and 150 aa.** Predictions with no GENCODE overlap drop to 1.1% / 4.3% / 3.3%, but MANE sensitivity given assembled falls to about 0.85. Per-bin view for StringTie and Trinity:
- No GENCODE overlap: 26-38% of the 100-199 aa predictions, about as noisy as 50-99 aa; 1-5% at 200 aa or more.
- These are chance ORFs on noncoding transcripts, mostly in the short-read assemblies.

**What limits precision is not the minimum length:**
1. Partial ORFs from fragmented transcripts can never match exactly (Trinity: 20k partials).
2. The MANE precision row counts correct alternative-isoform CDSs as misses, so it is low by construction. The GENCODE exact/compatible shares are the fair measure.
3. 100-199 aa ORFs on noncoding transcripts.

## Sensitivity restricted to expressed genes (2026-10-21)
- Expression: `stringtie -e --rf -G gencode.v50` on the same STAR BAM (`data/gm12878/expr/`). Gene TPM is the sum of its transcripts' TPMs.
- Script: `genome_annotation/scripts/expressed_sensitivity.py`. Table: `genome_annotation/results/gm12878_expressed_sensitivity.tsv`.
- 11,076 of the 19,229 MANE genes have TPM >= 1; 7,823 have TPM >= 10.

MANE Select, gene TPM >= 1 (>= 10):

| | assembled | exact | compatible | exact given assembled |
|---|---|---|---|---|
| LRAA | 0.544 (0.621) | 0.510 (0.580) | 0.655 (0.731) | 0.937 |
| StringTie | 0.796 (0.874) | 0.741 (0.814) | 0.879 (0.898) | 0.931 |
| Trinity GG | 0.688 (0.777) | 0.625 (0.707) | 0.892 (0.913) | 0.909 |

**Reading**
- For expressed genes, short-read RNA-seq from one cell line plus TransDecoder gives the exact MANE CDS for 74-81% (StringTie), and a correct full or partial CDS for about 90% (Trinity/StringTie).
- LRAA is lower because its long-read depth is lower; the TPMs come from the short reads.

**Caveat**
- For GENCODE (all isoforms), gene-level TPM does not tell which isoforms are expressed, so the GENCODE rows barely move.
- A transcript-level TPM filter would be the right one there, but isoform-level quantification is itself uncertain.

## Sample-matched Ribo-seq support (GM12878, Cenik et al. 2015, GSE65912)

Data: SRR1803147 and SRR1803148 (CHX, RNase A+T1). Pipeline:
- adapter trim to 25–34 nt;
- rRNA/tRNA removal with bowtie2: 37–48% of reads remain non-rRNA;
- STAR end-to-end, unique alignments only: **11.5M footprints in total** (modest depth).

Scripts:
- `genome_annotation/scripts/gm12878_riboseq_{download,process}.sh`
- `analyses/riboseq/riboseq_orf_support.py`
- `analyses/riboseq/expr_matched_support.awk`

**Calibration on MANE CDSs** (`results/gm12878_riboseq_calibration.tsv`):
- 98.8% of footprints are sense to the CDS.
- Periodicity is weak, as expected for A+T1 digestion. The best frame-0 fraction is about 0.55, for read lengths 27, 29 and 30 only (chance is 0.33).
- The support test is therefore conservative. Of MANE-exact predictions that have ≥10 P-sites, only 64% pass the frame test.

**Footprint support by prediction class** (distinct CDSs; "covered" = at least 10 P-sites; "supported" = frame-0 enriched, binomial p < 0.01):

| Assembly | Class | n | Covered | Supported |
|---|---|---|---|---|
| LRAA | GENCODE exact | 12,021 | 89% | 59% |
| LRAA | Novel (no GENCODE CDS overlap) | 885 | 4.3% | 1.5% |
| StringTie | GENCODE exact | 13,807 | 83% | 54% |
| StringTie | Novel | 4,127 | 1.0% | 0.4% |
| Trinity | GENCODE exact | 13,832 | 84% | 53% |
| Trinity | Other overlap | 30,509 | 69% | 41% |
| Trinity | Novel | 11,235 | 0.6% | 0.2% |

**Expression-matched comparison** (StringTie source-transcript TPM; `results/gm12878_riboseq_stringtie_by_expression.tsv`):
- At TPM 5–20, 89% of exact ORFs are covered, but only 3.4% of novel ORFs (235 ORFs).
- At TPM ≥ 100, 96% of exact ORFs are covered, but only 10% of novel ORFs (30 ORFs).
- So the lack of footprints is not explained by low RNA abundance.
- Most unannotated predictions (mostly 50–199 aa, many on lncRNA-like transcripts) show no ribosome engagement in this cell line. They are best treated as candidate ORFs, not proteins.
- This matches the minimum-length sweep: the precision cost of the 50 aa minimum falls mainly on this class.

Caveats:
- Depth is modest, and the A+T1 periodicity is weak.
- Only uniquely mapped footprints are counted, so repeat-derived ORFs are invisible.
- Single-exon transcripts assembled on the wrong strand would also lack sense footprints.
- A deeper, RNase I Ribo-seq dataset, or MS (Sinitcyn 2023 PXD024364, 534 GB raw), would be needed to call individual novel ORFs as translated.
