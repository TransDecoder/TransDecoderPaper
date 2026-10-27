# Simulated GM12878 Trinity benchmark: report (2026-10-26)

Evidence: `evidence/E30_simulated_trinity_benchmark.md`. Notebook: `notebooks/02_benchmark.html`, section
"Simulated GM12878 Trinity assembly". Figure: `figures/figS_simtrinity.{pdf,png}`. Tables: `data/simtrinity/`.
Not yet in the manuscript.

## Design
- GENCODE v50: coding transcripts with a complete annotated CDS (no `*_NF` tags; 252,930) plus all lncRNAs
  (191,362). Pre-mRNAs (retained intron etc.) excluded.
- Expression: salmon 1.10.3 on the GM12878 RNA-seq (ISR, `--seqBias --gcBias`).
- Reads: Polyester 1.38.0, 52.75M stranded (RF) pairs of 101 nt, `illumina5` errors, the rounded salmon
  counts of the kept transcripts (not scaled up). Two small Polyester patches (fragments < 50 nt redrawn;
  position-101 error profile = position 100).
- STAR (two-pass; run on methods02 for memory) and genome-guided Trinity 2.15 with `--min_kmer_cov 2`,
  `--SS_lib_type RF`; minimap2 alignment of the assembly.
- Each assembled transcript labelled from its alignment and the annotated CDS of compatible source
  transcripts: complete, 5' partial, 3' partial, internal, noncoding (lncRNA only), or excluded
  (ambiguous, short, unresolved). Projected CDSs check out: 33,561 of 33,939 translate to an exact
  substring of the source protein.

## Results
- Strand: R1 antisense to its source in 99.97% of uniquely aligned reads.
- Assembly: 84,972 transcripts. Scored truth 33,939 (complete 28,187; 5' partial 1,454; 3' partial 1,579;
  internal 2,719); negatives 20,815; unresolved 32.6% (mostly UTR-only pieces and isoform mixtures).
- The simulated assembly is more complete than the real GM12878 Trinity assembly (complete 83% vs 68% of
  scored; transcripts starting inside the CDS 7.0% vs 13.6%), so partial CDSs are under-represented.

| method | stop F1 | exact F1 | lncRNAs with a prediction |
|---|---|---|---|
| TransDecoder (release) | 0.916 | 0.844 | 8.1% |
| CPAT | 0.901 | 0.799 | 11.2% |
| RNAsamba | 0.900 | 0.787 | 11.7% |
| GeneMarkS-T | 0.873 | 0.687 | 20.8% |
| Borf | 0.857 | 0.723 | 19.7% |
| CPC2 | 0.844 | 0.781 | 7.4% |
| TransDecoder v6.0.0 | 0.808 | 0.671 | 17.5% |
| longest ORF | 0.779 | 0.571 | 70.6% |
| Prodigal | 0.630 | 0.453 | 51.8% |

- TransDecoder exact sensitivity by type: complete 0.86, 3' partial 0.75, internal 0.69, 5' partial 0.57.
  Prodigal is far better on 5' partials and internals (0.92-0.98) but predicts on half the lncRNAs;
  CPAT, CPC2 and RNAsamba almost never get 5' partials or internals exactly right (<= 0.015).
- By expression, TransDecoder exact sensitivity stays at 0.80-0.88 in every TPM bin; CPAT drops to 0.69
  below 1 TPM.
- TransDecoder misses 2,188 complete CDSs at stop level, 1,753 with no prediction (mostly < 200 codons).
- CodAn was not run (too slow).

## Next
Read-depth ablation (Brian, 2026-10-26): halve the reads repeatedly, reassemble and rescore, to obtain more
partial assemblies; stop when too few transcripts assemble to be informative.
