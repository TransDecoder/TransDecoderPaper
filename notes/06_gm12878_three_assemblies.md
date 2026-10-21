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
