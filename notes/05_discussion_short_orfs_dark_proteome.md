# Discussion threads: short ORFs, false positives, and the dark proteome (2026-10-20)

From Brian's direction; to be developed into Results (length sweep) and Discussion.

## Minimum ORF length (Results)
- TransDecoder's default minimum is 100 aa. Lowering it toward 50 aa recovers short annotated
  proteins (about 30% of the fully assembled MANE CDSs missed on GM12878 are < 100 aa), but
  false positives rise steeply as the minimum falls: short ORFs arise by chance far more often,
  and their coding scores carry less evidence.
- Planned analysis: sweep the minimum length over 100, 90, 80, 70, 60 and 50 aa.
  - Benchmark: complete.wRand and fragments.wRand inputs for all 8 organisms. Report sensitivity
    by true protein length, and false positives on the GC-matched random negatives.
  - Possibly also add GENCODE lncRNAs (human, mouse) as realistic noncoding negatives.
  - GM12878 (LRAA, StringTie): short MANE/GENCODE recovery against predictions with no
    annotated-CDS overlap, by predicted length.

## Experimental support for short ORFs (Discussion)
- Sequence-based prediction alone cannot tell translated short ORFs from chance ones at
  acceptable precision. Evidence of translation is needed: Ribo-seq (periodicity, initiation
  profiling), mass-spectrometry proteomics / immunopeptidomics, and emerging single-molecule
  peptide sequencing.
- TransDecoder's homology retention (`--retain-blastp-hits`, `--retain-pfam-hits`) already shows
  the route: external evidence overrides the length/score filters. Experimental translation
  evidence could be used the same way.

## The dark proteome (Discussion)
- Unannotated and noncanonical ORFs (sORFs/smORFs, uORFs, lncRNA-encoded and alternative-frame
  ORFs, novel isoform-specific coding regions) form a "dark proteome" relevant to basic biology.
- In cancer, noncanonical ORFs and tumor-specific transcripts (novel isoforms, fusions,
  aberrant splicing) can produce neoantigens and HLA-presented peptides. These are candidate
  immunotherapy targets: TCR therapies, vaccines, and tumor-specific antigens from
  noncanonical translation.
- Genome-wide CDS annotation of assembled transcripts (the GM12878 demonstration), with
  TransDecoder run at relaxed length thresholds, defines candidate search spaces for
  proteogenomics and immunopeptidomics databases. Precision trade-offs at short lengths then
  matter directly for database size and the false-discovery rate.
- References to add: Ribo-seq sORF catalogs and the GENCODE Ribo-seq ORF set (Mudge et al. 2022);
  noncanonical ORFs as cancer antigens (e.g., Laumont et al.; Ouspenskaia et al. 2022);
  long-read proteogenomics (Miller et al. 2022). Verify the citations before use.
