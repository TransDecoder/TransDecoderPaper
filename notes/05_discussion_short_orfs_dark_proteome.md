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

## Results of the minimum-length sweep (2026-10-20; PyTD f5e6a39, LongOrfs -m 100..50)

Tables: `work/bench_v3/scores/by_length.tsv` (score_by_length.py) and
`work/gm12878_td/eval/{lraa,stringtie}_m*.summary.tsv` (m100 = `*_kozak_add`).

**lncRNA negatives** (human complete + 33,668 GENCODE v50 lncRNAs; mouse similar):

| min aa | 50-99 aa sens | lncRNAs called | precision |
|---|---|---|---|
| 100 | 0.31 | 6,608 (20%) | 0.74 |
| 70 | 0.77 | 13,840 (41%) | 0.58 |
| 50 | 0.86 | 21,416 (64%) | 0.47 |

- At m50: 14,880 lncRNA calls of 50-99 aa vs 430 correct calls of that length (~35:1).
- At m100: 4,139 lncRNA calls of 100-149 aa vs 958 correct ones.

**Random negatives** (i.i.d., GC-matched; 8 organisms):
- Fraction of random sequences with a call, m100 -> m50:
  - human complete 0.006 -> 0.065;
  - S. pombe 0.005 -> 0.159;
  - Tetrahymena 0.014 -> 0.146.
- Short-protein (50-99 aa) sensitivity on fragments rises from 0.04-0.16 to 0.74-0.95.
- Random sequence is a much easier negative than lncRNA.

**GM12878** (m100 -> m50):

| | StringTie | LRAA |
|---|---|---|
| distinct CDSs | 32.7k -> 49.6k | 23.5k -> 26.5k |
| GENCODE-exact gained | +1,786 | +1,813 |
| MANE-exact gained | +112 | +108 |
| CDSs with no GENCODE CDS overlap | 11.8% -> 34.7% (about 13.4k new) | 3.6% -> 10.7% |

The StringTie assembly includes many noncoding / intronic transcripts, so the increase is larger.
