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

## Ribo-seq ORF catalog comparison (2026-10-20)

Catalogs: GENCODE phase 1 (Mudge et al. 2022; 7,264 ORFs, peptide-support tiers incl. HLA) and
phase 2 (TransCODE 2026; 28,359 ORFs), in `data/riboseq_orfs/`. Matching is by stop codon
(`analyses/riboseq/match_riboseq.py`).

**Catalog reach**
- Phase 1: median 34 aa; 29% are >= 50 aa and 5% are >= 100 aa.
- Phase 2: median 20 aa; about half start at near-cognate codons.
- Most translated sORFs are beyond any ATG-start, length-floored caller.

**Human lncRNA benchmark calls**
- Only 1-1.5% match a Ribo-seq ORF stop, at every length cutoff.
- Of the 495 catalog ORFs that are findable in principle (ATG start, >= 50 aa, inside the chosen
  transcript), TransDecoder calls 25% at -m 50 and 45% of the 83 that are >= 100 aa.

**GM12878 unannotated (no GENCODE CDS) predictions**
- Ribo-seq support: 6% (LRAA) vs 2% (StringTie).
- At -m 50, 148 (LRAA) / 157 (StringTie) predictions match ORFs with HLA peptide support.

## Short-ORF frame-contrast guard (PyTD 7e9ffa2, not pushed)

The rule: ORFs shorter than 100 codons need frame contrast >= 200 x the median per-base coding
contrast, unless they have homology support.
- Summary table: `work/bench_v3/scores/short_orf_summary.tsv`.

**Human lncRNA, -m 50 + guard 200 vs -m 100**

| | short sens | lncRNA calls | precision |
|---|---|---|---|
| complete | 0.31 -> 0.56 | 6,608 -> 6,943 | 0.737 -> 0.729 |
| fragments | 0.13 -> 0.52 (all-length sens 0.80 -> 0.87) | 6,508 -> 6,833 | 0.700 -> 0.706 |

Mouse is similar.

**Comparison with lowering the cutoff alone**
- Guard 150 reaches the short-protein sensitivity of -m 80 with about 2,900 fewer lncRNA calls.
- Random-negative calls stay near the -m 100 level for all 8 organisms.
- In AT-rich P. falciparum the guard is strict: short sensitivity 0.11 vs 0.03 at -m 100.

**GM12878, -m 50 + guard 200 vs -m 100**

| | GENCODE-exact gained | MANE-exact gained | no-GENCODE-overlap |
|---|---|---|---|
| LRAA | +679 | +44 | 3.6% -> 3.7% |
| StringTie | +1,001 | +66 | 11.8% -> 12.0% |

**Key point for the Discussion:** the guard also removes most Ribo-seq-supported lncRNA ORFs.
- Ribo-seq matches drop from 458 to 178 (LRAA) and 578 to 217 (StringTie).
- HLA-supported matches drop from 148 to 37 (LRAA).
- Translated noncanonical sORFs lack the codon-usage signal that composition models detect.
- So sequence-only prediction can either annotate canonical short proteins precisely or capture
  the dark proteome, not both. Translation evidence (Ribo-seq, proteomics) is required for the latter.

## lncRNA negatives: comparator tools (2026-10-20; final defaults)

Table: `work/bench_v3/scores/lnc_by_era.human.tsv` (`analyses/lnc_holdout/lnc_calls_by_era.py`).

Human complete.wLnc (33,668 lncRNAs):

| tool | lncRNAs called | sensitivity (stop) |
|---|---|---|
| TransDecoder | 20.6% | 0.98 |
| CPAT | 10.5% | 0.98 |
| CPC2 | 7.9% | 0.96 |
| RNAsamba | 9.1% | 0.99 |
| Borf | 23.3% | 0.96 |
| GeneMarkS-T | 24.6% | 0.98 |
| Prodigal | 35.3% | 0.98 |

**Training-leakage check:** lncRNA genes in GENCODE v28 vs genes added later.
- RNAsamba: 6.2% vs 10.9%. Its training lncRNAs inflate its apparent precision.
- CPAT: 10.9% vs 10.2%. CPC2: 7.9% vs 7.9%. Neither is affected.

**Fragments**
- TransDecoder has the highest sensitivity: 0.86, vs CPAT 0.79, RNAsamba 0.81, CPC2 0.42.
- It has a much higher exact F1: 0.68, vs 0.47 or less.

**Discussion point**
- Transcript-level classifiers trained on lncRNA vs mRNA reject noncoding transcripts better than
  TransDecoder's ORF-level Markov test.
- But they do not locate the CDS in partial transcripts.
- Future direction: add a transcript-level classifier, or features learned from lncRNAs, as an
  optional filter.
