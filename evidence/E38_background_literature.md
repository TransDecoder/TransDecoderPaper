# E38. Literature behind the Background, Abstract and Conclusions drafts (2026-10-31)

## Purpose
The sources and exact passages for each claim in the new Background paragraphs (main.tex), and how the
Abstract and Conclusions drafts were checked. The full texts were read through paperclip
(`/papers/<id>/content.lines`; L = line). Bibliographic details come from Crossref.

## Claims and sources
| Claim in Background | Source (bib key) | Passage |
|---|---|---|
| Long reads capture full-length transcripts | PardoPalacios2024 (PMC11543605) | L60: long-read sequencing captured full-length and novel transcripts |
| Annotation suites translate predicted coding regions before database searches; Trinotate's translation step is TransDecoder | Raghavan2022 (PMC8921630) | L229 (Trinotate: nucleotide inputs translated with TransDecoder; BLAST+ and HMMER/Pfam) |
| Phylogenomic studies of non-model eukaryotes start from predicted ORFs/proteins | MaurerAlcala2024 (PMC11631190) | L23 |
| Ciliates use alternative genetic codes | MaurerAlcala2024 | L14 (seven alternative codes among Ciliophora) |
| Described as one of the most widely adopted ORF predictors for RNA-seq data | MaurerAlcala2024 | L14, L23 ("few have become as well adopted as TransDecoder for RNA-seq based datasets") |
| Lower precision than newer tools across 32 eukaryotes | MaurerAlcala2024 | L24 (TransDecoder 85.0%, lowest of the tools tested) |
| Incomplete transcripts: sequence upstream of the first ATG may be UTR or truncated CDS | Signal2021 (bio_db2e704d6c94) | L12, L44 |
| Adapted to genome-guided assemblies (StringTie, Cufflinks) | Varabyou2023 (bio_9e81382f6bfd; published Nat Comput Sci) | L34 |
| 64% of GENCODE transcripts at MANE loci matched, against 88% for ORFanage | Varabyou2023 | L35-L37 (72,909 transcripts; longest ORF kept, partials removed) |
| Start codons often upstream of experimentally verified starts | Tang2015 (PMC4499116) | L39-L42 (Table 2: 52.2% exact; 570 longer predictions, 568 running off the 5' end) |
| CDSs predicted in StringTie transcripts close to NCBI CDS lengths | Freedman2025 (PMC12047660) | L56, L62 |
| Longest ORF more accurate than TransDecoder in another study | Grau2026 (PMC12910111) | L117 (performance drops with TransDecoder vs the longest-ORF CDS prediction of GeMoSeq; L169 describes it as a longest-ORF heuristic) |

Already in the bibliography and used as before: Grabherr2011, Pertea2015, Yu2026 (LRAA), Miller2022,
Tang2015, Signal2021. The Signal2021 entry's first author was corrected from "Brandon" to "Beth"
(Crossref).

## The "Here we describe" paragraph, Abstract and Conclusions
These restate the manuscript's own Results and Discussion:
- the 2.7-fold speed-up (E08);
- the decoy test's ~3.5-fold fewer lncRNA predictions (Discussion; E20, E22);
- median exact F1 0.912 complete and 0.822 fragments, and classifiers ≤ 0.486 on fragments (Table 1; E07, E25);
- simulated assemblies: TD first at every depth in human and down to 1/8 in A. thaliana, i.e. 10 of 12
  levels ab initio, and with homology at all 12 (E31, E32, E33);
- GM12878 Ribo-seq and MS (E01, E02, E25);
- the ranks (Fig. S4);
- the v7 changes (E23: decoy test, regenerated GC-dependent retention table, presets, deterministic start rule).

## Originality check
- Every new sentence was written from these notes, not adapted from the source text.
- An automated check compared the new Abstract, Background and Conclusions text with the full text of all
  eight cited papers. It found 0 shared runs of six or more words. The check is in the commit message of
  this record.
