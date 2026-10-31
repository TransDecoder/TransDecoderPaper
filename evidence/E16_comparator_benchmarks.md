# E16. How comparator papers and third-party studies benchmarked (literature review, 2026-10-23)

## Purpose
Background for the Background and Discussion sections and for the benchmark design. It records how the
comparator tools were evaluated, and which published studies ran an earlier TransDecoder version and
what they found. Nothing here is cited in the manuscript yet; every reference must be checked against
its paper before it is cited (E14).

## Sources
Four review agents on 2026-10-23 read full text through paperclip (PMC `content.lines`; L = line in
that file) and supplements, or WebFetch where noted. Their copies and downloads are in the main session
scratchpad (`papers/`, `dl/epmc/`, `codan_S2.tsv`; not in git).

## Comparator tools' own benchmarks
| Tool (paper) | Task scored | Data | Negatives | Match criterion | Statistics | TransDecoder compared? |
|---|---|---|---|---|---|---|
| CodAn (Nachtigall 2021, PMC8138839) | CDS | 34 RefSeq species, 2,000 full-length transcripts each; strand-blind; simulated NoStart/NoStop/both; Ribo-seq-validated starts (Lee 2012, Lim 2018) | 3'UTRs; Rfam ncRNA > 200 nt | stop on full-length sets, whole CDS on Ribo-seq and partial sets; per transcript P/R/F1 | mean ± SD across species; no tests | **yes, v5.5.0, defaults (no homology, -m 100, apparently no -S)** |
| GeneMarkS-T (Tang 2015, PMC4499116) | CDS | RefSeq NM_ (S. pombe, mouse, fly, Arabidopsis); partials trimmed by an empirical RGASP coverage profile; real fly RGASP assemblies; 1,392 mouse Ribo-seq-validated starts | none (argue none exist) | stop = TP; starts scored separately; "Sp" is precision | none; 50 resamples for training-size curves | **yes, SourceForge build (no version), defaults, min length swept 90-480 bp, stranded and not** |
| CPAT (Wang 2013, PMC3616698) | coding/noncoding class | human RefSeq + GENCODE/Cabili lncRNA | lncRNA | class only | chi-square model test; 10-fold CV | no |
| CPC2 (Kang 2017, PMC5793834) | class | 6 species, RefSeq/Swiss-Prot positives (CD-HIT 0.9), Ensembl ncRNA | ncRNA | class only | none | no |
| RNAsamba (Camargo 2020, PMC7671399) | class | published test sets, cross-species, sORF sets, truncated-CDS sets (MMseqs2 homology filtering train/test) | lncRNA | class only | means of 5 runs; Mann-Whitney for Kozak | no |
| Prodigal (Hyatt 2010, PMC2848648) | genes | prokaryotes; 2,443 experimentally verified starts | none | stop; start+stop | none | no (predates) |
| ANGLE (Shimizu 2006) | coding regions | simulated frameshift/substitution errors (abstract only; full text not read) | ? | MCC | ? | no (predates); ANGEL's fallback reimplements TD's longest-ORF + hexamer scheme (`src/DumbORF.py`; issue #2) |
| Borf (GitHub, no journal paper) | ORF | README: none | - | - | - | see bioRxiv below |

**CodAn results for TD.** Strand-specific stop-level P/R/F1 0.70/0.69/0.70 in vertebrates (Table 2; human 52.5% correct stops, Table S2_A); 0 correct on NoStart sets for every tool (Table S2_F), which suggests a scoring artifact for 5'-open output. The tables have internal errors: in Table 3, P = R = 0.18 with F1 0.93; Table S2_D contradicts Table 3; Table S2_E has scrambled columns. On 3'UTR negatives TD had the best specificity outside vertebrates (L180). Runtime on 2,000 transcripts: TD 9 s, CodAn 32 s.

**CodAn training species (Table S1_A).** VERT: human, mouse, chicken, anole, zebrafish, salmon. INV: C. elegans, **D. melanogaster**, honeybee, Nematostella. PLANTS: **A. thaliana** plus 5 crops. FUNGI: **S. pombe** plus 5. Five of our eight test organisms (human, mouse, fly, Arabidopsis, S. pombe) are CodAn training species. There are no protist or algal models.

**GeneMarkS-T results for TD.** Mouse Ribo-seq start test, minimum length 300 for all tools (Table 2): exact 5' and 3' 52.2% (GMST 68.6%, Prodigal 51.3%); 568 of TD's 570 "longer" predictions run off the 5' end. "TransDecoder adopts the 'longest-ORF' rule and selects the 5'-most AUG" (L73). On RGASP assemblies, TD had 2-3x more transcripts with more than one CDS (Table S2) and more FPs (Table 3).

## Third-party benchmarks that included TransDecoder
| Study | TD version / settings | Data and metric | Finding about TD |
|---|---|---|---|
| ORFanage, Varabyou 2023 Nat Comput Sci, doi:10.1038/s43588-023-00496-1 | defaults, longest ORF kept, partials removed | GENCODE v41 MANE loci; exact ORF | 64% exact vs 88%; no ORF in 16,484; reference-guided tool, "different niche" |
| TIdeS, Maurer-Alcala & Kim 2024 GBE, doi:10.1093/gbe/evae252 | v5.5.0 + DIAMOND Swiss-Prot; ciliate code | rnaSPAdes assemblies of 32 eukaryotes; frame and coverage via blastx | lowest precision (85.0%); full-length recall 55.8% vs 80.2% (TIdeS) and 73.9% (GMST) |
| Borf, Signal & Kahlke 2021 bioRxiv, doi:10.1101/2021.04.12.439551 | defaults, both strands | stranded Trinity (human, zebrafish, Arabidopsis, yeast) | up to 15% antisense picks; 5' extension mislabels ORF class; 11-23% of assembled transcripts lack an upstream in-frame stop |
| TranSuite, Entizne 2020 bioRxiv, doi:10.1101/2020.12.15.422989 | genome workflow, --single_best_only | AtRTD2 | 26.9% differ; > 9,000 long downstream ORFs past PTCs; 5,385 non-AUG 5'-extended starts |
| TCW v4, Soderlund bioRxiv, doi:10.1101/733311 | v5.5.0 + Swiss-Prot hits, single best | 15,000 rice/lotus transcripts with perfect TrEMBL hits | wrong starts 2,022 vs 740; 0 wrong ends |
| CodingQuarry, Testa 2015 BMC Genomics, doi:10.1186/s12864-015-1344-4 | Cufflinks transcripts + Pfam | S. pombe, S. cerevisiae; gene-level Sn/Sp | 80.2/73.5 vs 91.3/89.0; "performs extremely well" for its purpose |
| inGAP-CDG, Peng 2016 Genome Biol, doi:10.1186/s13059-016-1094-x | not stated | human Trinity vs CCDS; orthologs | higher Sn, lower Sp, redundant, shorter CDSs |
| GeMoRNA, Grau & Keilwagen 2026 NAR, doi:10.1093/nar/gkag091 | v5.7.0 -S, no homology | 7 species, StringTie/Scallop/Cufflinks; CDS F1 | worse than the plain longest ORF ("Surprisingly...", L103) |
| Freedman & Sackton 2025 Genome Res, doi:10.1101/gr.280377.124 | v5.5.0 + BLASTP | 21 species; BUSCO-type completeness, PSAURON | near-NCBI CDS lengths; ORF errors at transcript termini |
| TD2, Mao ... Haas ... Sommer 2025 bioRxiv, doi:10.1101/2025.04.13.648579 (not independent) | defaults | MANE, C. elegans, Arabidopsis, Plasmodium; LNCipedia negatives | TD2 higher sensitivity, lower lncRNA FP |

Recurring criticisms are 5'-end and start errors (extension to the transcript end), antisense and
multiple ORFs, longest-ORF bias in PTC isoforms, false positives at short lengths, and speed. No study
compared ORF predictors using Ribo-seq or MS evidence from real assemblies (GMST and CodAn used Ribo-seq
starts on reference transcripts only).

## Comparison with our design (E07)
- **Consistent with the field:** stop-level and exact start+stop matching, per CDS (Prodigal, GMST, CodAn); several species; complete and simulated partial inputs; strand-blind variants; runtime.
- **Beyond them:**
  - every tool run on the same inputs with its organism-appropriate options;
  - GC-matched random negatives and GENCODE lncRNA negatives;
  - ciliate code and GC extremes;
  - one ORF per transcript enforced uniformly;
  - Ribo-seq and MS evidence on real assemblies (GM12878).
- **Gaps relative to them:**
  - (i) no experimentally validated start set: GMST and CodAn used Ribo-seq TIS;
  - (ii) no uncertainty estimates (none of the papers has CIs either);
  - (iii) fragments are cut once at a uniform random point rather than by an empirical assembly profile (GMST), and real de novo assemblies are scored only for human (GM12878);
  - (iv) one transcript per gene, so non-canonical isoforms (PTC, retained intron; TranSuite's ldORF issue) are untested outside GM12878;
  - (v) CodAn's training overlaps 5 of our 8 test species;
  - (vi) operating points are fixed, whereas GMST swept the minimum length for every tool.
