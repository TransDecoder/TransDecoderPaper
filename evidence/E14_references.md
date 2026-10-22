# E14. Literature references added on 2026-10-22

## Claims supported
Citations added to the manuscript (`manuscript/references.bib`, 15 new entries). Each was checked
in the paperclip corpus (full text, with line numbers) where available. Bibliographic details were
checked against Crossref (`api.crossref.org/works/<doi>`) for every entry.

## Sources
| key | where cited | what it supports | checked in |
|---|---|---|---|
| Stover2006 | methods_benchmark (annotation sources) | the Tetrahymena Genome Database | paperclip PMC1347417 L1 |
| Sheng2020 | methods_benchmark | the complete *T. thermophila* MAC genome assembly. Ye et al. call it "the latest published MAC genome of Tetrahymena (TGD, http://ciliate.org)" (PMC11754650 L37, ref. 23 = L165) | Crossref; paperclip abstract-only record |
| Ye2025 | methods_benchmark | the gene models (UTRs for 26,687 genes); the download is `Tetrahymena_thermophila_mac_genome_v6_maintranscript.gff3`, attributed to Ye et al. in `benchmark_v2/annotation_downloads.tsv` | paperclip PMC11754650 L5, L137 |
| Haddrill2007 | results_coding_failures (*D. melanogaster* chr4/Y) | the no-crossover region (chiefly chromosome 4) "shows almost no evidence of codon usage bias" and has lower GC3 | paperclip PMC1852418 L17, L19 |
| Campos2012 | same | non-crossover regions show "lower codon usage bias, lower GC content"; the fourth chromosome shows the largest effects; recombination there "minimal or completely absent" | paperclip PMC3318434 L5, L12 |
| Singh2014 | same | reduced codon bias of Y-linked genes | Crossref only; full text not in paperclip. Supported second-hand by Ricchio et al. 2021 (PMC8623413 L12: "two distinctive features of Drosophila Y-linked genes—reduced codon bias and gigantic introns—most likely reflect the lack of recombination", refs 24–27, of which [25] = Singh 2014, L126) |
| Kozak1987 | methods_algorithm (Kozak score) | naming of the start-codon context | Crossref only |
| Pertea2015, Trapnell2010, Haas2003, Wu2005 | methods_algorithm (StringTie, Cufflinks, PASA, GMAP) | tool citations | Crossref only |
| Yu2026 | methods_algorithm (LRAA) | LRAA, bioRxiv preprint (Brian is co-corresponding author) | paperclip bio_1e10d9ed6e09 L17 |
| Martin2011, Langmead2012, Dobin2013 | methods_experimental_evidence (cutadapt, bowtie2, STAR) | tool citations | Crossref only |

## Text added
- results_coding_failures: "Both chromosomes have little or no crossing-over, and their genes have
  reduced codon bias and low GC3, consistent with less efficient selection on codon usage
  [Haddrill2007, Campos2012, Singh2014]." Campos 2012 excluded Y genes. The Y part rests on
  Singh2014, which was not read directly (see above).

## Gaps
- Not read in full: Singh2014 and Sheng2020 (abstract-only records, which could not be opened in
  paperclip), and the classic tool papers (not in the paperclip corpus).
- Some entries list the first author and "others". Fill in full author lists if the journal style
  requires them.
- Yu2026 uses the DOI given by paperclip (10.64898/2026.02.12.705617). Crossref does not resolve it.
- Kozak (1986, Cell) is not cited; add it if the journal prefers the original
  context paper.
- The LaTeX was not compiled locally (no TeX installation). All cite keys resolve to bib entries,
  checked by script.
