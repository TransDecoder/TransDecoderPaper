# E30. Simulated GM12878 Trinity assembly: a benchmark on assembled transcripts with known CDSs (2026-10-26)

## Question
The benchmark inputs (E07, E28) are annotated transcripts cut at random. Real assemblies also contain
isoform mixtures, chimeras, UTR-only pieces and assembly errors, and their fragment-length and
truncation profiles come from coverage, not from a cut model. How do the methods do on a real-style
genome-guided Trinity assembly where the true CDS of every contig is known?

Approach: simulate a stranded GM12878-like RNA-seq library from annotated transcripts only, at the
expression levels of the real reads. Then assemble it with the pipeline of the real GM12878 Trinity
assembly (E11), label each contig by where its source CDS falls, and benchmark all methods on it.

## Data
- **Annotation and genome:** GENCODE v50, GRCh38 primary assembly (`$P/data/gm12878/ref/`, as in E03/E11).
  - `gffread -w` of the GTF gives 646,577 transcript sequences: `$P/data/simtrinity/ref/gencode.v50.transcripts.fa`.
- **Real reads:** ENCSR000AED rep 1 (stranded dUTP/RF, PE 101), ENCFF001REK/REJ, 97,548,052 pairs (E11).
- **Transcripts simulated** (`01_gencode_tables.py`; `data/simtrinity/selection_counts.tsv`):
  - **Coding:** of 278,484 `protein_coding` transcripts:
    - 25,227 carry `cds_start_NF`, `cds_end_NF`, `mRNA_start_NF` or `mRNA_end_NF`;
    - 340 fail the complete-CDS check: no ATG 201 (3 selenoproteins), internal stop 100 (84 selenoproteins), no stop 13, length not a multiple of 3: 13.
    - **252,930 are kept.** Their translations equal GENCODE `pc_translations` for 252,929 of them.
  - **Negatives:** all 191,362 `lncRNA` transcripts.
  - Nothing else is simulated: no pre-mRNA, retained-intron, NMD or processed transcripts, pseudogenes or small RNAs.

## Design and commands
Scripts are in `pipeline/simtrinity/`. Large files are in `$P/data/simtrinity/` and `$P/work/simtrinity/`.

1. **Expression** (`02_salmon.sh`): salmon 1.10.3 (`$P/env/trinity215`).
   - Index of the full GENCODE v50 transcriptome, `-k 31 --keepDuplicates`. Without `--keepDuplicates`, salmon keeps only one of identical sequences, which could be a transcript of an excluded biotype.
   - `salmon quant -l ISR --seqBias --gcBias -p 6 --validateMappings`.
   - Detected library type ISR (71.0M ISR vs 1.2M SF / 1.7M SR fragments).
   - 73,798,725 of 97,548,052 pairs assigned (75.65%).
   - Fragment length mean 215.43, sd 24.09 (`aux_info/meta_info.json`; the 10,000-draw `fld.gz` histogram gives 215.44 / 25.22).
2. **Counts** (`03_make_counts.py`): the rounded salmon `NumReads` of each kept transcript.
   - 94,650 kept transcripts have at least 1 pair: 67,309 coding, 27,341 lncRNA.
   - **52,753,990 pairs:** 49,865,569 coding, 2,888,421 lncRNA.
   - Reads that salmon assigned to excluded biotypes are not simulated, so the depth is 71% of the assigned pairs. Per-transcript coverage matches the real library.
   - The transcripts are dealt to 5 chunks with equal read totals (10.45–10.70M pairs each).
3. **Reads** (`06_simulate.sh` → `04_polyester_chunk.R`, `05_polyester_to_fastq.py`):
   - Polyester 1.38.0 in a new environment, `$P/env/polyester`: `micromamba create -c conda-forge -c bioconda bioconductor-polyester r-base`; R 4.3.3, Biostrings 2.70.1.
   - `simulate_experiment_countmat` per chunk: one sample, `paired = TRUE`, `readlen = 101`, `distr = "normal"`, `fraglen = 215.43`, `fragsd = 24.09`, `error_model = "illumina5"`, `strand_specific = TRUE`, `seed = 20261025 + chunk`.
   - Two patches to Polyester internals, applied in the run script's namespace (deviation, see Caveats):
     - fragment lengths below 50 nt are redrawn (the normal draw can give lengths ≤ 0, which crashes `subseq`);
     - read position 101 uses the position-100 profile of `illumina5`, which covers positions 1–100 only (otherwise "NA in probability vector").
   - **Strand:** in Polyester's source (`get_reads`; `reverse_complement` is not called when `strand_specific = TRUE`), mate 1 is the first 101 nt of the fragment in transcript sense, and mate 2 the reverse complement of its last 101 nt. For the RF layout of the real library, **R1 = Polyester mate 2** and R2 = mate 1.
   - FASTQ with constant quality `I` (Polyester writes FASTA).
   - Read names are `c<k>_<n>:<source transcript>:<fragment start>-<end>`. They contain no `/`, so STAR keeps the source in the name.
   - Wall clock 1.2–1.4 h per chunk, 5 chunks in parallel.
4. **Alignment** (`07_star.sh`): STAR 2.7.11b with the real run's index (`$P/data/gm12878/star_stringtie/star_index`, GENCODE v50 junctions, overhang 100) and parameters:
   - `--twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS --limitBAMsortRAM 40000000000`.
   - 6 threads instead of 12.
   - **Run on methods02,** because methods01 had only 26–31 GB of free memory (peak STAR memory 37 GB). Launcher: `$P/work/simtrinity/run_star_methods02.sh`.
   - A first attempt with `--limitBAMsortRAM 8000000000` failed at BAM sorting: the reads are grouped by source transcript, so the sort needs about 30 GB. It was moved to `$P/_trash/simtrinity_star_failed_sortram_20261026/`.
   - 96.84% of the 52,753,990 pairs mapped uniquely (`sim_GM12878.Log.final.out`).
5. **Assembly** (`08_trinity_gg.sh`, methods01): Trinity 2.15.2 (`$P/env/trinity215`, with the ParaFly xargs drop-in of E11):
   - `--genome_guided_bam sim_GM12878.STAR.bam --genome_guided_max_intron 100000 --SS_lib_type RF --min_kmer_cov 2 --max_memory 50G --CPU 8`.
   - 28,494 of 28,494 partitions completed.
   - **84,972 transcripts in 46,339 Trinity genes.**
   - Then `minimap2 -t 4 -ax splice:hq -uf --secondary=no` to GRCh38, sorted and indexed: 84,838 have a primary alignment (99.84%), with 500 supplementary alignments.
   - Unlike the real run (E11), this run used `--min_kmer_cov 2`, as the design specified.
6. **Labelling** (`label_assembly.py`; the docstring has the full rules). Each contig with a primary alignment is compared with all simulated transcripts on its strand.
   - **Compatible:** inside the overlap of the two extents, the contig's exonic bases (alignment deletions bridged) equal the transcript's exonic bases. So the introns are identical, with no retained intron and no alternative splice site. The contig may extend beyond either end of the transcript.
   - **End tolerance:** only when no exactly compatible transcript exists, up to 10 nt of terminal overhang at the contig's own aligned ends is tolerated. These are bases that minimap2 aligned into the intron instead of across the junction.
   - **CDS projection:** the source CDS (GTF CDS + stop codon) is projected through the alignment.
     - `cds_start` is the first whole codon of the reference frame (the ATG if present).
     - `cds_end` is the end of the stop codon if present, otherwise the end of the last whole codon.
     - Type: complete, 5prime_partial (stop only), 3prime_partial (ATG only) or internal.
     - Where the CDS runs into at most 10 soft-clipped or overhanging bases at a contig end, it is extended in frame to the last whole codon of the contig. More than 10 such bases makes the contig unresolved.
   - **Resolution:** several compatible coding transcripts with the same projection (start, end, type) resolve; different projections make the contig **ambiguous**.
   - **Short:** fewer than 50 sense codons gives `short_complete` / `short_partial`.
   - **Frame check:** the projected CDS must translate with no internal stop, be a whole number of codons, and have the ATG and stop that its type requires. Failures become `unresolved:frame_check_failed`.
   - **Negatives (`noncoding`):** contigs compatible only with lncRNAs whose exons overlap no CDS of any GENCODE v50 transcript on the same strand.
   - **Everything else is `unresolved`,** with a reason: unaligned, chimeric (has a supplementary alignment), large_softclip (more than 30 nt), softclip_at_open_cds_end, no_overlap, incompatible, utr_only, lnc_overlaps_cds, frame_check_failed.
   - **Outputs** in `$P/data/simtrinity/labels/sim_trinity.*`: `labels.tsv` (all contigs), `truth.tsv`, `negatives.txt`, `excluded.txt`.
   - The truth file is in the format of `benchmark_v2/scripts/score_predictions.py`, plus columns `type`, `source_transcripts`, `gene_name`, `tpm` (sum of the salmon TPMs of the source transcripts that give the projection) and `cds_frac_present`.
   - Negatives are absent from the truth file, so every prediction on them is a false positive, as with the `lnc_*` and `rand_*` negatives of the benchmark.
7. **Benchmark input:** `$P/data/simtrinity/inputs/fragments_simtrinity.fa` (all 84,972 contigs, including unresolved and ambiguous; md5 cad8ac11…) and `fragments_simtrinity.truth.tsv` (md5 cf6de6e6…).
   - The `fragments` prefix makes `run_tool.sh` choose the partial-length models (RNAsamba partial weights; CodAn partial models, not used). The input is mostly partial: 5,752 of 33,939 truth CDSs are partial.
   - No `unstranded` in the name, so the input is treated as stranded: TD `-S`, Prodigal plus-strand only, GMST `--strand direct`.
   - Organism `human`: CPAT human model with cutoff 0.364.
8. **Runs** (`$P/work/simtrinity/launch_simtrinity_runs.sh`, 4 jobs at a time on methods01; run directories `$P/work/simtrinity/runs/<tool>/`):
   - `run_tool.sh` (unmodified) for pytd_v7 (PyTD 3a793ff), perl_td_v600, prodigal, gmst, cpat, cpc2, rnasamba and borf.
   - The longest-ORF baseline is `pytd_v7/longest_orf.preds.tsv`, as in E28.
   - CodAn was skipped (too slow). perl_td_v7 was skipped: it is identical to pytd_v7 on every other input (E08, E25).
9. **Scoring** (`score_simtrinity.py`):
   - Predictions on excluded contigs (unresolved, ambiguous, short) are removed. They count as neither true nor false positives and are reported separately.
   - The rest are scored with `score_predictions.py`: stop within 2 nt, exact when the start is also within 2 nt, one credit per truth CDS.
   - Transcript IDs are cut at the first whitespace, because GeneMarkS-T reports the whole FASTA header.
10. **Summaries** (`09_summarize.py`; `10_incompatible_breakdown.py`).
    - The real GM12878 Trinity GG assembly (E11) was labelled with the same script against the same transcript set, for comparison.
    - The E29 statistic was recomputed with `pipeline/derive/03_utr_gm12878_truth.py` on both assemblies.
11. **Notebook:** `notebooks/02_benchmark.Rmd`, section "Simulated GM12878 Trinity assembly", and `figures/figS_simtrinity.{pdf,png}`.

## Software
- salmon 1.10.3, Trinity 2.15.2, bowtie2 2.5.5, jellyfish 2.3.1 (`$P/env/trinity215`).
- STAR 2.7.11b, minimap2 2.31-r1302, samtools 1.24, gffread 0.12.9 (`$P/env/rnaseq_asm`).
- Polyester 1.38.0, R 4.3.3, Biostrings 2.70.1, IRanges 2.36.0 (`$P/env/polyester`).
- Python 3.11.17 (`$P/env/td`).
- Predictors as in E07/E28: PyTD 3a793ff, TransDecoder v6.0.0 (66d4712), Prodigal 2.6.3, GeneMarkS-T 5.1, CPAT 3.0.5, CPC2 1.0.1, RNAsamba 0.2.5, Borf 1.2.

## Results
**Strand.**
- Transcript level: 20,000 pairs from chunk 0 were mapped to their source transcripts (`minimap2 -ax sr`). Of the reads on the source transcript, R1 was reverse in 15,936/15,936 and R2 forward in 15,936/15,936 (`data/simtrinity/strand_check_transcripts.tsv`).
- Genome level: in a 1% sample of the STAR BAM, uniquely mapped R1 lies on the strand opposite its source transcript's gene in 510,542 of 510,704 reads (99.97%) (`strand_check_genome.tsv`).

**Labels** (`label_counts.tsv`):

| label | simulated | real GM12878 Trinity GG |
|---|---|---|
| complete | 28,187 | 15,489 |
| 5' partial | 1,454 | 2,283 |
| 3' partial | 1,579 | 2,329 |
| internal | 2,719 | 2,714 |
| **scored coding (truth)** | **33,939** | 22,815 |
| noncoding (negatives) | 20,815 | 13,316 |
| ambiguous | 590 | 805 |
| short_complete / short_partial | 498 / 1,464 | 303 / 2,016 |
| unresolved | 27,666 (32.6%) | 389,705 (90.9%) |
| total | 84,972 | 428,960 |

- Complete CDSs are 83.1% of the scored coding contigs in the simulated assembly vs 67.9% in the real one.
- **Simulated unresolved by reason:**
  - utr_only 13,447: compatible with coding transcripts but holding only UTR sequence;
  - incompatible 12,137;
  - no_overlap 595;
  - large_softclip 483;
  - chimeric 464;
  - softclip_at_open_cds_end 283;
  - unaligned 134;
  - lnc_overlaps_cds 92;
  - frame_check_failed 31.
- **The incompatible contigs** (`incompatible_breakdown.tsv`, by closest simulated transcript):
  - 6,033 have only annotated introns but in a combination no simulated transcript has (assembler isoform mixtures);
  - 3,322 have the right introns but an exon end more than 10 nt off;
  - 1,964 contain an intron absent from every overlapping simulated transcript;
  - 818 differ by at most 10 bases at places the end tolerance does not cover.
- **Strictness check:** without the end tolerance and the in-frame extension through short clips (the first labelling pass), the simulated assembly had 28,893 unresolved contigs. The counts were complete 28,119, 5' partial 1,351, 3' partial 1,390, internal 2,369, negatives 20,409 (`label_counts_simulated_no_end_tolerance.tsv`). With both rules: 27,666 unresolved; +68 complete, +103 / +189 / +350 partials and +406 negatives. Most remaining unresolved contigs are UTR-only pieces or genuine assembly recombinations, not labeller artefacts.
- **Truth check:** the projected CDS translates to a substring of the first source transcript's protein in 33,561 of 33,939 truth CDSs. In 169 more, at most 1% of amino acids differ (sequencing errors that survived assembly). In 209 (0.6%) it differs more; only the first-listed source was compared, and these cases were not examined further.

**E29 statistic** (contigs carrying a MANE Select stop with the reference structure to their 5' end; `e29_five_prime.tsv`):
- simulated: 1,116 / 16,001 = 7.0% begin inside the CDS (median missing fraction 0.78);
- real Trinity: 1,915 / 14,087 = 13.6% (reproduces E29).

The simulated assembly is more complete than the real one. It has uniform fragment starts, no pre-mRNA and no degradation, and its contigs are about 5× fewer.

**Complete CDSs by expression:**
- Scored contigs with source TPM < 1 / 1–5 / 5–20 / 20–100 / > 100 are complete in 70.4% / 95.5% / 96.2% / 95.1% / 82.6% (`type_by_expr.tsv`).
- Simulated coding transcripts with at least one complete contig: 19.8% / 71.6% / 84.2% / 88.3% / 82.9% (`complete_by_source_expr.tsv`; 44,621 / 14,096 / 5,671 / 2,349 / 572 transcripts).

**Accuracy, all scored contigs** (33,939 truth CDSs, 20,815 negatives; `scores_overall.tsv`, `negatives.tsv`):

| method | stop sens | stop prec | stop F1 | exact sens | exact prec | exact F1 | negatives with a prediction | predictions on excluded contigs |
|---|---|---|---|---|---|---|---|---|
| TransDecoder (pytd_v7) | 0.899 | 0.934 | **0.916** | 0.828 | 0.861 | **0.844** | 8.1% | 9,918 |
| TransDecoder v6.0.0 | 0.893 | 0.738 | 0.808 | 0.742 | 0.613 | 0.671 | 17.5% | 14,387 |
| longest ORF | 0.948 | 0.661 | 0.779 | 0.695 | 0.485 | 0.571 | 70.6% | 21,392 |
| Prodigal | 0.978 | 0.465 | 0.630 | 0.704 | 0.334 | 0.453 | 51.8% | 28,819 |
| GeneMarkS-T | 0.890 | 0.856 | 0.873 | 0.700 | 0.674 | 0.687 | 20.8% | 11,751 |
| CPAT | 0.888 | 0.914 | 0.901 | 0.788 | 0.811 | 0.799 | 11.2% | 10,260 |
| RNAsamba | 0.902 | 0.898 | 0.900 | 0.789 | 0.786 | 0.787 | 11.7% | 10,579 |
| CPC2 | 0.776 | 0.925 | 0.844 | 0.718 | 0.856 | 0.781 | 7.4% | 8,911 |
| Borf | 0.865 | 0.849 | 0.857 | 0.730 | 0.716 | 0.723 | 19.7% | 11,656 |

**Sensitivity by type**, stop / exact (`scores_by_type.tsv`):

| method | complete (28,187) | 5' partial (1,454) | 3' partial (1,579) | internal (2,719) |
|---|---|---|---|---|
| TransDecoder | 0.922 / 0.859 | 0.717 / 0.574 | 0.819 / 0.751 | 0.800 / 0.689 |
| TransDecoder v6.0.0 | 0.943 / 0.780 | 0.672 / 0.607 | 0.772 / 0.542 | 0.569 / 0.537 |
| longest ORF | 0.949 / 0.665 | 0.875 / 0.875 | 0.908 / 0.535 | 1.000 / 1.000 |
| Prodigal | 0.981 / 0.672 | 0.924 / 0.922 | 0.958 / 0.586 | 0.985 / 0.984 |
| GeneMarkS-T | 0.939 / 0.727 | 0.682 / 0.649 | 0.762 / 0.528 | 0.564 / 0.551 |
| CPAT | 0.943 / 0.901 | 0.576 / 0.011 | 0.868 / 0.827 | 0.505 / 0.013 |
| RNAsamba | 0.941 / 0.900 | 0.583 / 0.012 | 0.892 / 0.849 | 0.672 / 0.015 |
| CPC2 | 0.889 / 0.849 | 0.472 / 0.010 | 0.286 / 0.272 | 0.052 / 0.001 |
| Borf | 0.924 / 0.806 | 0.565 / 0.252 | 0.744 / 0.612 | 0.480 / 0.265 |

**Per-type precision and F1 (added 2026-10-27, E31).** E31 scores these at full depth too (level d1 in `data/simtrinity/depth/scores_by_type.tsv`), in two versions:
- precision over the predictions on truth transcripts of each type;
- the same, but also counting predictions on negatives whose own ORF is of that type (`*_wneg` columns).

TransDecoder's exact F1 with negatives counted, vs the best other method:

| type | TransDecoder | best other method |
|---|---|---|
| complete | 0.869 | CPAT 0.888 |
| 5' partial | 0.599 | GeneMarkS-T 0.455 |
| 3' partial | 0.776 | CPAT 0.796 |
| internal | 0.747 | Prodigal 0.872 |

**Sensitivity by source TPM**, stop / exact (`scores_by_expr.tsv`; 16,562 / 9,856 / 4,697 / 2,175 / 649 truth CDSs):

| method | < 1 | 1–5 | 5–20 | 20–100 | > 100 |
|---|---|---|---|---|---|
| TransDecoder | 0.884 / 0.801 | 0.904 / 0.840 | 0.937 / 0.880 | 0.925 / 0.874 | 0.860 / 0.818 |
| TransDecoder v6.0.0 | 0.853 / 0.721 | 0.931 / 0.763 | 0.953 / 0.783 | 0.924 / 0.752 | 0.807 / 0.606 |
| Prodigal | 0.975 / 0.753 | 0.978 / 0.665 | 0.986 / 0.652 | 0.980 / 0.629 | 0.992 / 0.669 |
| GeneMarkS-T | 0.849 / 0.685 | 0.928 / 0.713 | 0.948 / 0.731 | 0.923 / 0.719 | 0.814 / 0.623 |
| CPAT | 0.846 / 0.693 | 0.922 / 0.870 | 0.949 / 0.906 | 0.935 / 0.893 | 0.858 / 0.777 |
| RNAsamba | 0.870 / 0.693 | 0.926 / 0.873 | 0.950 / 0.904 | 0.938 / 0.886 | 0.895 / 0.790 |
| CPC2 | 0.693 / 0.618 | 0.858 / 0.812 | 0.886 / 0.850 | 0.842 / 0.810 | 0.624 / 0.592 |
| Borf | 0.818 / 0.667 | 0.909 / 0.784 | 0.933 / 0.816 | 0.898 / 0.794 | 0.755 / 0.670 |
| longest ORF | 0.946 / 0.743 | 0.944 / 0.660 | 0.964 / 0.647 | 0.945 / 0.614 | 0.948 / 0.635 |

**Predictions on excluded contigs** (`excluded_preds.tsv`), as the fraction with a prediction, TransDecoder:
- utr_only 1.6%;
- incompatible 68.0%;
- ambiguous 80.3%;
- short_complete 18.9%; short_partial 3.2%.

Prodigal predicts on 56.5% of all unresolved contigs and the longest-ORF baseline on most of them.

**Predicted ORF types** (scored contigs, from each prediction's own ends; `pred_types.tsv`), as complete / 5' / 3' / internal:
- TransDecoder 82.6 / 5.8 / 5.4 / 6.2%;
- Prodigal 55.6 / 34.2 / 4.9 / 5.4%;
- longest ORF 54.6 / 32.9 / 4.0 / 8.6%;
- CPAT, CPC2 and RNAsamba report no 5'-open ORFs.

**Observations**
- On this assembly TransDecoder has the highest stop F1 (0.916) and exact F1 (0.844). CPAT (0.901 / 0.799) and RNAsamba (0.900 / 0.787) follow.
- Prodigal and the longest-ORF baseline find the most stops: sensitivity 0.978 and 0.948, the best on 5' partials and internals. But they predict on 52% and 71% of the lncRNA-derived negatives, so their precision is 0.47 and 0.66. This repeats the random/lncRNA-negative result of E04/E07 on assembled sequence.
- **Misses on complete CDSs.** TransDecoder's complete-CDS stop sensitivity (0.922) is lower than on the benchmark's complete inputs. Of the 2,188 complete CDSs it misses at stop level:
  - 1,753 have no prediction, mostly short CDSs: 1,100 under 100 codons, 633 of 100–199;
  - 435 have a prediction in another frame or ORF.
  - GENCODE includes many short alternative-isoform CDSs.
- **Exact matches on partial CDSs.** On 5' partials and internals the ATG-requiring classifiers almost never get the open start (exact ≤ 0.015). TransDecoder leaves the 5' end open on 0.574 of 5' partials (0.689 internals); Prodigal on 0.922 (0.984). The same pattern appears in E28.

## Checks
- Strand: transcript-level and genome-level checks above, both RF.
- Truth: every truth CDS passed the frame check (whole codons, ATG/stop as the type requires, no internal stop). 99.4% of them match the source protein exactly or within 1% of amino acids.
- The labels partition the input: truth + negatives + excluded = 84,972 = FASTA count (asserted in `score_simtrinity.py`).
- E29 statistic on the real assembly reproduced exactly (1,915 / 14,087).
- PyTD worktree at 3a793ff (`pytd_version.txt`).
- After the notebook knit, the checks section of notebook 02 still shows 0 mismatches (see the notebook).

## Caveats
- **Depth** is that of the kept transcripts (52.75M pairs, 71% of salmon-assigned pairs). Reads that salmon assigned to excluded biotypes (pre-mRNA-like retained-intron, NMD, other) were not simulated. Real assemblies are noisier: 91% of the real contigs are unresolved against this transcript set.
- **Polyester:**
  - uniform fragment starts (no positional bias) and substitution-only errors;
  - two internal patches (fragment length ≥ 50 nt, error profile of position 100 reused at 101);
  - constant base qualities;
  - reads are written grouped by source transcript.
- **Expression bins** use the summed salmon TPM of the source transcripts that give the projection, not the reads in the contig.
- **Unresolved contigs (32.6%)** are excluded from the truth and from the false-positive count. Predictions on UTR-only contigs (13,447; TransDecoder predicts on 1.6%, Prodigal and the longest ORF on more) are not counted as false positives, although these contigs contain no annotated CDS.
- **The labeller's end tolerance (10 nt) and in-frame extension through ≤ 10 clipped or overhanging bases** were added after inspecting the first, strict labelling pass (see Results). Rule and counts are reported for both passes.
- **Host:** STAR ran on methods02 (memory), with the same index and parameters (6 threads). Everything else ran on methods01.
- CodAn and perl_td_v7 were not run (see Design).
- **Single simulation (one seed set).** No replicate simulations were made.
