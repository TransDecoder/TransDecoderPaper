# E32. Simulated Arabidopsis Trinity assembly and its read-depth ablation (2026-10-28)

## Question
Do the conclusions of the simulated human Trinity benchmark (E30) and its read-depth ablation (E31) hold for a plant? We use *A. thaliana*: a compact genome with short introns, few lncRNAs and the annotation of the main benchmark (TAIR12).

Approach: the same design as E30/E31, with Arabidopsis expression levels from a stranded Col-0 leaf library. Approved by Brian on 2026-10-27.

**What this tests (Brian, 2026-10-28).** The noncoding RNAs keep their real expression; they were not boosted. They carry 0.4% of the simulated reads, so only 1,372 negative contigs assemble at full depth, and fewer at lower depth. This assembly therefore tests **frame and end placement on assembled coding transcripts**, not the rejection of noncoding sequence. For noncoding rejection, see the human simulation (E30/E31) and the lncRNA inputs of the main benchmark (E04).

## Data
- **Genome and annotation:** TAIR12, the same files as the main benchmark (`$P/data/annot/athal/`, E07): `Athaliana_TAIR12_genomic_Chr_softmasked.fa.gz`, uppercased, and `TAIR12_genome_annotation_TAIR.gff3.gz`.
- **Transcript set** (`pipeline/simtrinity_athal/01_tair12_tables.py`; `data/simtrinity_athal/selection_counts.tsv`):
  - **All transcripts:** 54,713 GFF3 transcript features (mRNA, ncRNA, pseudogene transcripts, tRNA, rRNA, …). All are in the salmon index.
  - **Coding:** all 36,572 TAIR12 mRNAs pass the E30 rules (ATG start, stop codon end, length divisible by 3, no internal stop).
    - Translation verified for all 36,572.
    - All are on Chr1–5; TAIR12 has no organellar mRNA records.
    - The GFF3 CDS segments already include the stop codon; no stop was added from the genome.
    - TAIR12 has no `*_NF`-type partial flags.
  - **Negatives:** transcripts with ncRNA_class lncRNA (2,618) or antisense_RNA (3,347), keeping those whose exons overlap no CDS on the same strand.
    - Kept: 2,576 lncRNA (42 removed) and 3,342 antisense_RNA (5 removed), 5,918 in total.
    - antisense_RNA is kept and reported separately (Brian, 2026-10-28).
  - Not simulated: pseudogenes, tRNA, rRNA, sno/snRNA and miRNA precursors.
- **Real reads: ENA ERR4140056.**
  - Study PRJEB38272 (WOX transcriptome), sample SAMEA6828490: *A. thaliana* Col-0 wild-type leaf, replicate 3.
  - Library: poly(A), NEBNext Ultra Directional (dUTP, RF), HiSeq 2500, 2 × 100 nt, 30,525,183 pairs.
  - Downloaded from the ENA FTP; md5 matches ENA (R1 bbac8950…, R2 394b6671…).
  - Why this run: untreated wild type; dUTP stranded 2 × 100 HiSeq, the closest match to the GM12878 library (2 × 101 HiSeq dUTP); 30–60 M pairs; published. Candidates with 2 × 150 reads, treatments, mutants or proprietary libraries were not used.

## Design and commands
The E30/E31 pipeline is used with Arabidopsis inputs. New scripts are in `pipeline/simtrinity_athal/`, plus `pipeline/simtrinity/13b_depth_collect.py`, a generalised collector. Run on the E31 levels, it reproduces the committed E31 tables byte for byte. Large outputs are in `$P/data/simtrinity_athal/` and `$P/work/simtrinity_athal/`.
```bash
A=$T/pipeline/simtrinity_athal; S=$T/pipeline/simtrinity
python $A/01_tair12_tables.py $P/data/annot/athal/TAIR12_genome_annotation_TAIR.gff3.gz \
  $P/data/annot/athal/Athaliana_TAIR12_genomic_Chr_softmasked.fa.gz $P/data/simtrinity_athal/ref
bash $A/02_salmon_athal.sh          # salmon 1.10.3: index of all 54,713 transcripts (-k 31 --keepDuplicates); quant -l ISR --seqBias --gcBias -p 6
bash $A/03_star_index_athal.sh      # STAR 2.7.11b genomeGenerate, TAIR12 exon GTF, --sjdbOverhang 100 --genomeSAindexNbases 12
python $S/03_make_counts.py $P/data/simtrinity_athal/salmon/quant $R/kept.tsv $R/kept.fa $P/data/simtrinity_athal/sim 4
bash $A/06_simulate_athal.sh        # Polyester: E30 R script and patches (04_polyester_chunk.R), seeds 20261027+k; RF FASTQ
bash $A/07_star_athal.sh            # STAR two-pass, E30 parameters (--limitBAMsortRAM 40 GB), 6 threads, run locally (peak ~5 GB)
CPU=8 bash $A/run_level_athal.sh d1:1                  # Trinity GG + minimap2 (08_trinity_gg_athal.sh), label + 9 tools + score (12_label_and_run_athal.sh)
bash $A/run_level_athal.sh d2:0.5 / d4:0.25 / d8:0.125 / d16:0.0625 / d32:0.03125   # CPU 4; two lanes at a time
python $S/13b_depth_collect.py --data-root $P/data/simtrinity_athal --work-root $P/work/simtrinity_athal \
  $T/data/simtrinity_athal/depth d1 d2 d4 d8 d16 d32
python $A/14_negatives_by_class.py $T/data/simtrinity_athal/depth/negatives_by_class.tsv d1 d2 d4 d8 d16 d32
```
Details:
- **Expression:** 28,993,357 of 30,525,183 pairs assigned (95.0%).
  - Library type ISR, as expected: 27,262,836 ISR vs 0 ISF.
  - Fragment length mean 340.2, sd 88.1 (human: 215.4 / 24.1).
- **Simulation:** the rounded salmon counts of the kept transcripts give 28,853,869 pairs.
  - 28,728,866 come from 26,281 coding transcripts.
  - Only 125,003 (0.43%) come from 2,898 lncRNA/antisense transcripts.
  - Polyester settings and its two patches are as in E30 (paired 101 nt, `illumina5`, normal fragment lengths, `strand_specific = TRUE`; mates swapped to RF).
- **Alignment:** STAR uses the E30 parameters with a TAIR12 index (junctions from the TAIR12 exon GTF). 98.91% of the pairs map uniquely.
  - Strand check on a 1% sample: uniquely mapped R1 lies opposite its source transcript in 285,249 of 285,249 reads (`data/simtrinity_athal/strand_check_genome.tsv`).
- **Assembly:** genome-guided Trinity 2.15.2 with the E30 parameters (`--genome_guided_max_intron 100000 --SS_lib_type RF --min_kmer_cov 2 --max_memory 50G`), CPU 8 for d1 and 4 for the other levels; then `minimap2 -ax splice:hq -uf --secondary=no` to TAIR12.
  - The intron cap of 100,000 is kept as in E30, although Arabidopsis introns are short.
- **Labels:** `label_assembly.py` with the final E30 rules, on the TAIR12 tables.
  - The E29 statistic (5′-truncation of transcripts carrying a MANE Select stop) is not computed, because there is no MANE for Arabidopsis.
- **Methods:** the same nine methods, via the unmodified `run_tool.sh` with organism `athal`.
  - TransDecoder `-G Universal`; Prodigal table 1; GeneMarkS-T code 1, `--strand direct`.
  - CPAT, CPC2 and RNAsamba use their human models; CodAn not run.
  - Input `fragments_simtrinity.fa`: stranded, partial-length RNAsamba weights.
- **Scoring:** `score_simtrinity.py`, plus per-type precision with negatives counted (E31 rule) in the collector.
- **Depth series:** E31 subsampling (`samtools view --subsample`, seed 20261026) of the STAR BAM (`11_depth_subsample_athal.sh`).
  - The levels are nested: on Chr5, all 1,595,320 d4 read names are among d2's 3,190,247.
  - Observed fractions are within 0.0002 of target.
- **Stop:** after d32. The 5′ partial, 3′ partial and internal counts all fell with the total from d16 to d32 (2,484 → 1,940, 1,516 → 1,099, 9,426 → 8,330), which is the E31 rule that partial counts fall with the total. The type shares had not stabilised (d16 → d32: complete −8.3, internal +8.2 points).

## Software
As in E30/E31: salmon 1.10.3, Polyester 1.38.0 (R 4.3.3), STAR 2.7.11b, Trinity 2.15.2 (ParaFly xargs drop-in), minimap2 2.31-r1302, samtools 1.24, PyTD 3a793ff, TransDecoder v6.0.0 and the same comparator versions.

## Results
**Assembly and labels** (`depth/levels.tsv`, `depth/truth_counts.tsv`). Each class is shown as count (% of scored truth).

| level | pairs | Trinity transcripts | complete | 5′ partial | 3′ partial | internal | truth | negatives | unresolved |
|---|---|---|---|---|---|---|---|---|---|
| d1 | 28,853,433 | 23,947 | 17,317 (82.3) | 1,065 (5.1) | 593 (2.8) | 2,078 (9.9) | 21,053 | 1,372 | 1,139 |
| d2 | 14,429,353 | 22,942 | 15,335 (74.9) | 1,315 (6.4) | 793 (3.9) | 3,044 (14.9) | 20,487 | 990 | 1,021 |
| d4 | 7,217,805 | 22,400 | 12,699 (63.0) | 1,797 (8.9) | 1,031 (5.1) | 4,643 (23.0) | 20,170 | 638 | 1,007 |
| d8 | 3,610,788 | 22,655 | 9,216 (45.0) | 2,410 (11.8) | 1,456 (7.1) | 7,407 (36.2) | 20,489 | 321 | 1,058 |
| d16 | 1,803,816 | 21,228 | 5,780 (30.1) | 2,484 (12.9) | 1,516 (7.9) | 9,426 (49.1) | 19,206 | 154 | 1,065 |
| d32 | 901,520 | 16,147 | 3,170 (21.8) | 1,940 (13.3) | 1,099 (7.6) | 8,330 (57.3) | 14,539 | 74 | 847 |

- At d1, 21,010 Trinity genes assemble, and 23,938 transcripts have a primary alignment.
- **Unresolved is much lower than in human:** 4.8% of contigs at d1, against 32.6% (E30). The main reasons are incompatible structure (577) and UTR-only contigs (348).
- **The scored truth stays near 20,000 down to d8.** The Arabidopsis assembly fragments rather than losing genes.

**Negatives by class** (`depth/negatives_by_class.tsv`):
- At d1 there are 1,014 antisense_RNA, 349 lncRNA and 9 "both" contigs.
- Fraction of negatives with a prediction, antisense / lncRNA:

  | method | antisense | lncRNA |
  |---|---|---|
  | TransDecoder | 5.9% | 5.7% |
  | Prodigal | 16.9% | 40.4% |
  | RNAsamba | 54.5% | 51.6% |
  | longest ORF | 58.2% | 61.0% |
  | CPAT | 6.6% | 4.3% |
  | CPC2 | 2.7% | 2.3% |

- From about d8 there are too few negatives (321, then 154 and 74) for negative-based measures to mean much.

**Overall, stop F1 / exact F1 by level** (d1 / d2 / d4 / d8 / d16 / d32; `depth/scores_overall.tsv`):

| method | stop F1 | exact F1 |
|---|---|---|
| TransDecoder (pytd_v7) | 0.972 / 0.970 / 0.973 / 0.968 / 0.969 / 0.971 | 0.881 / 0.873 / 0.859 / 0.838 / 0.822 / 0.828 |
| TransDecoder v6.0.0 | 0.944 / 0.936 / 0.919 / 0.881 / 0.827 / 0.772 | 0.771 / 0.760 / 0.734 / 0.730 / 0.688 / 0.664 |
| longest ORF | 0.976 / 0.982 / 0.987 / 0.992 / 0.994 / 0.994 | 0.731 / 0.747 / 0.773 / 0.815 / 0.852 / 0.876 |
| Prodigal | 0.938 / 0.945 / 0.955 / 0.963 / 0.970 / 0.974 | 0.679 / 0.699 / 0.731 / 0.779 / 0.824 / 0.854 |
| GeneMarkS-T | 0.950 / 0.941 / 0.922 / 0.882 / 0.827 / 0.771 | 0.764 / 0.764 / 0.759 / 0.747 / 0.716 / 0.679 |
| CPAT (human model) | 0.932 / 0.918 / 0.887 / 0.825 / 0.755 / 0.689 | 0.856 / 0.805 / 0.718 / 0.573 / 0.438 / 0.347 |
| RNAsamba (human model) | 0.954 / 0.953 / 0.940 / 0.910 / 0.879 / 0.851 | 0.838 / 0.784 / 0.689 / 0.540 / 0.403 / 0.316 |
| CPC2 | 0.911 / 0.883 / 0.830 / 0.723 / 0.594 / 0.486 | 0.862 / 0.814 / 0.728 / 0.582 / 0.435 / 0.336 |
| Borf | 0.942 / 0.931 / 0.906 / 0.858 / 0.792 / 0.727 | 0.880 / 0.843 / 0.788 / 0.690 / 0.585 / 0.512 |

Negatives with a prediction (%), same level order:

| method | d1 | d2 | d4 | d8 | d16 | d32 |
|---|---|---|---|---|---|---|
| TransDecoder | 6.0 | 6.1 | 4.5 | 4.7 | 2.6 | 4.0 |
| longest ORF | 58.8 | 56.6 | 54.5 | 52.6 | 43.5 | 40.5 |
| RNAsamba | 53.6 | 52.2 | 49.7 | 47.0 | 40.3 | 39.2 |
| Prodigal | 23.0 | 24.3 | 20.8 | 21.2 | 15.6 | 10.8 |

**By type: exact F1 with negatives counted, TransDecoder vs the best other method** (`depth/scores_by_type.tsv`):

| type | d1 | d4 | d8 | d16 | d32 |
|---|---|---|---|---|---|
| complete | 0.918 vs RNAsamba 0.973 | 0.921 vs RNAsamba 0.982 | 0.916 vs RNAsamba 0.986 | 0.918 vs RNAsamba 0.989 | 0.907 vs RNAsamba 0.988 |
| 5′ partial | 0.617 vs Prodigal 0.893 | 0.677 vs Prodigal 0.941 | 0.688 vs longest ORF 0.967 (Prodigal 0.950) | 0.704 vs longest ORF 0.977 (Prodigal 0.955) | 0.743 vs longest ORF 0.967 (Prodigal 0.950) |
| 3′ partial | 0.773 vs RNAsamba 0.858 | 0.855 vs RNAsamba 0.931 | 0.852 vs RNAsamba 0.951 | 0.853 vs RNAsamba 0.948 | 0.847 vs RNAsamba 0.944 |
| internal | 0.736 vs longest ORF 0.989 (Prodigal 0.975) | 0.757 vs longest ORF 0.998 (Prodigal 0.985) | 0.783 vs longest ORF 0.999 | 0.787 vs longest ORF 1.000 | 0.815 vs longest ORF 1.000 (Prodigal 0.987) |

- d2 values are in the table file.
- With so few negatives, the truth-only values differ from these by at most 0.01 (TransDecoder d1: 0.919 / 0.622 / 0.784 / 0.738).
- TransDecoder exact sensitivity, d1 → d32:

  | type | d1 | d32 |
  |---|---|---|
  | complete | 0.904 | 0.891 |
  | 5′ partial | 0.586 | 0.708 |
  | 3′ partial | 0.729 | 0.809 |
  | internal | 0.699 | 0.796 |

**Human (E30/E31) vs Arabidopsis** (figure `figures/figS_simtrinity_athal.{pdf,png}`, notebook 02):
- **Full-depth composition is similar:** 83.1% complete in human vs 82.3% in Arabidopsis.
  - Arabidopsis leaves fewer contigs unresolved (4.8% vs 32.6%).
  - Its negatives are much fewer (1,372 vs 20,815).
- **TransDecoder, overall:**
  - Stop F1 is stable in both species: 0.916–0.929 in human, 0.968–0.973 in Arabidopsis.
  - Exact F1 at full depth is 0.844 in human and 0.881 in Arabidopsis.
  - It is the best overall exact F1 down to 1/8 in Arabidopsis, and at every level in human.
  - In Arabidopsis the longest-ORF baseline and Prodigal pass it at d16–d32 (0.852 / 0.824 vs 0.822 at d16; 0.876 / 0.854 vs 0.828 at d32). At these levels internal fragments, whose reference frame has no stop, are about half the truth, and only 154 / 74 negatives remain to penalise them.
- **Per type, the same pattern holds in both species:**
  - The ATG-requiring classifiers (CPAT, RNAsamba, CPC2) are best on complete and 3′ partial CDSs but cannot leave a 5′ end open: exact sensitivity ≤ 0.025 on 5′ partials and internals.
  - Prodigal and the longest-ORF baseline lead on 5′ partials and internals.
  - TransDecoder ranks between these two groups: below the classifiers on complete and 3′ partial CDSs, below Prodigal and the longest-ORF baseline on 5′ partials and internals, and above the other group in each case.
- **Differences from human:**
  - In Arabidopsis, Prodigal leads TransDecoder on 5′ partials at every level, even with negatives counted, because there are few negatives to cost it. In human, TransDecoder led there at d1–d4.
  - RNAsamba (human model) predicts on about half of the Arabidopsis negatives at every level; in human the figure was 12–20%.

## Checks
- The labels partition each level's FASTA; `score_simtrinity.py` asserts this.
- `neg_pred_types.tsv`: the per-type counts sum to `n_pred_on_negatives` for every level × tool. `negatives_by_class.tsv`: the classes sum to the level's negatives (d1: 1,014 + 349 + 9 = 1,372).
- `13b_depth_collect.py` run on the E31 human levels reproduces the committed `data/simtrinity/depth/` tables byte for byte.
- After the knit, notebook 02 still shows 0 of 133 mismatches.

## Caveats
- **Few negatives.** By design (real expression), the noncoding test in Arabidopsis is weak, and below d4 it is uninformative.
- **Models and settings carried over from human.** CPAT, CPC2 and RNAsamba use their human models. The Trinity intron cap (100 kb) and the STAR parameters are the human ones.
- **No 5′-truncation statistic.** The E29 statistic needs MANE Select, so it is not available for Arabidopsis.
- **Single library, single simulation and single subsample per level.** Differences of about 0.01 are not interpretable.
- **Earlier progress messages (corrections).** They gave TransDecoder's d1 and d4 negative rates as 6.1% and 4.6%. The values are 6.0% (83 / 1,372) and 4.5% (29 / 638).
