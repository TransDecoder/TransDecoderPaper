# E01. GM12878 Ribo-seq support of predicted ORFs

## Claims supported
Manuscript: Results, "Translation and protein evidence for predicted ORFs in GM12878"; Methods, "Ribosome profiling support".

1. Predictions matching annotated CDSs are broadly covered by ribosome footprints.
   - Exact GENCODE matches: 83–89% have ≥10 P-sites, and 53–59% pass the frame-0 test.
2. The frame test is conservative, because periodicity is weak: the frame-0 fraction is about 0.55 for the read lengths used.
3. Unannotated predictions are almost never covered: 4.3% / 1.0% / 0.6% (LRAA / StringTie / Trinity), and 0.2–1.5% pass the frame test.
4. This is not explained by low RNA abundance. At StringTie transcript TPM 5–20, 89% of exact predictions are covered vs 3.4% of novel ones.

## Data
- **Ribo-seq: Cenik et al. 2015, Genome Res 25:1610 (doi:10.1101/gr.193342.115); GEO GSE65912.**
  - GM12878_Rep1 = SRR1803147 (GSM1609378; 159,425,486 reads).
  - GM12878_Rep2 = SRR1803148 (GSM1609379; 151,278,196 reads).
  - 36-nt single-end reads; cycloheximide; RNase A+T1 digestion.
  - Downloaded 2026-10-21 from the NCBI SRA AWS open-data bucket (`https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>`) and checked with `vdb-validate`.
- **Genome and annotation:**
  - GRCh38 primary assembly with GENCODE v50: `$P/data/gm12878/ref/`.
  - STAR index with GENCODE v50 junctions (overhang 100): `$P/data/gm12878/star_stringtie/star_index`.
  - MANE Select v1.5: `$P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`.
- **rRNA/tRNA contaminant set:** 575 sequences. These are U13369.1 (45S rDNA) plus the GENCODE v50 rRNA, rRNA_pseudogene, Mt_rRNA and Mt_tRNA genes.
- **Predictions:** genome-projected PyTD CDSs, final defaults, ab initio (record E03). These are `$P/work/gm12878_td/{lraa,stringtie,trinity}/final_abinitio/transcripts.fa.transdecoder.genome.gff3`.
- **Prediction classes:** `eval_genome_cds.py` pred_status tables, in `$P/work/gm12878_td/eval/status/<asm>_final_abinitio.pred_status.tsv`.

## Software
cutadapt 5.2; bowtie2 2.5.5; STAR 2.7.11b; samtools 1.24; pysam 0.24.1 (Python 3.13); sra-tools 3.4.1.

Environments are `$P/env/riboseq` and `$P/env/sratools` (STAR from `$P/env/rnaseq_asm`). PyTD is PY_dev 48c870c.

## Commands
```bash
bash $T/genome_annotation/scripts/gm12878_riboseq_download.sh     # .sra -> fastq.gz
bash $T/genome_annotation/scripts/gm12878_riboseq_rrna_index.sh   # rRNA/tRNA bowtie2 index
bash $T/genome_annotation/scripts/gm12878_riboseq_process.sh      # trim, rRNA removal, STAR, merge -> footprints.genome.bam
```
Processing steps (in `gm12878_riboseq_process.sh`):
- `cutadapt -a AGATCGGAAGAGCACACGTCT -O 3 --discard-untrimmed --trim-n -m 25 -M 34`
- `bowtie2 --very-sensitive -L 15 --un-gz` against rRNA/tRNA
- `STAR --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1`

Scoring:
```bash
cd $P; S=work/gm12878_td/eval/status; G=work/gm12878_td
env/riboseq/bin/python $T/analyses/riboseq/riboseq_orf_support.py \
  --bam data/gm12878/riboseq_cenik2015/footprints.genome.bam \
  --mane data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
  --pred lraa=$G/lraa/final_abinitio/transcripts.fa.transdecoder.genome.gff3,$S/lraa_final_abinitio.pred_status.tsv \
  --pred stringtie=$G/stringtie/final_abinitio/transcripts.fa.transdecoder.genome.gff3,$S/stringtie_final_abinitio.pred_status.tsv \
  --pred trinity=$G/trinity/final_abinitio/transcripts.fa.transdecoder.genome.gff3,$S/trinity_final_abinitio.pred_status.tsv \
  -o work/gm12878_td/ribo_gm12878/gm12878
# expression-matched comparison (StringTie transcript TPM)
awk -f $T/analyses/riboseq/expr_matched_support.awk data/gm12878/star_stringtie/stringtie.gtf \
  work/gm12878_td/ribo_gm12878/gm12878.stringtie.orf_support.tsv | sort -k1,1 -k2,2
```
Method summary:
- **Calibration:** on the first 4,000 MANE Select CDSs, choose per read length the 5′ offset (0–20) that maximises the frame-0 fraction.
  - A length is used only if its frame-0 fraction is ≥ 0.5 and it has ≥ 1,000 P-sites.
  - Reads: MAPQ ≥ 10, primary alignments only.
- **Per-ORF scoring:**
  - Count P-sites by frame, excluding the first and last codons.
  - Supported = ≥10 P-sites, frame 0 the largest, and binomial P(X ≥ k | n, 1/3) < 0.01.
  - Only the distinct representative of each CDS is scored.

## Outputs
- `$T/genome_annotation/results/gm12878_riboseq_calibration.tsv`
- `$T/genome_annotation/results/gm12878_riboseq_support.tsv`, by assembly × class × short (<100 aa) / long
- `$T/genome_annotation/results/gm12878_riboseq_stringtie_by_expression.tsv`
- Per-ORF tables (not in git): `$P/work/gm12878_td/ribo_gm12878/gm12878.<asm>.orf_support.tsv`

## Results
| | Rep1 SRR1803147 | Rep2 SRR1803148 |
|---|---|---|
| Reads | 159.4M | 151.3M |
| Passing trim, 25–34 nt | 9.64M (6.0%) | 94.0M (62.2%) |
| Non-rRNA/tRNA | 3.56M | 45.5M |
| Uniquely mapped | 1.71M | 9.75M |

Merged: 11.46M footprints. On MANE CDSs, 98.8% of footprints are sense (1,464,108 sense vs 17,178 antisense).

Calibration: usable lengths and offsets are 27 nt (offset 2, frame-0 0.53), 29 nt (17, 0.555) and 30 nt (2, 0.553). Other lengths are ≤ 0.50.

| Assembly | Class | n | ≥10 P-sites | Supported | Supported / covered |
|---|---|---|---|---|---|
| LRAA | exact | 12,021 | 10,721 (89%) | 7,065 (58.8%) | 0.659 |
| LRAA | novel | 885 | 38 (4.3%) | 13 (1.5%) | 0.342 |
| StringTie | exact | 13,807 | 11,511 (83%) | 7,396 (53.6%) | 0.643 |
| StringTie | novel | 4,127 | 43 (1.0%) | 18 (0.4%) | 0.419 |
| Trinity | exact | 13,832 | 11,639 (84%) | 7,398 (53.5%) | 0.636 |
| Trinity | novel | 11,235 | 70 (0.6%) | 22 (0.2%) | 0.314 |

The compatible and other_overlap classes, and the short/long split, are in `gm12878_riboseq_support.tsv`.

Expression-matched (StringTie), fraction with ≥10 P-sites:

| TPM | Exact | Novel |
|---|---|---|
| <1 | 0.697 | 0.006 |
| 1–5 | 0.767 | 0.011 |
| 5–20 | 0.893 | 0.034 (n=235) |
| 20–100 | 0.963 | 0.042 (n=72) |
| ≥100 | 0.964 | 0.100 (n=30) |

## Checks
- Strand: 98.8% of footprints on MANE CDSs are sense.
- MANE-exact predictions act as the positive control for the frame test.
- The ≥1 P-site fractions follow the same pattern as ≥10 (Trinity novel 6.9% vs exact 97.1%).

## Caveats
- **Rep1 contributes little** (1.7M of 11.5M). Most of its reads fall outside 25–34 nt after trimming: 51.7% too short, 41.9% too long, i.e. no adapter within 36 nt.
- **Weak periodicity from the A+T1 digestion.** The frame test has limited power, so "covered" is the more informative measure.
- **Unique alignments only**, so ORFs in repeats and recent duplicates cannot be assessed.
- **Single-exon transcripts** assembled on the wrong strand would lack sense footprints.
- **Offsets** are chosen by frame-0 fraction. Offsets that are equal mod 3 give near-identical frame assignments, so the absolute offset values (2, 17) are not biologically interpreted.
