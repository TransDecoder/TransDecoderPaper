# E03. GM12878 three-assembly annotation vs MANE/GENCODE (final defaults)

## Claims supported
Manuscript: Results, application to GM12878 (genome-wide CDS annotation of assembled transcripts).

1. PyTD with the final defaults was applied to three GM12878 assemblies: LRAA (PacBio long reads), StringTie (STAR short reads) and Trinity genome-guided (same STAR alignments).
2. **Given a correctly assembled transcript, the exact annotated CDS is recovered 91–93% of the time for MANE and 88–89% for GENCODE**, whatever the assembler. Overall sensitivity is limited by assembly, not by CDS prediction.
3. Trinity GG has the highest GENCODE compatible sensitivity (0.446), but many of its transcripts are fragments (20,200 partial distinct CDSs), so its exact-match share is the lowest.
4. "MANE precision" is low by construction (see Metrics). The GENCODE exact/compatible shares are the fair precision measures.
5. Unannotated predictions (no same-strand overlap with any GENCODE CDS) are 3.7% (LRAA), 12.0% (StringTie) and 15.5% (Trinity) of distinct CDSs.

## Data
- **Short reads:** ENCODE ENCSR000AED rep 1 (CSHL, polyA+, stranded dUTP/RF, HiSeq 2000 PE 101 nt; GEO GSM2072351), ENCFF001REK (R1) + ENCFF001REJ (R2), 97,548,052 pairs, md5 verified. Details and the upstream assembly: record E11 and `$P/data/gm12878/README.md`.
- **Long reads:** ENCODE ENCSR838WFC rep 2 (UCI, PacBio Sequel II Iso-Seq; GEO GSM5331162), ENCFF694DIE (FLNC), 2,538,701 reads, md5 38c1eaea4283c60b6d8f0b3025d540cd verified.
- **Genome and annotation:** GENCODE release 50, GRCh38 primary assembly (`$P/data/gm12878/ref/`; md5 da1a1125… genome, 289b91e5… GTF, from GENCODE MD5SUMS).
- **MANE Select v1.5:** `$P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`. 19,229 MANE Select transcripts on the primary assembly are used; 64 MANE transcripts on fix/alt contigs are skipped.
- **Assemblies (inputs to PyTD):**

| Assembly | Transcripts | Structures | FASTA |
|---|---|---|---|
| LRAA de novo | 39,530 | `$P/data/gm12878/longread_lraa/lraa_denovo/LRAA.gtf` | `.../LRAA.transcripts.fa` |
| StringTie de novo (no `-G`) | 83,342 | `$P/data/gm12878/star_stringtie/stringtie.gtf` | `.../stringtie.transcripts.fa` |
| Trinity GG | 428,960 (428,223 mapped, primary) | `$P/data/gm12878/trinity/Trinity.minimap2.bam` | `$P/data/gm12878/trinity/Trinity.fasta` |

## Software
- **PyTransDecoder** PY_dev 48c870c (`$P/work/trees/PY_dev`; 48c870c037b5ef5e09ba91124318274cc7a200e9). Final defaults:
  - LongOrfs `-m 50` (`DEFAULT_MIN_PROTEIN_LENGTH = 50`);
  - short-ORF guard: ORFs < 100 codons need frame contrast ≥ 200 × the median per-base coding contrast of the training ORFs, unless homology-supported (`--short-orf-length 100`, `--short-orf-min-contrast 200`);
  - additive Kozak term in start refinement (`--start-kozak-weight 25`, `--start-contrast-alpha 0.4`).
- Python 3.11.17 (`$P/env/td`); DIAMOND 2.2.8 and HMMER 3.4 (homology run only).
- Assembly tools (`$P/env/rnaseq_asm`): STAR 2.7.11b, StringTie 3.0.3, minimap2 2.31-r1302, samtools 1.24, gffread 0.12.9. LRAA v0.44.3 (Apptainer image `lraa-core_0.44.3.sif`; see `$P/data/gm12878/README.md`). Trinity 2.15.2: record E11.

## Commands
Upstream assembly (scripts in `$T/genome_annotation/scripts/data_prep/`, listed in its README):
- StringTie: `star_index.sh`, then `star_align_stringtie.sh`:
  - `STAR ... --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS`
  - `stringtie GM12878_ENCSR000AED_rep1.STAR.bam --rf -p 4 -l STRG -o stringtie.gtf`
  - `gffread -w stringtie.transcripts.fa -g GRCh38.primary_assembly.genome.fa stringtie.gtf`
- LRAA: `longread_minimap2.sh`, `pull_lraa.sh`, `lraa_denovo.sh`, `lraa_post.sh`:
  - `minimap2 -t 4 -ax splice:hq -uf --secondary=no GRCh38.primary_assembly.genome.fa ENCFF694DIE.fastq.gz`
  - `LRAA --genome GRCh38.primary_assembly.genome.fa --bam ENCFF694DIE.minimap2.bam --HiFi --cpu_budget 4 -o GM12878_ENCSR838WFC_rep2.LRAA` (de novo, no `--gtf`)
- Trinity GG: `trinity_gg.sh` (record E11).

TransDecoder and evaluation (2026-10-20/21), launched with
```bash
bash $P/work/gm12878_td/overnight.sh     # = $T/genome_annotation/scripts/gm12878_overnight.sh
```
which calls `gm12878_final_defaults.sh <name> <fasta> <structures> <abinitio|homology> [threads]` for each assembly (the work-area copies are identical to the repo copies):
```bash
final_defaults.sh lraa      $D/longread_lraa/lraa_denovo/LRAA.transcripts.fa $D/longread_lraa/lraa_denovo/LRAA.gtf abinitio
final_defaults.sh stringtie $D/star_stringtie/stringtie.transcripts.fa       $D/star_stringtie/stringtie.gtf       abinitio
final_defaults.sh trinity   $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam abinitio
final_defaults.sh trinity   $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam homology 16
```
Steps inside `gm12878_final_defaults.sh` (`td` = `python -m pytransdecoder.cli`):
```bash
td longorfs -t transcripts.fa -S
td predict -t transcripts.fa [--retain-blastp-hits blastp.outfmt6 --retain-pfam-hits pfam.domtblout]
td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $structs -t transcripts.fa \
   -o transcripts.fa.transdecoder.genome.gff3
python eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 --structures $structs \
   --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
   --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${name}_final_$mode \
   -o $P/work/gm12878_td/eval/${name}_final_$mode
python $T/analyses/riboseq/match_riboseq.py --pred transcripts.fa.transdecoder.genome.gff3 \
   --phase1 $P/data/riboseq_orfs/phase1/Ribo-seq_ORFs.all.bed \
   --phase2 $P/data/riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.bed > $P/work/gm12878_td/ribo/${name}_final_$mode.ribo.tsv
```
Homology mode adds `diamond blastp --db sprot_minus_9606.dmnd --max-target-seqs 1 --evalue 1e-5 --outfmt 6` and `hmmsearch --cut_ga --noali --domtblout` against Pfam-A (databases: record E06).

## Metrics (`$T/genome_annotation/scripts/eval_genome_cds.py`)
A CDS is its set of genomic segments, including the stop codon (GTF `CDS` + `stop_codon` merged). Several assembled transcripts can give the same CDS, so predictions are collapsed to **distinct predicted CDSs** (identical chromosome, strand and segments; the `complete` type wins).

Match levels between a prediction P and a reference CDS R (same chromosome and strand):
- **exact:** identical segments.
- **compatible:** P is an in-frame, contiguous piece of R. P's internal boundaries are R's intron boundaries and P's codons are R's codons. This credits a correct partial CDS.

Reference sets: **MANE** = 19,229 MANE Select transcripts on the primary assembly; **GENCODE complete** = 253,709 GENCODE v50 `protein_coding` transcripts without `cds_start_NF`/`cds_end_NF`.

Reference-centric (per reference transcript):
- **assembled:** some assembled transcript contains R's whole CDS with all of R's CDS introns. This is the upper bound for any CDS caller on that assembly.
- **sens exact / sens compatible:** fraction of R with an exact / compatible prediction.
- **exact given assembled:** (R exact and assembled) / R assembled. This isolates CDS prediction from assembly completeness.

Prediction-centric (per distinct predicted CDS, split complete / partial / all):
- **MANE precision exact:** among predictions overlapping a MANE CDS on the same strand, the fraction identical to that MANE CDS. **Low by design:**
  - a correct CDS of an alternative isoform (GENCODE-annotated but not MANE Select) counts as a miss;
  - partial predictions can never be identical to a full CDS (pred_partial value is 0.0000 for all assemblies).
- **GENCODE exact / compatible frac:** fraction identical to / compatible with any GENCODE v50 CDS (all transcript types with a CDS, including partial-CDS transcripts). These are the fair precision measures.
- **no GENCODE CDS overlap frac:** fraction with no same-strand overlap with any GENCODE CDS ("novel"/unannotated).
- **pred_complete GENCODE-exact frac:** the same, for complete predictions only.

## Outputs
- `$T/data/gm12878/eval_summaries/{lraa,stringtie,trinity}_final_abinitio.summary.tsv`
- `$T/data/gm12878/eval_summaries/trinity_final_homology.summary.tsv` (copied from `$P/work/gm12878_td/eval/` on 2026-10-21 for this record)
- Not in git:
  - predictions `$P/work/gm12878_td/<asm>/final_<mode>/transcripts.fa.transdecoder.genome.gff3` (and `.bed`, `.pep`, `.cds`);
  - per-reference / per-prediction status `$P/work/gm12878_td/eval/<asm>_final_<mode>.ref_status.tsv`, `$P/work/gm12878_td/eval/status/<asm>_final_abinitio.pred_status.tsv`;
  - Ribo-seq catalogue matches `$P/work/gm12878_td/ribo/<asm>_final_<mode>.ribo.tsv`;
  - logs `$P/work/gm12878_td/logs/final_*.log`, `$P/work/gm12878_td/logs_overnight.out`.

## Results
Ab initio, final defaults (values from the summary tables):

| | LRAA | StringTie | Trinity GG |
|---|---|---|---|
| Assembled transcripts (with genome structure) | 39,530 | 83,342 | 428,223 |
| Transcript ORFs predicted / projected to genome | 34,592 / 34,592 | 40,666 / 40,666 | 97,088 / 92,464 |
| Distinct predicted CDSs | 23,922 | 34,404 | 72,323 |
| complete / partial | 18,807 / 5,115 | 27,034 / 7,370 | 52,123 / 20,200 |
| MANE assembled | 6,294 (0.327) | 8,919 (0.464) | 7,676 (0.399) |
| MANE exact (sens) | 5,882 (0.306) | 8,268 (0.430) | 6,951 (0.362) |
| MANE compatible (sens) | 7,604 (0.395) | 10,938 (0.569) | 10,519 (0.547) |
| MANE exact given assembled | 0.935 | 0.927 | 0.906 |
| GENCODE complete assembled | 54,830 (0.216) | 66,881 (0.264) | 60,502 (0.239) |
| GENCODE complete exact (sens) | 48,765 (0.192) | 59,231 (0.234) | 53,007 (0.209) |
| GENCODE complete compatible (sens) | 68,291 (0.269) | 89,168 (0.352) | 113,048 (0.446) |
| GENCODE exact given assembled | 0.889 | 0.886 | 0.876 |
| All predictions: GENCODE exact / compatible frac | 0.503 / 0.726 | 0.401 / 0.616 | 0.191 / 0.423 |
| Complete predictions: GENCODE exact frac | 0.639 | 0.511 | 0.265 |
| All predictions: no GENCODE CDS overlap | 3.7% | 12.0% | 15.5% |
| All predictions: MANE precision exact | 0.257 | 0.276 | 0.115 |
| Ribo-seq ORF catalogue matches, stop (exact), prediction models | 178 (33) | 217 (61) | 455 (130) |

Trinity GG with homology retention (DIAMOND vs Swiss-Prot minus human + Pfam): 84,745 distinct CDSs; MANE exact 7,063 (0.367), compatible 0.570, exact given assembled 0.920; GENCODE complete exact 0.215, compatible 0.481, exact given assembled 0.902; no GENCODE CDS overlap 20.8%. The homology series is record E06.

Notes on the table:
- **Trinity ORFs not projected:** 4,624 of 97,088 (4.8%). Of these, 4,611 have a CDS not entirely within the aligned part of the transcript (soft-clipped or split alignment), and 13 lie on transcripts with no genome structure (unmapped) (`genome_map.log`).
- **Ribo-seq matches** count prediction models (rows of `*.ribo.tsv`), not distinct CDSs or distinct catalogue ORFs. The matched catalogue ORFs are 120 (24 exact) for LRAA, 152 (46) for StringTie and 258 (69) for Trinity. Catalogues: GENCODE Ribo-seq ORFs phase 1 (Mudge et al. 2022) and phase 2.

## Checks
- MANE-exact predictions are the positive control for Ribo-seq (E01) and MS (E02) support.
- The assembled-transcript counts in the evaluation equal the assembly sizes: LRAA GTF 39,530; StringTie GTF 83,342; Trinity minimap2 primary mapped 428,223 of 428,960.
- `predicted_cds_models` equals the number of ORFs projected by genome-map for each assembly (34,592 / 40,666 / 92,464).
- The three ab initio summary tables in git are byte-identical to `$P/work/gm12878_td/eval/`.
- Every reference counted exact is also counted assembled (0 exceptions in all three `ref_status.tsv` files), as expected.

## Caveats
- **One sample per technology:** one short-read replicate and one long-read replicate. LRAA sensitivity is limited by long-read depth (2.5M FLNC reads).
- **"Assembled" needs every CDS intron in a single transcript.** A CDS split across Trinity fragments counts as not assembled, though compatible partial predictions are still credited.
- **GENCODE exact/compatible fractions count matches to any GENCODE CDS.** This includes NMD and partial-CDS transcripts, so they are not restricted to the GENCODE-complete reference set used for sensitivity.
- **The overnight run called `$P/work/gm12878_td/final_defaults.sh`.** That file is identical to `$T/genome_annotation/scripts/gm12878_final_defaults.sh`. The git commit of this repo at run time was not recorded.
