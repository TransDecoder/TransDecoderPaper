# E02. GM12878 MS peptide support of predicted ORFs

## Claims supported
Manuscript: Results, "Translation and protein evidence for predicted ORFs in GM12878"; Methods, "Mass-spectrometry support". Also Discussion (notes/05).

1. Peptide coverage of predictions that exactly match GENCODE:
   - ≥100 aa: 92–94% have ≥1 peptide, nearly all of them ≥2.
   - <100 aa: 81–87%.
2. 10–25% of unannotated predictions share a peptide with the identified set. These are almost entirely paralogous or pseudogene copies of known genes (e.g. NCF1B/C, HERC2P3, processed RPS/RPL pseudogenes), where peptides cannot distinguish copy from parent.
3. Non-reference peptides on unannotated predictions are nearly all GM12878 variant-allele peptides of parent genes, landing on duplicated copies (e.g. a PWP2 copy on 21p; NAPSB).
4. One credible exception: C3orf86P, with 3 independent peptides from 3 proteases. It was protein-coding in the Ensembl release used for the original search and is a transcribed unitary pseudogene in GENCODE v50.
5. Because the original search database did not contain our predictions, this analysis confirms known proteins but cannot discover new ones.

## Data
- **Proteome: Sinitcyn et al. 2023, Nat Biotechnol 41:1776 (doi:10.1038/s41587-023-01714-x).** MassIVE MSV000086944 / PXD024364. PRIDE has no record of PXD024364, so the data came from MassIVE over HTTPS, downloaded 2026-10-21.
- **Search output:** `search/GM12878.zip`, 5.52 GB, md5 f062982347bb27a12b42c52823275b08. This is the MaxQuant 1.6.1.0 `txt/` folder.
  - Only `GM12878/peptides.txt` was extracted (md5 8bf9a4449708b9ed79fb72e741e2f026), plus `parameters.txt`, `summary.txt` and `tables.pdf`.
  - Search settings: 6 proteases (AspN, chymotrypsin, GluC, LysC, LysN, trypsin), PSM and protein FDR 1%, minimum peptide length 7.
  - Search database: `cell_lines.mutations.fasta` (Ensembl ENSP plus GM12878 variant entries plus contaminants), in `sequence/`.
- **Reference proteomes** (`refdb/`):
  - GENCODE v50 `gencode.v50.pc_translations.fa.gz`: 382,428 sequences.
  - UniProt UP000005640 with isoforms: 169,651 sequences, including TrEMBL.
- **Local root:** `$P/data/gm12878/proteome_sinitcyn2023/`.
- **Predictions:** `$P/work/gm12878_td/<asm>/final_abinitio/transcripts.fa.transdecoder.pep`. Classes are from pred_status (E01, E03).

## Software
Python 3.11 with pyahocorasick and pandas (`$P/env/proteo`). Script: `$T/analyses/proteomics/ms_peptide_support.py`.

## Commands
```bash
cd $P; D=data/gm12878/proteome_sinitcyn2023; S=work/gm12878_td/eval/status; W=work/gm12878_td
env/proteo/bin/python -I $T/analyses/proteomics/ms_peptide_support.py \
  --peptides $D/search/GM12878/peptides.txt \
  --ref $D/refdb/gencode.v50.pc_translations.fa.gz --ref $D/refdb/uniprot_human_UP000005640_with_isoforms.fasta.gz \
  --pred lraa=$W/lraa/final_abinitio/transcripts.fa.transdecoder.pep,$S/lraa_final_abinitio.pred_status.tsv \
  --pred stringtie=$W/stringtie/final_abinitio/transcripts.fa.transdecoder.pep,$S/stringtie_final_abinitio.pred_status.tsv \
  --pred trinity=$W/trinity/final_abinitio/transcripts.fa.transdecoder.pep,$S/trinity_final_abinitio.pred_status.tsv \
  -o $W/proteomics
```
Method:
- Drop rows with Reverse or Potential contaminant set to '+'.
- Normalise I→L and deduplicate. Sum MS/MS counts, and keep the best score and lowest PEP.
- Exact substring matching (Aho–Corasick) against the predicted proteins and the two reference sets.
- Non-reference peptide = absent from both reference sets.

Locus classes of novel ORFs with peptides, relative to the parent gene of their best reference protein, plus overlap with GENCODE pseudogenes:
```bash
python3 -I $T/analyses/proteomics/novel_peptide_loci.py $W/proteomics/ms_peptide_support.novel_with_peptides.tsv \
  $P/data/gm12878/ref/gencode.v50.primary_assembly.annotation.gtf \
  lraa=$W/lraa/final_abinitio/transcripts.fa.transdecoder.genome.gff3 \
  stringtie=$W/stringtie/final_abinitio/transcripts.fa.transdecoder.genome.gff3 \
  trinity=$W/trinity/final_abinitio/transcripts.fa.transdecoder.genome.gff3 -o $W/proteomics/novel_peptide_loci
```

## Outputs
- Committed: `$T/data/gm12878/gm12878_ms_peptide_support.tsv`, `gm12878_ms_peptide_counts.tsv`, `gm12878_ms_novel_with_peptides.tsv`.
- Not in git: `$P/work/gm12878_td/proteomics/ms_peptide_support.{<asm>.model_hits,novel_specific_hits}.tsv`.

## Results
| Peptides | Count |
|---|---|
| `peptides.txt` rows | 768,554 |
| Reverse removed | 1,747 |
| Contaminants removed | 8,204 |
| Kept | 758,603 |
| Distinct after I/L normalisation | 758,225 |
| In GENCODE v50 or UniProt | 753,266 |
| Absent from both | 4,959 (4,913 flagged "Mutated" by MaxQuant) |

Predictions with ≥1 / ≥2 distinct peptides:

| Assembly | Class | Short (<100 aa) | Long |
|---|---|---|---|
| LRAA | exact | 281 / 269 of 323 | 11,041 / 10,892 of 11,698 |
| LRAA | novel | 26 / 23 of 150 | 182 / 147 of 735 |
| StringTie | exact | 305 / 284 of 378 | 12,546 / 12,347 of 13,429 |
| StringTie | novel | 179 / 149 of 874 | 670 / 547 of 3,253 |
| Trinity | exact | 485 / 455 of 565 | 12,479 / 12,309 of 13,267 |
| Trinity | novel | 234 / 155 of 3,396 | 784 / 585 of 7,839 |

Full table, including compatible and other_overlap: `gm12878_ms_peptide_support.tsv`.

Novel ORFs sharing peptides with known proteins:

| | LRAA | StringTie | Trinity |
|---|---|---|---|
| Novel ORFs with peptides | 208 | 849 | 1,018 |
| ...≥2 peptides from one reference protein | 169 | 691 | 732 |
| ...contained in a reference protein | 29 | 144 | 108 |
| ...overlapping a GENCODE pseudogene | 172 (83%) | 759 (89%) | 777 (76%) |

Location relative to the parent gene (`novel_peptide_loci.py`; `gm12878_ms_novel_peptide_loci.summary.tsv`):

| Class | LRAA | StringTie | Trinity |
|---|---|---|---|
| Overlaps the parent gene span (nested copy; all also overlap a pseudogene) | 43 | 102 | 104 |
| Same chromosome, ≤2 Mb (tandem/segmental duplication) | 39 | 131 | 176 |
| Same chromosome, >2 Mb | 53 | 151 | 264 |
| Other chromosome (e.g. processed pseudogenes) | 73 | 462 | 473 |
| Parent gene name not in GENCODE | 0 | 3 | 1 |

These agree with the one-off check made during the analysis to within ±2 per class; the differences come from the choice of nearest parent locus.

Novel-specific peptides: 37 novel-class models (3 / 12 / 22) at 11 loci. Nearly all are variant-allele peptides of a parent gene on a duplicated copy. Examples:
- **PWP2 D25N** on a chr21:5.13 Mb copy. GNLNFTCNGNSVLSPVGNR, 5 PSMs, PEP 2.7e-37.
- **NAPSA I40T** on NAPSB (chr19:50.33 Mb).
- **AMZ2 N30D** at chr17:64.97 Mb, consistent with AMZ2P1.
- **DTX2 G384E** on a chr7:77 Mb copy.
- **EIF4B L471Q** on a chr7:104.67 Mb copy.

Non-variant, non-reference peptides:
- **C3orf86P** (chr3:44.42 Mb; ENSG00000225873). StringTie STRG.10275.2.p1 and Trinity TRINITY_GG_18087_c0_g1_i3.p1 (101 aa) carry 3 peptides:
  - SRGQFGQGQE (GluC, PEP 0.024);
  - KADAPAVSPESPQ (LysN, 2 PSMs, PEP 0.0011);
  - DSGSRNPCPSAPTWARPKPEE (AspN, PEP 1.5e-4).
- **USP41P** (chr22:20.36 Mb): ambiguous. Effectively one region, and the peptide also matches a USP18 variant entry.
- **STRG.14334.1.p1:** a 7-mer the authors assigned to ACD, likely a chance match.

## Checks
- **Reproduction (2026-10-21):** rerunning the command above into `$P/work/gm12878_td/proteomics_repro` gave an identical summary table (175 s, 2.2 GB RAM).
- The exact-match classes act as the positive control.
- Trinity's exact, compatible and other_overlap classes contain non-reference (variant) peptides. LRAA and StringTie do not, because Trinity assembles GM12878 alleles from reads while the others take sequence from the reference genome. This is a consistency check.

## Caveats
- The published search contained only Ensembl and variant proteins, so this cannot detect new proteins. That would need a re-search of the raw data (534 GB).
- "≥2 distinct peptides" counts nested or overlapping peptides from different proteases separately.
- Tryptic status means little for the other five proteases. The experiments (protease and fragmentation) are listed per peptide.

## Release (v7) update (2026-10-25)
- **Claim 2 correction (development runs):** the share of unannotated predictions with a peptide was
  9-24%, not 10-25% (`notebooks/07_gm12878_evidence.Rmd`, earlier-run checks).
- **Release runs** (PyTD 3a793ff, release defaults; E25; tables in `data/benchmark/final_v7/` and
  `data/gm12878/final_v7/`): exact-match predictions >= 100 aa with >= 1 peptide 94-95% (development
  93-94%); short exact-match (< 100 aa) 85-90% (81-87%); unannotated with any peptide 30-49%
  (167 of 341, 748 of 1,687, 765 of 2,528 for LRAA, StringTie, Trinity; development 9-24% of
  885 / 4,127 / 11,235). The rise comes from the smaller denominator: the q-value test removed many
  short unannotated ORFs without peptides. Unannotated CDSs with peptides overlapping a GENCODE
  pseudogene: 86-93% (development 76-89%).
- **C3orf86P:** the 101-aa ORF is predicted only in the StringTie `--sensitive` run, not in any release
  default run. The peptide evidence above is unchanged.
