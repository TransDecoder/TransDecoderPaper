# E06. Homology restriction series (clade-excluded Swiss-Prot), GM12878

## Claims supported
Manuscript: Results, homology integration and the GM12878 application (to be drafted).

1. **Homology evidence gives modest gains** over ab initio predictions on assembled transcripts.
   - MANE exact given assembled, Swiss-Prot minus human: +1.2 (LRAA) / +1.1 (StringTie) / +1.4 (Trinity) points.
   - MANE compatible: +0.9 / +2.3 / +2.1 points.
2. **The Swiss-Prot gain shrinks with evolutionary distance** (human → primates → mammals → vertebrates → all animals removed). Excluding all animals still gives a small gain.
3. **Pfam alone recovers most of the gain, independent of the BLAST database.**
   - Pfam HMMs are built from all organisms, including human, so this is not distance-controlled.
4. **Retention artifact from an Alu-derived segment.**
   - One curated protein with an Alu-derived segment drives most homology-retained unannotated ORFs: chimp UTY, Swiss-Prot Q6B4Z3 (UTY_PANTR).
   - Its C-terminal ~85 aa (residues ~995–1079) is Alu-encoded (…QWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQAHLELLTSGG…).
   - Minus-human vs minus-primates (BLAST only): 4,692 of the 4,981 Trinity CDSs found only in the minus-human run are unannotated. Of those, 4,468 have UTY_PANTR as best hit, 95% are <100 aa, and 185 more hit slow-loris LINE-1 ORF1 (LIN1_NYCCO).
   - StringTie: 1,232 of 1,305.
   - GM12878 is female, so UTY is not expressed. These are Alu-containing transcripts.

## Data
- Swiss-Prot 2026_03 (`$P/data/homology/uniprot_sprot.fasta.gz`).
- NCBI taxonomy dump (`$P/data/taxonomy/nodes.dmp`, downloaded 2026-10-21).
- Pfam-A 38.2 (30,134 families).
- Clade-excluded databases:

| Clade removed | Taxid | Entries removed | Entries kept |
|---|---|---|---|
| Human | 9606 | 20,432 | 555,316 (existing `sprot_minus_9606.dmnd`: 555,317) |
| Primates | 9443 | 27,203 | 548,545 |
| Mammals | 40674 | 67,970 | 507,778 |
| Vertebrates | 7742 | 87,072 | 488,676 |
| Metazoa | 33208 | 110,685 | 465,063 |

  405 entries have taxids missing from nodes.dmp and are kept.

## Software
- PyTD PY_dev 48c870c.
- DIAMOND (env/td): `blastp --max-target-seqs 1 --evalue 1e-5 --outfmt 6`, default sensitivity.
- HMMER `hmmsearch --cut_ga`: chunked single-threaded for LRAA and StringTie; `--cpu 16` for Trinity, run by `gm12878_final_defaults.sh`.

## Commands
```bash
python3 -I $T/genome_annotation/scripts/make_sprot_clade_excluded.py $P/data/homology/uniprot_sprot.fasta.gz \
  $P/data/taxonomy/nodes.dmp <taxid> > sprot_minus_clade_<taxid>.fasta
diamond makedb --in sprot_minus_clade_<taxid>.fasta -d sprot_minus_clade_<taxid>
bash $T/genome_annotation/scripts/homology_series_driver.sh   # 33 configs; calls homology_series.sh
python3 $T/analyses/homology_series_retention.py <A.pred_status> <B.pred_status> <A blastp outfmt6>
```
- The series reuses each assembly's final_abinitio LongOrfs output, then runs predict, genome-map and `eval_genome_cds.py`.
- Configurations per assembly:
  - BLAST only × 5 clades;
  - Pfam only;
  - BLAST + Pfam × 5 clades.

## Outputs
- `$T/genome_annotation/results/gm12878_homology_series.tsv` (summary table)
- `$T/genome_annotation/results/eval_summaries/hseries/*.summary.tsv`
- `$T/genome_annotation/results/gm12878_homology_primate_retention.txt`
- Not in git: `$P/work/gm12878_td/<asm>/hseries/`, `$P/work/gm12878_td/eval/hseries/*status.tsv`

## Results (Trinity; LRAA and StringTie are in the table)
| Config | Distinct CDS | MANE exact | MANE compatible | Exact given assembled | GENCODE compatible | Complete, no GENCODE overlap |
|---|---|---|---|---|---|---|
| Ab initio | 72,323 | 0.3615 | 0.5470 | 0.9055 | 0.4456 | 12.2% |
| BLAST minus human | 84,216 | 0.3671 | 0.5680 | 0.9196 | 0.4795 | 15.5% |
| BLAST minus primates | 79,545 | 0.3671 | 0.5677 | 0.9196 | 0.4793 | 12.4% |
| BLAST minus mammals | 75,096 | 0.3648 | 0.5569 | 0.9138 | 0.4623 | 11.9% |
| BLAST minus vertebrates | 73,802 | 0.3635 | 0.5522 | 0.9105 | 0.4552 | 12.0% |
| BLAST minus metazoa | 73,306 | 0.3630 | 0.5506 | 0.9095 | 0.4519 | 12.1% |
| Pfam only | 76,151 | 0.3673 | 0.5626 | 0.9200 | 0.4647 | 12.4% |
| BLAST minus human + Pfam | 84,745 | 0.3673 | 0.5704 | 0.9201 | 0.4808 | 15.6% |
| BLAST minus metazoa + Pfam | 76,373 | 0.3669 | 0.5631 | 0.9192 | 0.4659 | 12.4% |

## Checks
- **Human exclusion by taxonomy** matches the earlier OX=9606 database to within one entry.
- **Minus-human vs minus-primates** differ almost only in unannotated short ORFs. The sensitivity columns are unchanged (0.3671 / 0.5680 vs 0.3671 / 0.5677). This localises the artifact.

## Caveats and follow-ups
- **No retention filter.** TransDecoder retains any ORF with a hit at E ≤ 1e-5, regardless of hit coverage or repeat origin. Options to evaluate:
  - a subject-coverage minimum;
  - masking TE/Alu-derived segments in the database;
  - reporting TE-homology-only ORFs separately, as candidate neo-ORFs (Brian 2026-10-21).
- **Pfam is not distance-controlled**, because its HMMs include human sequences.
- **Default DIAMOND sensitivity** may under-detect distant homologs; `--more-sensitive` has not been tested.
