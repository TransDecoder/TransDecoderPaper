# E06. Homology restriction series (clade-excluded Swiss-Prot), GM12878

## Claims supported
Manuscript: Results, "Homology evidence and evolutionary distance" and "Homology to transposable elements" (results_gm12878_parameters.tex), Fig. 7c,d, Supplementary Table S7. Since 2026-10-30 the manuscript uses the release rerun (section "Release rerun with PyTD a5b24c3" at the end); the numbers in the claims below are the development series.

1. **Homology evidence gives modest gains** over ab initio predictions on assembled transcripts.
   - MANE exact given assembled, Swiss-Prot minus human: +1.8 (LRAA) / +1.1 (StringTie) / +1.4 (Trinity) points (0.9345 → 0.9528, 0.9270 → 0.9380, 0.9055 → 0.9196; ab initio = `<asm>_final_abinitio.summary.tsv`, minus human = `hseries/<asm>_blast_9606_nopfam.summary.tsv`). Corrected 2026-10-22: LRAA was given as +1.2.
   - MANE compatible: +0.9 / +2.3 / +2.1 points.
2. **The Swiss-Prot gain shrinks with evolutionary distance** (human → primates → mammals → vertebrates → all animals removed). Excluding all animals still gives a small gain.
3. **Pfam alone recovers most of the gain, independent of the BLAST database.**
   - Pfam HMMs are built from all organisms, including human, so this is not distance-controlled.
4. **Retention artifact from an Alu-derived segment.**
   - One curated protein with an Alu-derived segment drives most homology-retained unannotated ORFs: chimp UTY, Swiss-Prot Q6B4Z3 (UTY_PANTR).
   - Its C-terminal ~85 aa (residues ~995–1079) is Alu-encoded (…QWCDLSSLQPPPPGFKRFSHLSLPNSWNYRHLPSCPTNFCIFVETGFHHVGQAHLELLTSGG…).
   - Minus-human vs minus-primates (BLAST only): 4,692 of the 4,981 Trinity CDSs found only in the minus-human run are unannotated. Of those, 4,465 have UTY_PANTR as best hit (4,468 have a chimpanzee protein, PANTR), 95% are <100 aa, and 185 more hit slow-loris LINE-1 ORF1 (LIN1_NYCCO).
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
- `$T/data/gm12878/gm12878_homology_series.tsv` (summary table)
- `$T/data/gm12878/eval_summaries/hseries/*.summary.tsv`
- `$T/data/gm12878/gm12878_homology_primate_retention.txt`
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

## Decision (Brian, 2026-10-21)
No change to TransDecoder's homology retention.
- TE-homology-retained ORFs are kept and presented as candidate neo-ORFs. TE-derived peptides are of interest as tumor antigens for immunotherapy.
- Users can screen them out easily, by best-hit protein or by repeat overlap.
- Paper framing: these ORFs are retained on homology when their coding signal alone would not retain them. Most are short (95% <100 aa).

## Caveats
- Pfam is not distance-controlled, because its HMMs include human sequences.
- Default DIAMOND sensitivity may under-detect distant homologs; `--more-sensitive` has not been tested.

## Release rerun with PyTD a5b24c3 (2026-10-30)
Approved by Brian (2026-10-30): redo the series with the release code and use it for Fig. 7c,d, the
Results paragraphs and a supplementary table. The development series above (48c870c, fixed
short-ORF test, no start assist) and its files are kept for comparison.

### Code and runs
- **PyTD a5b24c3** (`$P/work/trees/PY_hstart`; `pytd_version.txt` checked in every run directory by
  the scripts): release defaults (decoy q-value test q 0.05, 50-codon minimum, start refinement), and
  with hits the homology-assisted start refinement (on by default; E33 "Adoption in PyTD"). Same
  invocation as E33's GM12878 Trinity run (`13_gm12878_trinity_v7h_v1.sh`): fresh `longorfs -S`,
  `predict --retain-blastp-hits ... --retain-pfam-hits ...`, `genome-map`, `eval_genome_cds.py`
  (MANE v1.5, GENCODE v50).
- **Scripts** (`genome_annotation/scripts/`): `homology_series_v7h_prep.sh` (per assembly: one fresh
  LongOrfs, DIAMOND per clade, Pfam file), `homology_series_v7h.sh` (one configuration: fresh run
  directory, copy of the fresh LongOrfs work directory, Predict, genome-map, eval),
  `homology_series_v7h_driver_v1.sh` and `_v2.sh` (3 configurations at a time).
  - v1 ran niced; after Brian asked (2026-10-30) for normal priority, v1 was stopped by PID after
    prep and 9 Trinity configurations; the 3 Trinity runs in progress (~6 min old) were moved to
    `$P/_trash/hseries_v7h_niced_killed/` and rerun fresh by v2 (no nice), which ran the remaining
    27 configurations. Logs `$P/work/gm12878_td/logs/hseries_v7h/{driver_v1,driver_v2,prep.*,*}.log`:
    36 configurations, 0 failed. (An earlier v1 start, killed within a minute because its log went to
    /dev/null, is in `$P/_trash/hseries_v7h_aborted_start/`.)
- **Outputs (not in git):** `$P/work/gm12878_td/<asm>/hseries_v7h/{_longorfs,<config>}`,
  `blastp.minus_<clade>.outfmt6`, `pfam.domtblout`; `$P/work/gm12878_td/eval/hseries_v7h/`.
- **12 configurations per assembly:** ab initio (no hits), BLAST × 5 clades, Pfam only, BLAST + Pfam ×
  5 clades; databases as above.

### Search inputs and their verification
- **LongOrfs:** the fresh a5b24c3 `longest_orfs.pep` equals (`cmp`) the release ab initio run's
  (`final_v7_abinitio`, 3a793ff) for all three assemblies, and also equals the md5 of the
  `final_abinitio` LongOrfs used by the development series (LRAA 9b1baeee..., StringTie fabc3bc3...,
  Trinity 85695c77...; 117,560 / 233,507 / 633,606 peptides).
- **DIAMOND:** run fresh (blastp `--max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads 4`, default
  sensitivity) for all 15 assembly × clade pairs. All 15 outputs are byte-identical to the
  development series' `hseries/blastp.minus_<clade>.outfmt6`, and Trinity minus human is identical to
  `trinity/final_v7_homology/blastp.outfmt6`.
- **Pfam (reused, verified):** hmmsearch `--cut_ga` outputs on the identical peptide files were copied
  (`pfam.domtblout.source`): Trinity from `trinity/final_v7_homology` (`--cpu 8`, 2026-10-25; its
  hits are identical, ignoring comment lines, to the independent 2026-10-21 search in
  `final_homology`, 214,744 domain lines); LRAA and StringTie from the development series'
  12-chunk searches (each chunk ends `# [ok]`). Spot check: a fresh single-threaded `hmmsearch
  --cut_ga` of 2,000 random peptides per assembly (`random.Random(6)`) gave the same domain hits as
  the reused files on every field except E-values (1,526 LRAA and 835 StringTie domain lines;
  target, model, scores, domain numbers, HMM, alignment and envelope coordinates identical).
  E-values differ because they scale with the number of sequences searched (chunked); PyTD uses
  neither E-value for selection or start decisions (hits pass the gathering threshold; R3 uses the
  envelope start), so the predictions do not depend on it.

### Checks
- **Ab initio** a5b24c3 runs are byte-identical to the E25 release runs (final GFF3, pep, genome GFF3;
  all three assemblies), so the ab initio row equals E25's.
- **Trinity Swiss-Prot minus human + Pfam** is byte-identical (GFF3, pep, genome GFF3; eval summary)
  to E33's `trinity/final_v7h_homology` run.
- Notebook 08: all homology and TE checks pass (new expected values below; the development values are
  kept in comments). The other `ok = FALSE` rows in its check table (q-value, minimum-length
  paragraphs) predate this work and are unchanged.

### Scoring and data (as above, with `_v7h` names)
- `pipeline/derive/08_homology_series_v7h.sh` -> `data/gm12878/gm12878_homology_series_v7h.tsv`
  (same columns as `gm12878_homology_series.tsv`), `data/gm12878/eval_summaries/hseries_v7h/`,
  `data/gm12878/gm12878_homology_primate_retention_v7h.txt` (`homology_series_retention.py`, minus
  human vs minus primates, BLAST only; now also LRAA).
- `pipeline/derive/08_homology_series_v7h_table.py` -> `data/gm12878/gm12878_homology_series_v7h.tex`
  (rows of Supplementary Table S7, `tab:hseries`).
- Fig. 7c,d: notebook 08 (`figures/fig7_parameters.{pdf,png}`).

### Results: development (48c870c) -> release (a5b24c3)

**lraa** (development 48c870c -> release a5b24c3)

| config | distinct CDS | MANE exact given assembled | gain (points) | MANE compatible | complete, no GENCODE overlap |
|---|---|---|---|---|---|
| abinitio | 23,922 -> 22,769 | 0.9345 -> 0.9291 |  | 0.3954 -> 0.3918 | 3.9% -> 1.4% |
| blast_9606_nopfam | 24,854 -> 24,300 | 0.9528 -> 0.9728 | +1.8 -> +4.4 | 0.4039 -> 0.4067 | 4.1% -> 2.4% |
| blast_9443_nopfam | 24,777 -> 24,203 | 0.9527 -> 0.9712 | +1.8 -> +4.2 | 0.4037 -> 0.4060 | 3.8% -> 1.9% |
| blast_40674_nopfam | 24,356 -> 23,590 | 0.9446 -> 0.9474 | +1.0 -> +1.8 | 0.3999 -> 0.3983 | 3.8% -> 1.6% |
| blast_7742_nopfam | 24,161 -> 23,253 | 0.9411 -> 0.9382 | +0.7 -> +0.9 | 0.3979 -> 0.3948 | 3.8% -> 1.5% |
| blast_33208_nopfam | 24,095 -> 23,129 | 0.9382 -> 0.9357 | +0.4 -> +0.7 | 0.3967 -> 0.3937 | 3.9% -> 1.5% |
| pfam | 24,478 -> 23,959 | 0.9533 -> 0.9480 | +1.9 -> +1.9 | 0.4025 -> 0.3993 | 3.8% -> 1.9% |
| blast_9606_pfam | 24,898 -> 24,397 | 0.9538 -> 0.9730 | +1.9 -> +4.4 | 0.4040 -> 0.4067 | 4.1% -> 2.5% |
| blast_9443_pfam | 24,824 -> 24,306 | 0.9536 -> 0.9712 | +1.9 -> +4.2 | 0.4038 -> 0.4059 | 3.8% -> 2.1% |
| blast_40674_pfam | 24,615 -> 24,099 | 0.9536 -> 0.9549 | +1.9 -> +2.6 | 0.4032 -> 0.4013 | 3.8% -> 1.9% |
| blast_7742_pfam | 24,534 -> 24,026 | 0.9534 -> 0.9496 | +1.9 -> +2.0 | 0.4025 -> 0.3996 | 3.8% -> 1.9% |
| blast_33208_pfam | 24,520 -> 24,002 | 0.9533 -> 0.9488 | +1.9 -> +2.0 | 0.4025 -> 0.3994 | 3.9% -> 1.9% |

**stringtie** (development 48c870c -> release a5b24c3)

| config | distinct CDS | MANE exact given assembled | gain (points) | MANE compatible | complete, no GENCODE overlap |
|---|---|---|---|---|---|
| abinitio | 34,404 -> 30,846 | 0.9270 -> 0.9149 |  | 0.5688 -> 0.5630 | 10.3% -> 3.8% |
| blast_9606_nopfam | 38,012 -> 36,088 | 0.9380 -> 0.9599 | +1.1 -> +4.5 | 0.5919 -> 0.5980 | 12.3% -> 7.3% |
| blast_9443_nopfam | 36,877 -> 34,660 | 0.9384 -> 0.9595 | +1.1 -> +4.5 | 0.5916 -> 0.5974 | 10.5% -> 5.0% |
| blast_40674_nopfam | 35,380 -> 32,520 | 0.9335 -> 0.9371 | +0.6 -> +2.2 | 0.5801 -> 0.5800 | 10.2% -> 4.2% |
| blast_7742_nopfam | 34,895 -> 31,763 | 0.9316 -> 0.9273 | +0.5 -> +1.2 | 0.5752 -> 0.5722 | 10.3% -> 4.1% |
| blast_33208_nopfam | 34,707 -> 31,438 | 0.9299 -> 0.9221 | +0.3 -> +0.7 | 0.5724 -> 0.5682 | 10.4% -> 4.0% |
| pfam | 35,801 -> 33,328 | 0.9392 -> 0.9334 | +1.2 -> +1.8 | 0.5856 -> 0.5819 | 10.4% -> 4.7% |
| blast_9606_pfam | 38,169 -> 36,335 | 0.9396 -> 0.9613 | +1.3 -> +4.6 | 0.5938 -> 0.5996 | 12.3% -> 7.4% |
| blast_9443_pfam | 37,047 -> 34,928 | 0.9396 -> 0.9604 | +1.3 -> +4.5 | 0.5934 -> 0.5988 | 10.6% -> 5.2% |
| blast_40674_pfam | 36,126 -> 33,770 | 0.9388 -> 0.9428 | +1.2 -> +2.8 | 0.5876 -> 0.5874 | 10.4% -> 4.8% |
| blast_7742_pfam | 35,905 -> 33,495 | 0.9395 -> 0.9370 | +1.2 -> +2.2 | 0.5865 -> 0.5844 | 10.5% -> 4.8% |
| blast_33208_pfam | 35,857 -> 33,414 | 0.9389 -> 0.9340 | +1.2 -> +1.9 | 0.5859 -> 0.5829 | 10.5% -> 4.8% |

**trinity** (development 48c870c -> release a5b24c3)

| config | distinct CDS | MANE exact given assembled | gain (points) | MANE compatible | complete, no GENCODE overlap |
|---|---|---|---|---|---|
| abinitio | 72,323 -> 59,514 | 0.9055 -> 0.8933 |  | 0.5470 -> 0.5336 | 12.2% -> 3.1% |
| blast_9606_nopfam | 84,216 -> 75,523 | 0.9196 -> 0.9363 | +1.4 -> +4.3 | 0.5680 -> 0.5662 | 15.5% -> 6.8% |
| blast_9443_nopfam | 79,545 -> 70,516 | 0.9196 -> 0.9360 | +1.4 -> +4.3 | 0.5677 -> 0.5656 | 12.4% -> 4.2% |
| blast_40674_nopfam | 75,096 -> 64,164 | 0.9138 -> 0.9143 | +0.8 -> +2.1 | 0.5569 -> 0.5495 | 11.9% -> 3.2% |
| blast_7742_nopfam | 73,802 -> 62,087 | 0.9105 -> 0.9037 | +0.5 -> +1.0 | 0.5522 -> 0.5415 | 12.0% -> 3.2% |
| blast_33208_nopfam | 73,306 -> 61,283 | 0.9095 -> 0.9001 | +0.4 -> +0.7 | 0.5506 -> 0.5386 | 12.1% -> 3.2% |
| pfam | 76,151 -> 65,802 | 0.9200 -> 0.9160 | +1.5 -> +2.3 | 0.5626 -> 0.5564 | 12.4% -> 3.8% |
| blast_9606_pfam | 84,745 -> 76,278 | 0.9201 -> 0.9367 | +1.5 -> +4.3 | 0.5704 -> 0.5690 | 15.6% -> 7.0% |
| blast_9443_pfam | 80,127 -> 71,341 | 0.9201 -> 0.9363 | +1.5 -> +4.3 | 0.5700 -> 0.5685 | 12.5% -> 4.3% |
| blast_40674_pfam | 77,241 -> 67,431 | 0.9197 -> 0.9213 | +1.4 -> +2.8 | 0.5654 -> 0.5607 | 12.3% -> 3.9% |
| blast_7742_pfam | 76,525 -> 66,377 | 0.9190 -> 0.9161 | +1.4 -> +2.3 | 0.5635 -> 0.5576 | 12.4% -> 3.9% |
| blast_33208_pfam | 76,373 -> 66,146 | 0.9192 -> 0.9155 | +1.4 -> +2.2 | 0.5631 -> 0.5571 | 12.4% -> 3.8% |

**What changed.**
- Gains over ab initio are about three times larger with the release code: Swiss-Prot minus human,
  MANE exact given assembled +4.4 / +4.5 / +4.3 points (LRAA / StringTie / Trinity; development
  +1.8 / +1.1 / +1.4); MANE compatible +1.5 / +3.5 / +3.3 (development +0.9 / +2.3 / +2.1). Two reasons:
  the release ab initio baseline is lower (the q-value test; E25), and the hits also refine start
  codons. Trinity with Swiss-Prot minus human + Pfam: ab initio 0.8933, the same hits for selection
  only (E33 `final_v7_homology`, 3a793ff) 0.9201 (+2.7), with the start assist 0.9367 (+4.3).
- Minus all animals: +0.7 / +0.7 / +0.7 (development +0.4 / +0.3 / +0.4).
- **Pfam alone no longer recovers nearly all of the gain**: 94.8% / 93.3% / 91.6% exact given
  assembled vs 97.3% / 96.0% / 93.6% with Swiss-Prot minus human (ab initio 92.9 / 91.5 / 89.3), i.e.
  +1.9 / +1.8 / +2.3 of +4.4 / +4.5 / +4.3 points (development: Pfam alone equal to minus human).
  The BLAST trim rule (R2: the subject's N-terminus at the ATG) has no Pfam counterpart.
- With the start assist, BLAST minus vertebrates / minus animals without Pfam end slightly above ab
  initio, while the Pfam configurations stay at about +2 points.
- **Retention artifact** (`gm12878_homology_primate_retention_v7h.txt`): Trinity 5,127 CDSs only with
  minus human (development 4,981); 4,791 unannotated (4,692); 4,506 < 100 aa, 88% of 5,127 (4,752,
  95%); best hit UTY_PANTR 4,555 (4,465), LIN1_NYCCO 201 (185; Swiss-Prot P08548, "LINE-1 reverse
  transcriptase homolog", slow loris; the development text called it an ORF1 protein). MANE exact
  0.3738 vs 0.3737 (0.3671 both). StringTie: 1,532 only with minus human, 1,426 unannotated, 1,353
  UTY_PANTR (1,401 / 1,305 / 1,232). LRAA: 150 / 108 / 95.

### TE screen recipe (PyTD documentation)
- PyTD main **d3892b5** (docs only, on top of a5b24c3): `docs/homology_te_screen.md` and a README
  paragraph. Recipe: ORFs predicted only with homology (Predict again without hits, compare IDs);
  best hit from the first `outfmt6` line and Swiss-Prot names; flag by best-hit name (LINE-1,
  transpos*, reverse transcriptase, endogenous retrovirus, retrovirus-related) and by repeat overlap
  (>= 50% of CDS bases in RepeatMasker SINE/LINE/LTR/DNA/RC/Retroposon intervals, bedtools coverage
  on the genome GFF3's CDS segments); remove with `grep -v -w -F -f` / `seqkit grep -v`.
- Tested on the series with the commands of the document: `pipeline/derive/08_te_screen_check_v7h.sh`
  -> `data/gm12878/gm12878_te_screen_v7h.tsv` (Swiss-Prot minus human, BLAST only, vs ab initio).
  UCSC hg38 `rmsk.txt.gz` downloaded 2026-10-30 (`$P/data/repeats/ucsc_hg38_rmsk/`, 5,683,690 rows,
  gzip integrity checked; UCSC has no md5 file for it); bedtools 2.31.1 in a new env `$P/env/bedtools`.

| assembly | predicted ORFs | ab initio | homology only | flagged by name | by overlap | either | UTY_PANTR best hit | of these genome-mapped | flagged by overlap |
|---|---|---|---|---|---|---|---|---|---|
| LRAA | 35,361 | 33,028 | 2,545 | 60 | 184 | 185 | 116 | 116 | 107 |
| StringTie | 42,315 | 36,479 | 6,229 | 467 | 2,066 | 2,068 | 1,534 | 1,534 | 1,466 |
| Trinity | 101,681 | 80,797 | 21,646 | 2,153 | 7,106 | 7,396 | 6,743 | 5,098 | 4,970 (97%) |

  (Counts are ORF IDs on transcripts, not distinct genomic CDSs. "Homology only" here is ab initio vs
  minus human, a wider set than the minus human vs minus primates comparison above.) The name screen
  misses UTY_PANTR (a host protein), the overlap screen catches 92-97% of its genome-mapped hits (LRAA 92%, StringTie 96%, Trinity 97%);
  1,645 Trinity UTY_PANTR ORFs are on transcripts without a genome mapping and are not screened.

### Manuscript (2026-10-30)
- results_gm12878_parameters.tex: both homology paragraphs rewritten with the release numbers, the
  start-assist share (Trinity, E33), Pfam "about half", the recipe pointer and the 97% screen result;
  the todo removed. Fig. 7 caption (panels c,d: released PyTD a5b24c3; Supplementary Table S7).
  Supplementary Table S7 (`tab:hseries`). Todos in methods_algorithm.tex and main.tex no longer list the
  homology series as development code.

### Caveats (release rerun)
- The Pfam search outputs for LRAA/StringTie come from chunked searches; only E-values (unused) differ
  from a single search.
- Pfam remains not distance-controlled.
