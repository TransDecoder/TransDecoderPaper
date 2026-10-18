# Benchmark v2: coding-region prediction accuracy on reference transcripts

An update of the 2017–18 benchmark (TransDecoder/TransDecoder_benchmarking_data) with:
- the most current reference annotation for each organism,
- more organisms, including a non-standard genetic code and GC extremes,
- an unstranded setting,
- current comparator tools,
- seeded, scripted generation of every input.

## Truth sets (`scripts/build_truth_from_annotation.py`, `scripts/build_all_truth.sh`)
Built from each organism's most current reference genome annotation (genome FASTA + GFF3/GTF; `truth_sources.tsv`, download URLs in `data/annot/sources.tsv`). One transcript per gene: the tagged canonical transcript where one exists (MANE Select, Ensembl canonical), otherwise the longest eligible CDS.

Eligible CDS:
- nuclear, protein-coding; CDS contained in the exons and contiguous in the transcript;
- ATG start and a terminal stop of the organism's translation table;
- length divisible by 3, no internal stop (excludes selenoproteins and read-through);
- transcript sequence is ACGT only.

**Missing UTRs.** Transcripts annotated without a 5′ (or 3′) UTR get one. Its length is drawn from that organism's empirical 5′ (3′) UTR lengths among kept transcripts that have one, and it is filled with the genomic sequence flanking the transcript (unspliced). `utr5_source` / `utr3_source` in the truth table record `annotated`, `added` or `none`, so results can be stratified.

UTRs are trimmed to ≤250 nt, as in 2018 (`truth/`); `truth_fullutr/` keeps them whole.

| org | annotation | selection | code | genes | 5′ / 3′ UTR added |
|---|---|---|---|---|---|
| human | MANE v1.5 (RefSeq RS_2025_08 / Ensembl 116), GRCh38 | MANE Select (64 on alt/fix patches excluded) | 1 | 19,167 | 351 / 335 |
| mouse | GENCODE M39 (Ensembl 116), GRCm39 | Ensembl canonical, no `*_NF` tags | 1 | 21,409 | 834 / 814 |
| dmel | FlyBase FB2026_03 (r6.69) | longest CDS | 1 | 13,884 | 517 / 414 |
| athal | TAIR12 (TAIR GFF3 via AraBase mirror) | longest CDS | 1 | 26,867 | 1,142 / 912 |
| spombe | PomBase (2026-10-17) | longest CDS | 1 | 5,100 | 599 / 383 |
| tthermophila | TGD v6 (Ye et al. 2025), main transcripts | one per gene | 6 (ciliate) | 26,456 | 279 / 316 |
| pfalciparum | RefSeq RS_2026_09 (PlasmoDB-derived), 3D7 | longest CDS | 1 | 5,281 | 875 / 947 |
| creinhardtii | RefSeq v5.5 (latest NCBI; v6.1 only via JGI login) | longest CDS | 1 | 17,621 | 0 / 0 |

Notes:
- TGD's "v6 assembly" download does not match its v6 GFF3; the genome-browser FASTA (`1-upd-Genome-assembly.fasta`) does and is used.
- NCBI's copy of TAIR12 (GCA_978657495.1) lacks the CDS for 6,562 isoforms, so the TAIR-produced GFF3 is used.

The superseded RefSeq-gbff truth sets (`scripts/build_truth_from_gbff.py`, `refseq_sources.tsv`) are kept for reference; 85% of the *T. thermophila* RefSeq models there had no UTRs.

## Input variants (`scripts/make_benchmark_inputs.py`, seed 1)
Every organism gets eight variants: {`complete`, `fragments`} × {stranded, `.unstranded`} × {without, with `.wRand`}.
- **complete**: full transcripts, CDS on the + strand.
- **fragments**: each transcript is cut once inside the CDS at a fraction ~ N(0.5, 0.1).
  - If the cut is 3′ of the CDS midpoint, the 3′ piece is kept (5′-partial truth).
  - Otherwise the 5′ piece is kept (3′-partial truth).
  - Partial truth coordinates are trimmed to whole codons.
- **.unstranded**: a random half of the transcripts are reverse-complemented. Tools predict on both strands.
- **.wRand**: adds one random negative per transcript. It has the same length and the same G+C fraction as its transcript, with i.i.d. bases. Any prediction on a negative is a false positive.

Stranded inputs are run strand-specifically:
- TransDecoder: `-S`
- Prodigal: + strand calls only
- CodAn: `--strand plus`

## Tools (`scripts/run_tool.sh`)
| label | what |
|---|---|
| pytd | PyTransDecoder, current main (default settings) |
| pytd_norefine | PyTransDecoder `--no-refine-starts` |
| perl_td | Perl TransDecoder, current master (with the 2026 fixes) |
| perl_td_v600 | Perl TransDecoder v6.0.0 (before the fixes) |
| perl_td_v571 | Perl TransDecoder v5.7.1 |
| prodigal | Prodigal 2.6.3, single mode, `-n`, `-g <table>`, as in 2018 |
| codan | CodAn 1.1 with its pre-trained clade model (`*_full` for complete, `*_partial` for fragments) |
| gmst | GeneMarkS-T 5.1 (self-training; licensed, installed outside the repo with the key in `~/.gm_key`). It supports genetic codes 1/4/11 only, so Tetrahymena runs with the standard code |
| cpat | CPAT 3.0.5 (logistic regression on hexamer/Fickett/ORF features; the ORF step of SQANTI3). Prebuilt human/mouse/fly model and its published cutoff (human for organisms without one); `--antisense` when unstranded; ciliate stops via `--stop TGA` |
| cpc2 | CPC2 1.0.1 (`--ORF`; longest ATG ORF of transcripts labelled coding; `-r` when unstranded; used with CPAT by IsoformSwitchAnalyzeR) |
| borf | Borf 1.2, defaults (`-s` when unstranded) |
| rnasamba | RNAsamba 0.2.5 with the authors' human weights (full-length / partial). Forward strand only, so stranded inputs only; proteins mapped back to transcript coordinates |
| angel | ANGEL 3.0 (PacBio; deprecated upstream), self-training workflow from its README (`dumb_predict` → `angel_make_training_set --random` → `angel_train` → `angel_predict --output_mode best`, ≥100 aa) |
| longest_orf | Baseline: the longest LongOrfs ORF per transcript |
| all_orfs | Baseline: every LongOrfs ORF (≥100 aa) |

Only TransDecoder, Prodigal and CPAT take the ciliate genetic code; GeneMarkS-T, CodAn, CPC2, Borf, RNAsamba and ANGEL run *T. thermophila* with the standard code.
ORFanage (reference-guided) is not included: it needs a reference annotation of the same genome. NCBI ORFfinder and orfipy enumerate ORFs like the `all_orfs` baseline.

CodAn has no ciliate or apicomplexan model; the INV model is used for *T. thermophila* and *P. falciparum*. It also cannot use the ciliate genetic code.

## Scoring (`scripts/score_predictions.py`, `scripts/score_all.py`)
Matching follows the 2018 benchmark:
- **Stop level**: same strand, and the 3′ end is within 2 nt.
- **Exact level**: additionally, the 5′ end is within 2 nt.

Each truth CDS is credited at most once.

Metrics:
- Sensitivity = TP / truth.
- Precision = TP / all predictions, including those on random negatives. (2018 called this "specificity".)
- F1.

## Running
```
scripts/build_all_truth.sh data/annot data/bench_v3
scripts/make_all_inputs.sh data/bench_v3
JOBS=16 scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs  # local, 16 runs at a time
scripts/score_all.py  data/bench_v3/inputs work/bench_v3/runs > results/summary.tsv
```
