# Benchmark v2: coding-region prediction accuracy on reference transcripts

An update of the 2017–18 benchmark (TransDecoder/TransDecoder_benchmarking_data) with:
- current RefSeq annotation,
- more organisms, including a non-standard genetic code and GC extremes,
- an unstranded setting,
- current comparator tools,
- seeded, scripted generation of every input.

## Truth sets (`scripts/build_truth_from_gbff.py`)
Built from RefSeq `*_rna.gbff`, one transcript per gene (the longest eligible CDS).

Eligible CDS:
- a single, complete, + strand span;
- no exception, ribosomal slippage or transl_except;
- ATG start and a terminal stop of the record's translation table;
- no internal stop;
- translation matches `/translation`;
- sequence is ACGT only.

UTRs are trimmed to ≤250 nt, as in 2018.

| org | assembly | filter | code | genes |
|---|---|---|---|---|
| human | GRCh38.p14 | NM_ + MANE Select | 1 | 19,222 |
| mouse | GRCm39 | NM_ | 1 | 20,990 |
| dmel | Release 6 | NM_ | 1 | 13,870 |
| athal | TAIR10.1 | NM_ | 1 | 27,330 |
| spombe | ASM294v3 | NM_ | 1 | 5,076 |
| tthermophila | JCVI-TTA1-2.2 | XM_ (model) | 6 (ciliate) | 26,190 |
| pfalciparum | 3D7 | XM_ | 1 | 5,281 |
| creinhardtii | v5.5 | XM_ | 1 | 17,621 |

Config: `organisms.tsv`. Raw downloads go in `data/refseq/<org>/` (not in git).

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
| longest_orf | Baseline: the longest LongOrfs ORF per transcript |
| all_orfs | Baseline: every LongOrfs ORF (≥100 aa) |

CodAn has no ciliate or apicomplexan model; the INV model is used for *T. thermophila* and *P. falciparum*. It also cannot use the ciliate genetic code.

GeneMarkS-T needs a licence key and is not yet included.

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
JOBS=16 scripts/run_all.sh data/bench_v2/inputs work/bench_v2/runs  # local, 16 runs at a time
scripts/score_all.py  data/bench_v2/inputs work/bench_v2/runs > results/summary.tsv
```
