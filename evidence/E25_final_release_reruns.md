# E25. Final reruns at the release commits: benchmark and GM12878 (2026-10-24; core runs)

Code: PyTD main 3a793ff (`work/trees/PY_v7`), Perl devel 41fc5ee (`work/trees/TD_v7`). All runs fresh (new run directories; `run_tool.sh` tools `pytd_v7`, `perl_td_v7`, `pytd_v7_precise`, `pytd_v7_sensitive`, `pytd_v7_homology`; GM12878 `genome_annotation/scripts/gm12878_final_v7.sh`). Release defaults: decoy q-value test at q 0.05, GC-table long-ORF accept, 50-codon minimum, start refinement (alpha 0.4, Kozak 17/bit).

Status: phase 1 (pytd_v7, perl_td_v7 on all 84 inputs) and the 7 GM12878 ab initio runs are done and scored here. Phase 2 (`--precise`/`--sensitive`) is done and scored below. Phase 3 (`pytd_v7_homology`, 84 inputs) finished on 2026-10-26 and is scored below (section "Homology"). The GM12878 Trinity homology run (`trinity_final_v7_homology`) has finished but is not yet used in the manuscript.

## Claims supported
- **Perl and Python give identical predictions on all 84 benchmark inputs** at the release commits (`sort preds.tsv | cmp`; 84/84 identical). This extends E24's 7 parity cases to the full benchmark, with the q-value test active.
- **Release medians reproduce the development runs of E20/E22** (`decoy_q0.05_longgc`, earlier RNG) to within 0.001 in sensitivity and 0.07 points in negative-call rate:

  | input (median over organisms) | release (pytd_v7) stop / exact / % negatives called | development q 0.05 | development fixed short-ORF test (cliff_c200) |
  |---|---|---|---|
  | complete + random (8) | 0.964 / 0.899 / 0.29 | 0.965 / 0.899 / 0.36 | 0.957 / 0.897 / 0.48 |
  | fragments + random (8) | 0.872 / 0.772 / 0.23 | 0.872 / 0.773 / 0.21 | 0.815 / 0.708 / 0.11 |
  | complete + lncRNA (human, mouse) | 0.963 / 0.893 / 5.04 | 0.964 / 0.894 / 4.95 | 0.976 / 0.905 / 18.93 |
  | fragments + lncRNA (human, mouse) | 0.858 / 0.756 / 5.89 | 0.857 / 0.755 / 5.83 | 0.858 / 0.753 / 18.65 |

  (`negatives_scores.tsv`, `score_negatives.py`.)
- **`--precise` and `--sensitive` (phase 2, 40 runs, none failed) also reproduce the development runs** (stop / exact / % negatives called; release vs development):

  | setting | complete + random | fragments + random | complete + lncRNA | fragments + lncRNA |
  |---|---|---|---|---|
  | `--precise` (q 0.01) | 0.955 / 0.888 / 0.29 (dev 0.956 / 0.888 / 0.33) | 0.837 / 0.729 / 0.17 (0.836 / 0.727 / 0.19) | 0.950 / 0.881 / 4.25 (0.950 / 0.881 / 4.09) | 0.794 / 0.703 / 4.56 (0.790 / 0.698 / 4.52) |
  | `--sensitive` (q 0.2) | 0.971 / 0.905 / 0.40 (0.970 / 0.905 / 0.42) | 0.895 / 0.798 / 0.45 (0.896 / 0.799 / 0.43) | 0.974 / 0.904 / 7.80 (0.974 / 0.904 / 7.67) | 0.906 / 0.795 / 10.23 (0.909 / 0.797 / 10.51) |
- **All 84 inputs, F1 against the previous main-benchmark default (`pytd_dev_v2`, fixed short-ORF test)**, medians (`all_scores.tsv`, `score_predictions.py`):

  | inputs | stop F1 | exact F1 | stop sensitivity | stop precision |
  |---|---|---|---|---|
  | complete (34) | 0.977 vs 0.973 | 0.907 vs 0.906 | 0.962 vs 0.963 | 0.994 vs 0.992 |
  | fragments (34) | 0.921 vs 0.890 | 0.812 vs 0.784 | 0.865 vs 0.819 | 0.989 vs 0.993 |
  | complete, full UTR (16) | 0.968 vs 0.968 | 0.897 vs 0.900 | 0.946 vs 0.945 | 0.990 vs 0.986 |

  Exact F1 improved on 61 inputs and fell on 17 (by at most 0.64 points: T. thermophila complete −0.64, mouse complete_fullutr −0.61, human complete_fullutr.unstranded −0.48). Largest gains: human complete.wLnc +9.8, human fragments.wLnc +9.2, mouse complete.wLnc +6.4 points.
- **Run time** (shared host, single thread; indicative): Perl/PyTD ratio median 2.71 (2.07–3.28) over 84 inputs. PyTD v7 took a median 1.75× the time of `pytd_dev_v2`, mainly the decoy step (LongOrfs + scoring of one shuffled copy of every transcript).

## Homology (phase 3; `pytd_v7_homology` vs `pytd_v7`)
Swiss-Prot 2026_03 without the organism's own entries (DIAMOND, one target, E ≤ 1e-5) and Pfam-A 38.2 (hmmsearch, gathering thresholds), on all LongOrfs ORFs; otherwise the release defaults. 84/84 runs finished, 0 failed. Scored with `score_v7.sh` (rows appended to `data/benchmark/final_v7/all_scores.tsv`); per-input differences from `compare_homology.py data/benchmark/final_v7/all_scores.tsv $P/data/bench_v3/inputs $P/work/bench_v3/runs pytd_v7 pytd_v7_homology > data/benchmark/homology_vs_abinitio.tsv`. The development comparison it replaces is kept as `homology_vs_abinitio.dev_v2.tsv` (E07).

Claims (Results, "Homology evidence recovers true coding regions, mostly in fragments"; Fig. 5d,e; checked in notebook 05):
- Complete-type inputs (complete, complete_fullutr, with and without random negatives; 48): stop sensitivity −0.29 to +3.72 points, exact F1 −0.35 to +2.09 points. Largest: human/mouse complete_fullutr(.unstranded); smallest: *P. falciparum* and *C. reinhardtii* complete_fullutr.unstranded.
- Fragments (32): stop sensitivity +0.42 to +9.20 points, exact F1 +0.11 to +4.65 points (largest: mouse fragments.unstranded(.wRand)).
- Stranded fragments, exact F1: mouse +3.08, human +2.58, *A. thaliana* +2.34; *T. thermophila* +0.65, *S. pombe* +0.56, *D. melanogaster* +0.39, *C. reinhardtii* +0.22, *P. falciparum* +0.13.
- Stop precision on inputs without lncRNAs: −0.51 to +1.32 points. Calls on random negatives: +0 to +2 per input in seven organisms, +12 to +39 in *T. thermophila*.
- lncRNA negatives (complete.wLnc): lncRNAs with a prediction, human 1,639 → 3,268, mouse 1,640 → 2,452 (fragments.wLnc: 1,928 → 3,522 and 1,904 → 2,698). Exact F1: complete.wLnc human −2.11, mouse −0.10; fragments.wLnc human +0.21, mouse +2.30 points.

Comparison with the development runs (E07, `pytd_dev_v2` with the fixed short-ORF test): the fragment gains are smaller at release (stranded fragments, exact F1 +3.8 to +4.7 in human/mouse/*A. thaliana* before, +2.3 to +3.1 now; +0.6 to +5.3 over all fragments before), and the lncRNA cost is larger in relative terms. The ab initio release run predicts on far fewer lncRNAs (1,639 vs 6,943 in human), so the ~1,600 lncRNAs added by homology double the count, whereas at development they added 884 to an already larger number. Homology-supported ORFs bypass the q-value test, so this added share is expected once that test removes most chance ORFs. Which lncRNA ORFs gain a homology hit (e.g. pseudogene-derived or paralogous sequence) was not examined.

## GM12878 (ab initio, release default; vs the earlier `final_abinitio` runs with the fixed short-ORF test)
Same transcripts as `final_abinitio`. Eval: `eval_genome_cds.py` (`eval.*.summary.tsv`).

| assembly | distinct CDSs | MANE exact (of assembled) | no GENCODE CDS overlap | unannotated ("novel") CDSs |
|---|---|---|---|---|
| LRAA | 22,769 vs 23,922 | 5,848 vs 5,882 (0.929 vs 0.935) | 1.5% vs 3.7% | 341 vs 885 |
| StringTie | 30,846 vs 34,404 | 8,160 vs 8,268 (0.915 vs 0.927) | 5.5% vs 12.0% | 1,687 vs 4,127 |
| Trinity | 59,514 vs 72,323 | 6,857 vs 6,951 (0.893 vs 0.906) | 4.2% vs 15.5% | 2,528 vs 11,235 |

Settings (LRAA / StringTie), distinct CDSs and no-overlap share: precise 21,984 / 28,947, 1.2% / 4.1%; default as above; sensitive 23,468 / 33,838, 2.3% / 9.6%. MANE exact: precise 5,737 / 8,043; sensitive 5,907 / 8,264.

**Protein matches** (`score_protein_matches.py`, GENCODE v50 translations; matching = identical or shared exact 30-mer; `gm12878_v7_protein_matches.tsv`):

| run | matching | non-matching | E22 development run at same q |
|---|---|---|---|
| LRAA default | 32,751 | 277 | 32,763 / 278 |
| LRAA precise | 31,672 | 189 | 32,011 / 193 |
| LRAA sensitive | 33,443 | 560 | 33,508 / 579 |
| StringTie default | 35,294 | 1,185 | 35,391 / 1,195 |
| StringTie precise | 33,647 | 698 | 33,869 / 749 |
| StringTie sensitive | 36,903 | 2,879 | 36,903 / 2,832 |
| Trinity default | 78,181 | 2,616 | — |

Earlier fixed short-ORF test: LRAA 33,587 / 1,005; StringTie 36,547 / 4,119.

**Ribo-seq frame support** (E01 method, same BAM and calibration procedure; `gm12878_v7_riboseq_support.tsv`): supported fraction of covered CDSs, all predictions: LRAA 0.657 (earlier 0.652), StringTie 0.615 (0.610), Trinity 0.612 (0.605); exact class 0.663 / 0.647 / 0.642. Unannotated CDSs with ≥ 10 P-sites remain rare (21 / 29 / 38).

**MS peptides** (E02 method; `gm12878_v7_ms_peptide_support.tsv`): unannotated CDSs with ≥ 1 peptide: LRAA 167 of 341 (49%; earlier 208 of 885, 23.5%), StringTie 748 of 1,687 (44%; earlier 849 of 4,127, 21%), Trinity 765 of 2,528 (30%; earlier 1,018 of 11,235, 9%). Exact class 93.6–94.5% (earlier 93.1–94.2%). Per E02, peptides on unannotated CDSs mostly come from paralogous or pseudogene copies of known genes, so these rates measure paralog/pseudogene content as much as translation; novel-specific peptides: 2 / 8 / 17.

**Cost: catalogued Ribo-seq ORFs** (GENCODE phase 1+2 catalogue matches by stop codon, `match_riboseq.py`; `gm12878_v7_riboseq_orf_catalog_matches.tsv`): release default 101 (LRAA), 95 (StringTie), 186 (Trinity), against 178 / 217 / 455 with the fixed short-ORF test and 184 / 221 with a 100-codon minimum (E04). Precise 79 / 81; sensitive 144 / 165. The q-value test removes more catalogued translated small ORFs than the length rules did; this must be stated in Results/Discussion.

## Commands
- Benchmark launcher: `work/bench_v3/launch_v7.sh` (log `runs/logs/launch_v7.out`: 170 done, 0 failed at scoring time).
- Scores: `score_negatives.py data/bench_v3/inputs complete.wRand,fragments.wRand,complete.wLnc,fragments.wLnc pytd_v7=… precise=… sensitive=… dev_q0.05=rules/…/decoy_q0.05_longgc … cliff_c200=…`; per-input F1 via `score_v7.sh` (wraps `score_predictions.py`).
- GM12878: `gm12878_final_v7.sh` (predict, genome map, `eval_genome_cds.py`, `match_riboseq.py`); then `riboseq_orf_support.py` and `ms_peptide_support.py` with the E01/E02 arguments on the v7 GFF3/pep files (`work/gm12878_td/final_v7_scoring/`).

## Caveats
- GM12878 MANE exact recovery falls 0.5–1.2 points with the q-value test (exact given assembled: LRAA −0.5, StringTie −1.2, Trinity −1.2 points), the price of 2.4- to 4.4-fold fewer unannotated CDSs.
