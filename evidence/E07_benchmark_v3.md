# E07. Cross-species benchmark (bench_v3)

**Status: complete except CodAn scores.** CodAn has 71/80 runs finished as of 2026-10-22 11:40. Nine are still running: 6 *T. thermophila* variants (complete.unstranded, complete_fullutr.unstranded, complete.wRand, complete.unstranded.wRand, fragments.wRand, fragments.unstranded.wRand), athal and mouse `complete.unstranded.wRand`, and mouse `fragments.unstranded.wRand`. The committed `summary.tsv` holds 62 CodAn rows. The homology-assisted final-default runs (`pytd_dev_v2_homology`, 84/84) were added on 2026-10-22; see Results.

## Claims supported
Manuscript: Methods, "Benchmark data", "Compared methods", "Accuracy measures" (`manuscript/sections/methods_benchmark.tex`); Results, benchmark (not yet drafted).

1. A benchmark was built on 8 eukaryotes from each organism's current reference annotation.
   - The organisms include the ciliate genetic code (*T. thermophila*) and GC extremes: *P. falciparum* 22% and *C. reinhardtii* 68% transcript GC.
   - The inputs are complete transcripts, full-UTR transcripts and fragments; stranded and unstranded; with and without GC-matched random negatives.
2. TransDecoder (Perl and PyTD, final defaults) was compared with:
   - its own earlier releases;
   - Prodigal, CodAn, GeneMarkS-T, CPAT, CPC2, Borf and RNAsamba;
   - the longest-ORF baseline. An all-ORFs (>= 100 codons) baseline was also scored (in summary.tsv) but is not reported (Brian, 2026-10-23: leave it out).
   - The Perl and Python final defaults give identical predictions on all 84 inputs (see E08).
3. Matching follows the 2018 evaluation: the stop within 2 nt for a stop-level match, and the start within 2 nt as well for an exact match. Each truth CDS is credited once. Precision counts predictions on random negatives.
4. Headline accuracy is given under Results below. The minimum-length sweep and lncRNA negatives are in E04, start refinement in E09, and the long-ORF rule in E10.

## Data
- **Annotations:** URLs and releases are in `$T/benchmark_v2/annotation_downloads.tsv` (identical to `$P/data/annot/sources.tsv`). They were downloaded 2026-10-18 to `$P/data/annot/<org>/` (`download.done`: Sun Oct 18 09:13 EDT 2026).

  | org | annotation | selection | code | genes | 5′ / 3′ UTR added | transcript GC% (complete.fa) |
  |---|---|---|---|---|---|---|
  | human | MANE v1.5 on GENCODE v50 GRCh38 primary | MANE Select | 1 | 19,167 | 351 / 335 | 52.3 |
  | mouse | GENCODE M39 (Ensembl 116), GRCm39 | Ensembl canonical, no `*_NF` tags | 1 | 21,409 | 834 / 814 | 51.7 |
  | dmel | FlyBase FB2026_03 (r6.69) | longest CDS | 1 | 13,884 | 517 / 414 | 50.9 |
  | athal | TAIR12 (AraBase mirror of TAIR GFF3) | longest CDS | 1 | 26,867 | 1,142 / 912 | 42.2 |
  | spombe | PomBase 2026-10-17 | longest CDS | 1 | 5,100 | 599 / 383 | 37.9 |
  | tthermophila | TGD v6 main transcripts | one per gene | 6 (ciliate) | 26,456 | 279 / 316 | 25.9 |
  | pfalciparum | RefSeq GCF_000002765.6 (RS_2026_09) | longest CDS | 1 | 5,281 | 875 / 947 | 22.1 |
  | creinhardtii | RefSeq GCF_000002595.2 v5.5 | longest CDS | 1 | 17,621 | 0 / 0 | 68.3 |

  Gene counts and UTR additions are from `$T/benchmark_v2/truth_logs_v3/<org>.log` and match `benchmark_v2/README.md`. Builder arguments per organism are in `truth_sources.tsv`.

- **Truth sets:**
  - `$P/data/bench_v3/truth/` has UTRs trimmed to ≤250 nt. `truth_fullutr/` keeps full UTRs.
  - Files are `<org>.fasta`, `<org>.truth.tsv` and `<org>.log`.
  - md5 of `truth/<org>.truth.tsv`: athal 3f8bcf3e…, creinhardtii 626ca1eb…, dmel 556b8b2e…, human e06135fa…, mouse c94939be…, pfalciparum 591d56bb…, spombe bf47c928…, tthermophila c7d60b9c…
- **Inputs:** `$P/data/bench_v3/inputs/<org>/<variant>.{fa,truth.tsv}`, made with seed 1. There are 10 variants per organism:
  - complete, complete.unstranded, complete.wRand, complete.unstranded.wRand;
  - fragments, fragments.unstranded, fragments.wRand, fragments.unstranded.wRand;
  - complete_fullutr, complete_fullutr.unstranded.

  `.wRand` inputs have exactly 2× the transcripts: one i.i.d. GC-matched negative per transcript, named `rand_*`.

  Human and mouse also have `complete.wLnc` and `fragments.wLnc` (GENCODE lncRNA negatives; see E04).
- **Homology databases** (`pytd_homology` only), in `$P/data/homology/`:
  - Swiss-Prot 2026_03 (02-Sep-2026) without the organism's own taxid: `sprot_minus_<taxid>.dmnd`;
  - Pfam-A 38.2 (2026-01).

## Software
Benchmark runs are in `$P/work/bench_v3/runs/<org>/<variant>/<tool>/`. Labels, trees and commits are below; run times are the `preds.tsv` mtimes.

| label | code | commit at run time | runs |
|---|---|---|---|
| `pytd_dev_v2` (**final defaults**) | PyTD `work/trees/PY_dev` (branch dev/frame-contrast-starts) | 48c870c (see Caveats) | 84, 2026-10-20 12:41–15:13 |
| `perl_td_dev_v2` (**final defaults**) | TransDecoder `work/trees/TD_fix` (branch fix/start-refinement) | de3f57e (see Caveats) | 84, 2026-10-20 12:45–15:32 |
| `pytd_dev_v2_homology` (final defaults + DIAMOND/Pfam) | PyTD `work/trees/PY_dev`; 5 *T. thermophila* runs re-predicted from `work/trees/PY_48c870c` (see Caveats) | 48c870c | 84, launch17, finished 2026-10-22 |
| `pytd`, `pytd_norefine`, `pytd_homology` | PyTD `work/trees/PY_fix` (dev/paper; = origin/main) | 27bfb5e | 80/80/56, 10-18 09:37 – 10-19 11:16 |
| `perl_td` | `TD_fix` | 351634a (tree HEAD on 10-18; the tree has since moved to de3f57e) | 80, 10-18 |
| `perl_td_v600` | `work/trees/TD_base`, v6.0.0 | 66d4712 | 80 |
| `perl_td_v571` | `work/trees/TD_v571`, v5.7.1 | 5791fde | 80 |
| prodigal | Prodigal V2.6.3 (Feb 2016), `$P/env/td` | – | 84 |
| codan | CodAn 1.1 (bioconda codan-1.1; Python 3.6.15, Biopython 1.69), `$P/env/codan`; models `$P/tools/CodAn/models` (git 43da738) | – | 68/80 |
| gmst | GeneMarkS-T (`gmst.pl --version`: "GeneMarkS version 5.1 March 2014"), in `~/software/genemarks_t/dist`. Licensed; it needs a GeneMark licence key in the user's home directory. Neither the key nor the software is distributed here | – | 84 |
| cpat | CPAT 3.0.5 (`$P/env/cpat`); prebuilt Human/Mouse/Fly models in `$P/tools/models/cpat` | – | 84 |
| cpc2 | CPC2 1.0.1 (`$P/env/cpc2`) | – | 84 |
| borf | Borf 1.2 (`$P/env/borf`) | – | 84 |
| rnasamba | RNAsamba 0.2.5 (TensorFlow 1.14), human weights in `$P/tools/models/rnasamba` | – | 44 (stranded inputs only) |
| angel | ANGEL 3.0 (`$P/tools/ANGEL`, git b805554; `$P/env/angel`) | – | 3 (dropped) |
| longest_orf, all_orfs | baselines from the `pytd` LongOrfs output | 27bfb5e | 80 each |

- **Common environment:**
  - Python 3.11.17 with Biopython 1.88 and Perl 5.32.1 (`$P/env/td`); DIAMOND 2.2.8; HMMER 3.4.
  - `PERL_HASH_SEED=0`.
- **Paper repo:** `$T` at 774429e when this record was written.
- **Development runs:** the intermediate labels (`pytd_dev`, `pytd_dev_pwm`, `pytd_dev_single`, `pytd_dev_kozak`, `pytd_dev_kozak_add`, `perl_td_dev`, and `pytd_dev_m*`) are scored in the same table. They are covered in E09 and E04.

## Commands
All scripts are in `$T/benchmark_v2/scripts/`. Each tool's exact command line is in `run_tool.sh`.
```bash
cd $P
bash TransDecoderPaper/benchmark_v2/scripts/download_annotations.sh data/annot     # sources.tsv = annotation_downloads.tsv
PYTHON=env/td/bin/python TransDecoderPaper/benchmark_v2/scripts/build_all_truth.sh data/annot data/bench_v3
#   -> build_truth_from_annotation.py --add-utr --max-utr {250|-1} [builder_args from truth_sources.tsv]
PYTHON=env/td/bin/python TransDecoderPaper/benchmark_v2/scripts/make_all_inputs.sh data/bench_v3
#   -> make_benchmark_inputs.py --seed 1 (complete/fragments x stranded/unstranded x +/-wRand; complete_fullutr from truth_fullutr)
bash TransDecoderPaper/benchmark_v2/launch_v3/launch.sh     # pytd pytd_norefine perl_td perl_td_v600 perl_td_v571 prodigal codan gmst; pytd_homology
bash TransDecoderPaper/benchmark_v2/launch_v3/launch2.sh    # second runner, same tools minus CodAn
bash TransDecoderPaper/benchmark_v2/launch_v3/launch3.sh    # cpat cpc2 borf; rnasamba on stranded variants
bash TransDecoderPaper/benchmark_v2/launch_v3/launch4.sh    # angel (stopped, see Caveats)
bash TransDecoderPaper/benchmark_v2/launch_v3/launch14.sh   # comparators on the wLnc inputs
bash TransDecoderPaper/benchmark_v2/launch_v3/launch16.sh   # final defaults pytd_dev_v2 / perl_td_dev_v2 on all inputs
bash TransDecoderPaper/benchmark_v2/launch_v3/launch17.sh   # pytd_dev_v2_homology (JOBS=3 THREADS=4)
R=work/bench_v3/runs/tthermophila                           # Predict re-run at 48c870c, reusing the searches
TransDecoderPaper/benchmark_v2/scripts/repredict_homology.sh work/trees/PY_48c870c \
  $R/{complete.wRand,complete.unstranded.wRand,fragments.unstranded,fragments.unstranded.wRand,fragments.wRand}/pytd_dev_v2_homology
env/td/bin/python TransDecoderPaper/benchmark_v2/scripts/score_all.py data/bench_v3/inputs work/bench_v3/runs \
  > work/bench_v3/scores/summary.tsv          # copied to $T/benchmark_v2/results_v3/summary.tsv
```
- **Launch scripts:** each `launch*.sh` in `launch_v3/` is a copy of `$P/work/bench_v3/launch*.sh`. Each runs `run_all.sh` locally with `xargs -P JOBS`. `launch5/9/10–13/15` are development runs (E04, E09, E10).
- **Before launch16,** the earlier final-default runs were moved to `$P/work/bench_v3/rerun_old16/` (168 runs; their summary is `rerun_old16/summary.pre16.tsv`). This step is not in the launch script (see Gaps).
- **Settings (from `run_tool.sh`):**
  - **Strandedness:**
    - stranded inputs: TD `-S`, Prodigal `--plus-only` filter, CodAn `-s plus`, GMST `--strand direct`;
    - unstranded: CPAT `--antisense`, CPC2 `-r`, Borf `-s`, ANGEL `--use_rev_strand`.
  - **Genetic code:** TD uses `-G` from `organisms.tsv` (Tetrahymena for tthermophila) and Prodigal `-g` (6). GMST is given `--gcode 1`; CPAT gets `--stop TGA` for table 6. CodAn, CPC2, Borf, RNAsamba and ANGEL use the standard code.
  - **CodAn models:** clade `*_full` / `*_partial` (VERT, INV, PLANTS, FUNGI). INV is used for tthermophila and pfalciparum. Input is line-wrapped at 80 nt.
  - **CPAT cutoffs:** Human 0.364, Mouse 0.44, Fly 0.39. The Human model is used for organisms without their own model.

## Outputs
- `$T/benchmark_v2/results_v3/summary.tsv`: 1,957 scored runs, long format. Columns: n_truth, n_pred, tp_stop, tp_exact, fp, fp_random, fp_wrong_strand, fp_dup, fn, sensitivity/precision/F1 at stop and exact level, and runtime_sec.
  - It is byte-identical to `$P/work/bench_v3/scores/summary.tsv` (written 2026-10-20 15:37).
- `$T/benchmark_v2/results_v3/homology_vs_abinitio.tsv`: `pytd_dev_v2_homology` minus `pytd_dev_v2` per input (84 rows), from `compare_homology.py` on the interim rescoring `$P/work/bench_v3/scores/summary.interim_20261022.tsv` (fixed scorer, all runs as of 2026-10-22 11:30; not committed, as CodAn is incomplete).
- `results_v3/by_length.tsv`, `short_orf_summary.tsv` and `lnc_by_era.human.tsv` are covered in E04.
- Per-run outputs (not in git): `$P/work/bench_v3/runs/<org>/<variant>/<tool>/{preds.tsv,runtime_sec.txt,*.log}`. Logs are in `runs/logs/`.

## Results
**Ranges across the 8 organisms** (`summary.tsv`). F1 values are min / median / max. Tools that cannot use the ciliate code score near zero on *T. thermophila*: GMST, CPC2, Borf and RNAsamba, F1 < 0.02. That is the source of the low minima marked *.

| variant | tool | stop F1 | exact F1 |
|---|---|---|---|
| complete | PyTD/Perl final defaults | 0.940 / 0.978 / 0.988 | 0.742 / 0.910 / 0.928 |
| complete | pytd (27bfb5e, PWM starts) | 0.946 / 0.972 / 0.985 | 0.484 / 0.728 / 0.906 |
| complete | perl_td_v600 (= v571) | 0.948 / 0.972 / 0.985 | 0.567 / 0.779 / 0.912 |
| complete | prodigal | 0.856 / 0.937 / 0.978 | 0.395 / 0.643 / 0.884 |
| complete | gmst | 0.008* / 0.973 / 0.987 | 0.001* / 0.770 / 0.925 |
| complete | codan (7 org.) | 0.710 / 0.974 / 0.983 | 0.378 / 0.830 / 0.887 |
| complete | cpat | 0.884 / 0.965 / 0.989 | 0.786 / 0.940 / 0.969 |
| complete | rnasamba | 0.014* / 0.990 / 0.993 | 0.002* / 0.953 / 0.982 |
| complete | longest_orf | 0.954 / 0.976 / 0.984 | 0.469 / 0.689 / 0.897 |
| complete.wRand | PyTD/Perl final defaults | 0.930 / 0.975 / 0.985 | 0.740 / 0.908 / 0.928 |
| complete.wRand | prodigal | 0.502 / 0.791 / 0.927 | 0.232 / 0.552 / 0.838 |
| complete.wRand | cpat | 0.762 / 0.946 / 0.971 | 0.612 / 0.903 / 0.945 |
| fragments | PyTD/Perl final defaults | 0.853 / 0.895 / 0.924 | 0.723 / 0.789 / 0.812 |
| fragments | pytd (27bfb5e) | 0.821 / 0.869 / 0.886 | 0.622 / 0.715 / 0.804 |
| fragments | prodigal | 0.869 / 0.936 / 0.978 | 0.612 / 0.782 / 0.925 |
| fragments | codan (7 org.) | 0.588 / 0.912 / 0.951 | 0.465 / 0.761 / 0.827 |
| fragments | cpat | 0.573 / 0.783 / 0.915 | 0.322 / 0.440 / 0.506 |
| fragments.wRand | PyTD/Perl final defaults | 0.852 / 0.891 / 0.923 | 0.723 / 0.789 / 0.811 |
| fragments.wRand | prodigal | 0.619 / 0.803 / 0.934 | 0.435 / 0.676 / 0.886 |

Per organism, stop F1 / exact F1 for the final defaults (`pytd_dev_v2`):

| org | complete | complete.wRand | fragments | fragments.wRand |
|---|---|---|---|---|
| athal | 0.974 / 0.897 | 0.971 / 0.894 | 0.853 / 0.723 | 0.852 / 0.723 |
| creinhardtii | 0.968 / 0.742 | 0.967 / 0.740 | 0.908 / 0.760 | 0.907 / 0.760 |
| dmel | 0.972 / 0.922 | 0.971 / 0.921 | 0.878 / 0.794 | 0.877 / 0.794 |
| human | 0.988 / 0.922 | 0.985 / 0.919 | 0.924 / 0.812 | 0.923 / 0.811 |
| mouse | 0.984 / 0.908 | 0.982 / 0.906 | 0.911 / 0.802 | 0.911 / 0.802 |
| pfalciparum | 0.985 / 0.928 | 0.984 / 0.928 | 0.895 / 0.806 | 0.893 / 0.806 |
| spombe | 0.982 / 0.912 | 0.979 / 0.909 | 0.890 / 0.784 | 0.889 / 0.784 |
| tthermophila | 0.940 / 0.840 | 0.930 / 0.830 | 0.896 / 0.762 | 0.887 / 0.755 |

What can be said without over-reading:
- At the stop level on complete transcripts, TransDecoder and most comparators fall in a narrow band (median 0.94–0.99). The differences are in:
  - start accuracy (exact F1);
  - fragments;
  - precision when random negatives are added. Prodigal's stop F1 falls from 0.937 to 0.791 (median) and the final defaults' from 0.978 to 0.975.
- **Fragments, stop level:** Prodigal is higher than TransDecoder in 7 of 8 organisms; *T. thermophila* is the exception (0.869 vs 0.896). This is consistent with its lower minimum gene length (90 nt) against the -m 50/100-codon floor; see E04 for length-stratified numbers. No significance testing was done.
- **Classifier-type tools** (CPAT, RNAsamba, CPC2) have high exact F1 on complete transcripts but low exact F1 on fragments (CPAT median 0.44).
- **CodAn** stop F1 is 0.710 on *P. falciparum*, which has no clade model (INV used).
- **Homology** (`pytd_dev_v2_homology` vs `pytd_dev_v2`, `homology_vs_abinitio.tsv`). Swiss-Prot without the organism's own entries plus Pfam, at the final defaults:
  - Complete-type inputs (complete and complete_fullutr, 48 inputs): stop sensitivity +0.0008 to +0.0246, exact F1 −0.0013 to +0.0180. Gains are smallest in *T. thermophila* and *C. reinhardtii* (exact F1 ≤ +0.0021).
  - Fragments (32 inputs): stop sensitivity +0.013 to +0.099, exact F1 +0.006 to +0.053. Largest in human, mouse and athal (exact F1 +0.040 to +0.053); smallest in *T. thermophila* and *C. reinhardtii* (+0.006 to +0.008).
  - Precision changes by −0.006 to +0.018 (stop level), so the gains come from recovered true CDSs. Predictions on GC-matched random negatives barely move: +0 to +3 per input in seven organisms, +12 to +36 in *T. thermophila* (e.g. complete.unstranded.wRand 1,059 → 1,095).
  - lncRNA negatives: lncRNAs with a prediction rise from 6,943 to 7,827 (human complete.wLnc) and 5,423 to 5,921 (mouse), so precision falls: exact F1 −0.0074 (human complete.wLnc) and +0.0020 (mouse); on fragments.wLnc the sensitivity gain still dominates (+0.021, +0.033).
  - This supersedes the earlier `pytd_homology` runs at 27bfb5e (stop sensitivity +0.0004 to +0.0089 on complete inputs).
- **Runtime** (wall clock under shared load; athal complete): PyTD final 173 s, Perl final 411 s, Prodigal 45 s, GMST 34 s, CodAn 35,114 s, ANGEL 62,503 s.

## Checks
- **Rescoring:** `score_all.py` was re-run on 2026-10-21 over the current runs directory (output in the session scratchpad). It reproduces every row of the committed `summary.tsv` in columns 1–20. The only difference is 6 CodAn runs that finished after 10-20 15:37: athal and human `fragments.unstranded.wRand`, mouse `fragments.wRand`, and tthermophila `complete`, `complete_fullutr` and `fragments`.
- **Final-default reproducibility:** `run_tool.sh pytd_dev_v2|perl_td_dev_v2` was re-run at the committed trees (PY_dev 48c870c, TD_fix de3f57e, both clean) on spombe/complete and creinhardtii/complete. The resulting `preds.tsv` is byte-identical to the benchmark's for both tools.
- **Perl = Python at final defaults:** `preds.tsv` coordinates are identical in 84/84 inputs (E08).
- **Truth counts:** every input's `truth.tsv` row count equals the organism's kept-gene count in its truth log, and `.wRand` FASTAs have 2× the transcripts.
- **v6.0.0 vs v5.7.1:** `perl_td_v600` and `perl_td_v571` `preds.tsv` files are byte-identical in 80/80 runs, so Predict output did not change between these releases on these inputs.

- **Homology runs at one code version:** after the re-prediction, all 84 `pytd_dev_v2_homology` GFF3 files carry the 48c870c attribute `upstream_frame_contrast` and none carries the 44e417a attribute `upstream_coding_score`. In the 5 re-predicted runs, 7–30 ORFs (0.03–0.13%) differ from the 44e417a predictions they replaced (kept in `predict_44e417a/`).

## Gaps
- **Builder for `sprot_minus_<taxid>.dmnd`:** the per-organism Swiss-Prot taxid-exclusion databases used by `pytd_homology` have no builder in this repo. Only the later clade builder, `genome_annotation/scripts/make_sprot_clade_excluded.py`, is committed. The commands that made the per-taxid databases were not saved.
- **Move to `rerun_old16/`:** the move of the earlier final-default runs is not scripted (`launch16.sh` assumes an empty target). Nor are the deletions or moves before relaunching other labels.
- **ANGEL decision:** dropped as too inefficient. It was stopped on 2026-10-19 with Brian's agreement, after 3 completed athal runs (complete, complete_fullutr, complete_fullutr.unstranded; 17.4–21.9 h each, against minutes for the others). Three more athal runs had started and are partial. This decision is not yet written in `notes/` or `benchmark_v2/README.md`, where ANGEL is still listed as a tool.
- **CodAn scores:** 12 runs pending; `summary.tsv` must be regenerated when they finish.
- **Early benchmark observations in notes/04** (e.g. "TD 0.96–0.985 vs Prodigal 0.92–0.94"; human exact 14,204 → 13,099) come from earlier partial runs (bench_v2 truth, `$P/work/bench_v2`), not bench_v3. They should not be quoted.

## Caveats
- **Homology runs and the shared tree:** launch17 ran from `work/trees/PY_dev`, which moved to 44e417a (Kozak in bits) on 2026-10-22 while it was running. The Predict step of 5 *T. thermophila* runs (complete.wRand, complete.unstranded.wRand, fragments.unstranded, fragments.unstranded.wRand, fragments.wRand) used 44e417a. Their Predict step was re-run with `repredict_homology.sh` from a worktree pinned at 48c870c, reusing their DIAMOND and Pfam results. LongOrfs and the searches were unaffected, as 44e417a changes only start refinement. `runtime_sec.txt` of those runs is the original run's.
- **Run code predates its commits.** The final-default runs ended at 15:13 (PyTD) and 15:32 (Perl) on 2026-10-20, a few minutes before 48c870c (15:38:00) and de3f57e (15:37:58) were committed. The 2-input rerun above matches the commits, but the other 82 inputs per tool were not rerun.
- **Commits on remotes:** 48c870c is on PyTransDecoder origin/main and de3f57e on TransDecoder origin/devel (verified 2026-10-21 with git rev-list against the fetched remotes).
- **Labels track trees, not commits.** `perl_td` ran from `TD_fix` on 10-18 (351634a). Rerunning that label today would use de3f57e.
- **`fp` column bug:** `score_predictions.py` counts each non-duplicate FP twice in `fp` (`c[call] += 1` with `call == "FP"`, then `c["FP"] += 1`). For example, tthermophila CodAn has n_pred 3,665, tp_stop 229 and fp 6,872 = 2 × 3,436. `fp_random`, `fp_dup`, sensitivity, precision and F1 are unaffected (precision uses `n_pred`). **Fixed 2026-10-21** in `score_predictions.py`; a toy case (2 TP, 1 dup, 1 random FP, 2 other FPs) gives fp = 4 = n_pred − tp_stop. `summary.tsv` will be regenerated with the fixed scorer when CodAn finishes; until then do not quote `fp` from the committed table.
- **GeneMarkS-T on *T. thermophila*:** it supports only codes 1/4/11, so it was run with code 1; TAA/TAG (Gln in ciliates) are read as stops. The same applies to CPC2, Borf, RNAsamba and CodAn. Their *T. thermophila* scores reflect the code, not the method.
- **Runtimes** are wall clock with 6–16 concurrent jobs on a 28-core host (load ~40). They are indicative only.
- **Manuscript out of date:** `methods_benchmark.tex` still describes the superseded RefSeq `rna.gbff` truth (NM_/XM_ records, "longest CDS per gene") and four input types. It also omits complete_fullutr, CPAT/CPC2/Borf/RNAsamba and GeneMarkS-T (a `\todo`), and the GC values are `\todo` (22.1% and 68.3% here).
- **Why the random negatives matter (2026-10-23):** on stranded complete transcripts the longest-ORF baseline has a stop-level F1 of 0.954-0.984 (median 0.976, against 0.978 for the final defaults), so stop-level accuracy on complete coding transcripts alone does not separate methods. With random negatives added, its stop precision falls to a median of 0.803 (0.499 *T. thermophila*, 0.567 *C. reinhardtii*), against 0.993 for the final defaults. Its median exact F1 on complete inputs is 0.689. Source: `$P/work/bench_v3/scores/summary.interim_20261022.tsv` (these cells are final; only CodAn cells are pending). The manuscript (results_benchmark.tex) now presents the noncoding test before start codons for this reason.
- *C. reinhardtii* v5.5 transcripts all have annotated UTRs (0 added). v6.1 is only available through JGI login.
