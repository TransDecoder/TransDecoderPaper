# E28. Balanced mixed inputs: complete, 5' partial, 3' partial and internal CDSs (2026-10-25)

**Status:** all 225 runs are done and scored, including CodAn's 21 (finished 2026-10-27; `launch_mixed.out`: `CODAN_LANE_DONE`, 0 failed). CodAn (partial model) has stop / exact F1 0.957 / 0.788 overall; by type, exact F1 complete 0.598, 3' partial 0.590, 5' partial 0.971 (best of all methods), internal 0.960. Its trial range of exact F1 is at most 0.016. Rescoring added rows only: existing rows of `mixed_scores.tsv`, `partial_types.tsv` and `five_prime_by_distance.tsv` are unchanged. Table 1 and Table S9 now include CodAn; the Results text names CodAn as best on 5' partials.

## Claims supported
Manuscript (2026-10-25): Methods "Benchmark data" (balanced mixed inputs); Table 1 "mixed" columns;
Results "Balanced types of coding region" (two paragraphs and the 5'-open distance sentence);
Supplementary Table S9 (`tab:mixedtypes`); Discussion (start-codon paragraph, Limitations). Values are
from `data/benchmark/final_v7/mixed_medians.tsv` (`pipeline/derive/02_mixed_table.py`: median over
organisms within a trial, mean over trials). CodAn cells are marked as runs in progress. This record backs the notebook section "Balanced mixed inputs" in `notebooks/02_benchmark.Rmd`. It answers two questions:
- how the methods do when complete transcripts and all three partial types, including internal fragments (no start and no stop), are mixed in one input;
- how the decision to leave a 5' partial open depends on the distance to the first in-frame ATG.

Background:
- The `fragments` inputs (E07) contain only 5' and 3' partials, about half each, and no internal fragments.
- Their cut is at N(0.5, 0.1) of the CDS, so 5' partials cut close to the start codon are rare.
- Methods that train on their input could be skewed by an input of one partial type (Brian, 2026-10-25). These methods are TransDecoder (coding and start models), GeneMarkS-T and Prodigal.
- So mixed inputs were built in place of separate single-type inputs.

## Data
- **Source:** `$P/data/bench_v3/inputs/<org>/complete.{fa,truth.tsv}`, the same transcripts that `fragments` is cut from (UTRs of at most 250 nt; E07). All 8 organisms.
- **Inputs:** `$P/data/bench_v3/inputs/<org>/fragments_mixed_t{1,2,3}.{fa,truth.tsv}`. There are 24 inputs. They are stranded and have no negatives.
  - The names start with `fragments` on purpose. `run_tool.sh` keys on that prefix to choose CodAn's `*_partial` models and RNAsamba's partial-length weights, as for `fragments`.
- **Design** (`pipeline/benchmark/make_mixed_inputs.py`; the docstring has the full specification):
  - **Categories.** In trial *t* the transcripts are shuffled with `random.Random(t)` (seeds 1, 2, 3). They are then assigned in turn to complete, 5' partial, 3' partial and internal, so each category gets 25% (counts differ by at most 1).
  - **Cuts.** The same RNG draws the cuts. N is the number of codons, including the ATG and the stop.
    - 5' partial: the cut falls in codon k, drawn uniformly from 2 to min(floor(0.8N), N−50), at a random phase. The 3' piece is kept: the rest of the CDS plus the 3' UTR.
    - 3' partial: the last whole codon kept, m, is drawn uniformly from max(ceil(0.2N), 50) to N−1. Then 0–2 extra bases are kept. The 5' UTR plus the 5' piece are kept.
    - internal: the window is a fraction f of the CDS, with f ~ N(0.5, 0.1) clipped to [0.05, 0.95]. The window length is clipped to [152, N·3−6] nt and placed uniformly at random strictly between the start and stop codons. No UTR remains, and 152 nt hold at least 150 nt of whole codons in any phase.
    - At least 150 nt of CDS (whole codons) must remain. A transcript that cannot meet its category is dropped.
  - **Truth.**
    - cds_start is the first whole codon of the original frame. cds_end is the end of the stop codon (complete, 5' partial) or of the last whole codon (3' partial, internal).
    - `type` is one of complete / 5prime_partial / 3prime_partial / internal.
    - Extra columns:
      - `cut_codon_5p` and `cut_codon_3p`: codons removed wholly or in part at each end; the stop counts at the 3' end.
      - `dist_first_inframe_atg`: nt from the fragment start to the first in-frame ATG within the reference CDS. It is set for 5' partials and internals and is NA otherwise.
- **Counts** (references per input / dropped, by category: complete, 5', 3', internal):

  | org | source | t1 refs | t2 refs | t3 refs | dropped t1 (5'/3'/int) | t2 | t3 |
  |---|---|---|---|---|---|---|---|
  | athal | 26,867 | 26,787 | 26,764 | 26,777 | 22/32/26 | 40/30/33 | 32/14/44 |
  | creinhardtii | 17,621 | 17,621 | 17,621 | 17,621 | 0/0/0 | 0/0/0 | 0/0/0 |
  | dmel | 13,884 | 13,789 | 13,789 | 13,778 | 38/26/31 | 28/32/35 | 45/31/30 |
  | human | 19,167 | 19,127 | 19,137 | 19,131 | 14/8/18 | 11/8/11 | 13/7/16 |
  | mouse | 21,409 | 21,367 | 21,357 | 21,362 | 15/11/16 | 13/22/17 | 18/13/16 |
  | pfalciparum | 5,281 | 5,275 | 5,277 | 5,278 | 1/1/4 | 0/1/3 | 0/1/2 |
  | spombe | 5,100 | 5,086 | 5,084 | 5,084 | 5/4/5 | 5/5/6 | 5/6/5 |
  | tthermophila | 26,456 | 26,377 | 26,362 | 26,372 | 26/22/31 | 26/33/35 | 33/30/21 |

  No complete transcript is ever dropped. The dropped ones have a CDS of fewer than about 52 codons. md5 (first 8 hex digits) of the `truth.tsv` files, t1/t2/t3:
  - athal de1e2a1b/f9a1c3b7/ef4681a5
  - creinhardtii 8bffea42/33739c98/aa30344d
  - dmel 619caa59/e6e5a1d3/eedf1440
  - human b92882ad/d9bf42b9/5cc0b2d3
  - mouse d4cfbd90/d2eab4d5/ef0ca77a
  - pfalciparum 86be074a/88ef8ab2/97d51a34
  - spombe 7a7c38f0/fc7bf258/413eaa47
  - tthermophila e2a79925/3258f87d/7f92df24
- **Validation** (scratch script, not committed): each truth CDS was checked against the reference transcript. Every input had 0 failures on each of these checks:
  - the CDS is the reference CDS restricted to the fragment, in frame;
  - start and stop are present exactly as the type implies;
  - there is no internal stop under the organism's table;
  - there are at least 150 nt of whole codons;
  - the cut columns and the first in-frame ATG are consistent.

  The only flags were 5 *C. reinhardtii* internals inside tandem repeats, where the check's `str.find` located an earlier copy of the fragment.

## Software
- The tools and settings are those of the `fragments` input (E07, E25), via `benchmark_v2/scripts/run_tool.sh`, which was not modified.
  - `pytd_v7`: PyTD 3a793ff. `perl_td_v7`: Perl 41fc5ee. `perl_td_v600`: v6.0.0.
  - prodigal, gmst, codan (`*_partial` model, confirmed in the process list, e.g. `PLANTS_partial -s plus`), cpat, cpc2, rnasamba (partial weights) and borf.
  - longest ORF: `pytd_v7/longest_orf.preds.tsv`. This is the pytd_v7 LongOrfs output, not the old `pytd` run that the `fragments` baseline uses.
- GeneMarkS-T, CodAn, CPC2, RNAsamba and Borf were not run on *T. thermophila*, which uses the ciliate code (as in E07).
- Run directories: `$P/work/bench_v3/runs/<org>/fragments_mixed_t<k>/<tool>/`. Logs: `runs/logs/<org>.fragments_mixed_t<k>.<tool>.log`.

## Commands
```bash
cd $P
env/td/bin/python TransDecoderPaper/pipeline/benchmark/make_mixed_inputs.py data/bench_v3/inputs     # seeds 1,2,3
nohup setsid bash work/bench_v3/launch_mixed.sh > work/bench_v3/runs/logs/launch_mixed.out 2>&1 &
#   two lanes, 4 jobs max: CodAn (xargs -P 2, longest inputs first) and the other tools (xargs -P 2)
env/td/bin/python TransDecoderPaper/pipeline/benchmark/score_mixed.py data/bench_v3/inputs work/bench_v3/runs
#   -> data/benchmark/final_v7/mixed_scores.tsv, five_prime_by_distance.tsv
env/td/bin/python TransDecoderPaper/pipeline/derive/02_partial_types.py data/bench_v3/inputs work/bench_v3/runs \
  > TransDecoderPaper/data/benchmark/final_v7/partial_types.tsv      # now also the 3 mixed variants
TransDecoderPaper/scripts/knit_notebooks.sh 02_benchmark.Rmd
```
**Rescoring when CodAn finishes.** Check that `launch_mixed.out` shows `CODAN_LANE_DONE` and 21 `codan` lines. Then rerun the last three commands; CodAn rows are picked up automatically. In the rows of `partial_types.tsv` for the existing `fragments*` variants, the regenerated file was byte-identical to the previous one.

## Outputs
- `data/benchmark/final_v7/mixed_scores.tsv`: one row per org × trial × tool (228 without CodAn). It has the columns of `summary.interim_20261025.tsv`, including `runtime_sec`.
- `data/benchmark/final_v7/partial_types.tsv`: per type, now including `fragments_mixed_t*` (types complete, 5prime_partial, 3prime_partial, internal).
- `data/benchmark/final_v7/five_prime_by_distance.tsv`: 5' partials per org × trial × tool × stratum × bin. The strata are `dist_first_inframe_atg` (0-29, 30-59, 60-119, 120-299, 300+, none) and `cut_frac` (0.1-wide bins).

## Results
Values are the median over organisms (7 for the n/a tools) in each trial, given as the mean over the 3 trials, with the min–max over trials in parentheses. They come from notebook 02, section "Balanced mixed inputs". No method or measure has a range over trials above 0.008.

| method | stop sens | stop prec | stop F1 | exact sens | exact prec | exact F1 |
|---|---|---|---|---|---|---|
| TransDecoder (pytd_v7) | 0.944 (0.942–0.945) | 0.997 | 0.970 (0.969–0.971) | 0.818 | 0.870 | 0.843 (0.841–0.844) |
| TransDecoder v6.0.0 | 0.853 | 0.998 | 0.916 | 0.712 | 0.832 | 0.768 |
| longest ORF | 0.979 | 0.980 | 0.979 | 0.839 | 0.839 | 0.839 |
| Prodigal | 0.971 | 0.941 | 0.952 | 0.815 | 0.787 | 0.796 |
| GeneMarkS-T (7) | 0.852 | 0.996 | 0.918 | 0.747 | 0.865 | 0.809 |
| CPAT | 0.768 | 0.994 | 0.867 | 0.428 | 0.556 | 0.493 |
| RNAsamba (7) | 0.887 | 0.959 | 0.920 | 0.485 | 0.523 | 0.496 |
| CPC2 (7) | 0.573 | 0.990 | 0.726 | 0.326 | 0.582 | 0.418 |
| Borf (7) | 0.823 | 0.989 | 0.890 | 0.571 | 0.695 | 0.623 |
| CodAn (7) | 0.921 | 0.994 | 0.957 | 0.767 | 0.804 | 0.788 |

**Sensitivity by type**, given as stop / exact (mean over trials of the median over organisms):

| method | complete | 5' partial | 3' partial | internal |
|---|---|---|---|---|
| TransDecoder | 0.959 / 0.893 | 0.936 / 0.776 | 0.936 / 0.869 | 0.940 / 0.783 |
| TransDecoder v6.0.0 | 0.946 / 0.756 | 0.819 / 0.721 | 0.847 / 0.638 | 0.807 / 0.735 |
| longest ORF | 0.988 / 0.701 | 0.953 / 0.953 | 0.984 / 0.686 | 1.000 / 1.000 |
| Prodigal | 0.977 / 0.666 | 0.963 / 0.961 | 0.971 / 0.660 | 0.976 / 0.972 |
| GeneMarkS-T | 0.956 / 0.771 | 0.820 / 0.784 | 0.849 / 0.660 | 0.798 / 0.772 |
| CPAT | 0.935 / 0.909 | 0.663 / 0.015 | 0.830 / 0.782 | 0.638 / 0.015 |
| RNAsamba | 0.981 / 0.952 | 0.778 / 0.018 | 0.952 / 0.929 | 0.827 / 0.018 |
| CPC2 | 0.934 / 0.915 | 0.660 / 0.014 | 0.400 / 0.380 | 0.253 / 0.006 |
| Borf | 0.957 / 0.894 | 0.733 / 0.318 | 0.851 / 0.763 | 0.774 / 0.349 |
| CodAn | 0.957 / 0.594 | 0.957 / 0.952 | 0.866 / 0.578 | 0.931 / 0.928 |

**5' partials, exact sensitivity (ORF left open), by distance from the fragment start to the first in-frame ATG:**

| method | 0-29 | 30-59 | 60-119 | 120-299 | 300+ | no ATG |
|---|---|---|---|---|---|---|
| TransDecoder | 0.668 | 0.687 | 0.762 | 0.836 | 0.913 | 0.725 |
| TransDecoder v6.0.0 | 0.698 | 0.651 | 0.689 | 0.808 | 0.974 | 0.379 |
| Prodigal | 0.973 | 0.968 | 0.971 | 0.974 | 0.982 | 0.829 |
| GeneMarkS-T | 0.687 | 0.773 | 0.807 | 0.869 | 0.997 | 0.395 |
| longest ORF | 0.964 | 0.959 | 0.967 | 0.973 | 0.990 | 0.756 |

Stop-level sensitivity of TransDecoder in the same bins is 0.931 / 0.945 / 0.943 / 0.958 / 0.994 / 0.725. So most 5' partials that are not left open still have the right stop; the start is placed at a downstream ATG. By the fraction of the CDS removed, TransDecoder's exact sensitivity is flat (0.74–0.79). Over the same range GeneMarkS-T falls from 0.86 to 0.58 and v6.0.0 from 0.76 to 0.56.

**Observations**
- **Which methods can represent internal partials.** TransDecoder, Prodigal, GeneMarkS-T and (by design) CodAn report ORFs that are open at both ends. The longest-ORF baseline reaches 1.000 on internals by construction: the reference frame of an internal fragment has no stop, so it is the longest ORF.
  - The ATG-requiring classifiers (CPAT, CPC2, RNAsamba) essentially never give the open 5' end (exact ≤ 0.018). Their stop-level hits on internals and 5' partials start at an internal ATG.
  - Borf reports some 5'-open ORFs (exact 0.32–0.35).
- **TransDecoder on internals and 5' partials** has high stop-level sensitivity (0.94). It starts at an in-frame ATG in about 16% of these fragments where the reference is open. This happens mostly when an ATG lies near the fragment start (exact 0.67 at 0-29 nt vs 0.91 at 300+ nt). E29 studies this decision; this record only measures it.
- **Mixed vs single-type inputs.** On `fragments` (E07; single type per sequence, cut at about half the CDS) TransDecoder's medians are stop / exact 0.871 / 0.735 on 5' partials and 0.876 / 0.816 on 3' partials. In the mixed inputs they are 0.936 / 0.776 and 0.936 / 0.869. The fragment-length distributions differ (mixed 3' partials keep 20–100% of the CDS), so these are not a controlled comparison of training composition.
- **Run time** (mean per run, shared host): pytd_v7 3.3 min, perl_td_v7 8.7, v6.0.0 5.2, all others ≤ 1 min.

## Checks
- `02_partial_types.py` summed over types reproduces `tp_stop` and `tp_exact` of `mixed_scores.tsv` (same matcher, one credit per reference).
- After the knit, the notebook's checks section shows 0 of 133 values differing from the manuscript.

## Caveats
- CodAn results are pending (see Status).
- Superseded inputs from earlier in the same session remain on disk. They are not part of this record and are not scored. `fragments_internal` and `fragments_5p` were made before the plan changed to mixed inputs. Their files are `$P/data/bench_v3/inputs/<org>/fragments_internal.*` and `fragments_5p.*`.
  - Run directories: `fragments_internal/pytd_v7` for 6 organisms (finished), and killed CodAn directories for athal on both inputs (no `preds.tsv`).
  - Their generator scripts were moved out of the repo to the session scratchpad.
  - `score_all.py` walks every run directory, so the next full rescoring would add `fragments_internal` rows. A `run_all.sh` call without `VARIANTS` would also run tools on these inputs, and on the new mixed inputs, which changes `N_INPUTS` in notebook 02 if they enter `summary.tsv`. Either move these inputs and runs aside, or filter them out.
- The longest-ORF baseline here comes from pytd_v7 LongOrfs, not from the older `pytd` run used for the other inputs.

## Housekeeping (coordinator, 2026-10-25)
- Superseded `fragments_internal` / `fragments_5p` inputs and runs moved to
  `$P/_trash/superseded_internal_20261025/`.
- `score_all.py` walks every run directory: before the next full rescore, filter out
  `fragments_mixed_t*` (scored separately by `pipeline/benchmark/score_mixed.py`) so that the
  84-input summary and `N_INPUTS` in notebook 02 are unchanged.
