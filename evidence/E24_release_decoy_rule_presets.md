# E24. Release implementation: decoy-FDR rule and bundled presets in Perl and PyTD, with identical output (2026-10-24)

## Claims supported
- The decoy-FDR ORF rule (E20) with named q-value settings (E22) and presets (E21) is the release ORF acceptance rule in both TransDecoder (Perl, `devel` 41fc5ee) and PyTransDecoder (`main` 3a793ff).
  - Without homology support, an ORF needs in-frame score > 0, frame contrast > 0 and q ≤ `--max-qvalue`. The default is 0.05; `--precise` is 0.01 and `--sensitive` is 0.2.
  - ORFs at or above the GC-table length are kept on length alone.
  - The 100-codon short-ORF contrast test and the experimental rule shapes (E19) are removed.
- **Perl and PyTD give identical results.** Their decoys, q-values and predicted coordinates are byte-identical on all 7 parity cases (`data/benchmark/release/parity_perl_pytd.tsv`). Together they cover:
  - complete, fragment and lncRNA inputs;
  - stranded and unstranded inputs;
  - the standard and Tetrahymena genetic codes;
  - the default, `--precise` and `--sensitive` settings;
  - preset mode;
  - a small input.
  
  | case | q-values | predictions |
  |---|---|---|
  | human complete + lncRNA | 64,936 | 20,148 |
  | T. thermophila fragments + random (Tetrahymena code) | 32,214 | 21,424 |
  | D. melanogaster fragments, unstranded + random | 22,854 | 12,079 |
  | C. reinhardtii complete + random, `--precise` | 18,434 | 16,811 |
  | human, 200 sequences, no preset | 353 | 154 |
  | human, 200 sequences, `--preset` | — | 97 |
  | human, 200 fragments, `--preset --sensitive` | — | 95 |

  Two more checks also matched:
  - an S. pombe complete + random run (scratch), with identical decoy q-values and 4,984 identical predictions;
  - all 42,700 preset threshold lookups (7 bundled presets × 5 q values × 1,220 strata, including the nearest-stratum fallbacks).
  
  Other output files differ only in formatting that was already there (FASTA line width, float digits in `hexamer.scores`, extra GFF attributes written by PyTD).
- **Seven bundled presets** ship in both packages (`pytransdecoder/presets`, `PerlLib/presets`; about 0.8 MB each, gzip):

  | name | organism | genetic code | GC % | auto-accept (nt) |
  |---|---|---|---|---|
  | mammal | human (+ lncRNA) | Universal | 49.8 | 726 |
  | insect | D. melanogaster | Universal | 50.8 | 804 |
  | plant | A. thaliana | Universal | 42.1 | 663 |
  | fungus | S. pombe | Universal | 37.8 | 615 |
  | apicomplexan | P. falciparum | Universal | 22.0 | 519 |
  | ciliate | T. thermophila | Tetrahymena | 25.9 | none |
  | green_alga | C. reinhardtii | Universal | 68.2 | 1347 |

  Each is built from all complete coding transcripts plus their coding fragments (and GENCODE lncRNAs for human), following E21's finding that presets need fragments.
- **Smoke test** (`bundled_preset_smoke.tsv`): for each organism, 200-sequence inputs of 100 coding + 100 random.
  - With the bundled preset: 0–1 of 100 random sequences get calls.
  - Trained on the input: 2–70 of 100 do (median 51.5).
  - Stop sensitivity changes by +1 to −10 points with the preset.
  
  These samples overlap the preset references, so this is a functional check only. The held-out evaluation is E21.
- Both implementations warn below 1,000 input transcripts when no preset is given.

## Implementation
**Shared design**
- Portable random numbers: xoshiro128** on 32-bit words, seeded by FNV-1a of `"<seed>:<transcript>:<copy>#<i>"`. Perl uses `PerlLib/Decoys.pm`; PyTD uses `core/prng.py`.
- Altschul–Erickson dinucleotide shuffle with Wilson's last-exit tree. Edges are visited in sorted order; Fisher–Yates shuffle.
- Strata use exact integer arithmetic. ORF-length bins come from powers by repeated multiplication; transcript-length bins are `bit_length(len // 500)`.
- Selection uses q-values at the precision written to `decoy_qvalues.tsv` (`%.5g`).
- Golden values are shared: `tests/test_decoys.py` and `__testing/test_decoys.pl`. A bulk check (2,000 random sequences including N, bins for 1–20,000 codons, 3,000 q-values) matched byte for byte.

**Perl**
- `util/decoy_qvalues.pl` writes the shuffled transcripts, runs `TransDecoder.LongOrfs` on them with the settings in `longorfs.params.json` (now written by LongOrfs), and scores the decoy ORFs with `score_CDS_likelihood_all_6_frames.pl`. The decoy ORFs are therefore exactly what LongOrfs finds.
- In preset mode it applies the thresholds instead (`preset_qvalues.tsv`).
- `select_best_ORFs_per_transcript.pl` takes `--qvalues` and `--max_qvalue`.
- `TransDecoder.Predict` adds `--max_qvalue/--decoy_q`, `--precise`, `--sensitive`, `--decoy_shuffles`, `--decoy_seed` and `--preset`, which supplies the training ORFs, hexamer model, Kozak model, auto-accept length and genetic-code check. The `TransDecoder` wrapper passes these through.

**Runtime of the decoy step**, human + lncRNA (53k transcripts): about 11 min in Perl and 3 min in PyTD, against 7 min to score the candidate ORFs.

**Tests**
- PyTD: 67 tests pass.
- Perl `__testing`: 5 targets pass. The `single_best` fixtures were refreshed in 20a0f88 because they predated the coding-score-first ranking (f7e75a1); the committed script already failed them.
- The Perl wrapper ran end to end on `sample_data/simple_transcriptome_target` (921 Trinity transcripts):
  - the small-input warning was printed;
  - automatic shuffles gave 18 per transcript;
  - genome mapping worked;
  - `--preset mammal` ran.

**Commits**
- PyTD `main` 44e417a..3a793ff includes the experimental commits of E18–E22, then:
  - be23241: the rule becomes the default and the experimental rules are removed;
  - c78845d: q rounding;
  - 3a793ff: presets.
- Perl `devel` ab38ca6..41fc5ee: 20a0f88, 4342436 and 41fc5ee.
- Release worktrees for the final runs: `work/trees/PY_v7` (3a793ff) and `work/trees/TD_v7` (41fc5ee).

**Bundled preset builds**
- `make_bundled_preset_refs.py`, then `data/benchmark/release/bundled_preset_run_ref.sh`: a fresh PyTD LongOrfs -S, Predict and build-preset.
- The PyTD working tree used had the content of c78845d; `pytd_version.txt` says be23241 because the q-rounding change was not yet committed.

## Caveats
- The bundled presets are built from the benchmark annotations, so they must not be scored on the benchmark inputs. Preset accuracy claims come from the held-out E21 test.
- The human preset parity cases used the E21 `humanmix` preset directory, built with the earlier random-number generator. That tests the preset code path, not the bundled files.
- Final benchmark and GM12878 reruns at 3a793ff / 41fc5ee are in progress (`run_tool.sh` tools `pytd_v7*` and `perl_td_v7*`; `genome_annotation/scripts/gm12878_final_v7.sh`); they will be recorded in E25.
