# E21. Presets for small or partial inputs (2026-10-24)

## Claims supported (planned; Brian: presets are in this paper)
- **Training on the input fails when the input is small and noisy.** With 50–500 sequences (half coding, half random), both the release rule and the decoy rule call ORFs on 8–82% of random sequences. The training set (top 500 longest ORFs) is then partly chance ORFs.
- **Clean small inputs are less affected.** All-coding sets of 50–200 transcripts keep high precision (0.96–1.0), but exact-start sensitivity drops: complete, 50 transcripts, 0.80–0.81 trained on the input vs 0.90 with a same-species preset.
- **The threshold sits between 500 and 1,000 sequences.** At 1,000 sequences (half random), training on the input works again: 0.8–0.9% random calls on complete transcripts, 3.7–3.8% on fragments.
- **A preset removes the small-input failure.** It is built from a decoy-rule run on a large held-out reference of the same species. With it, random-sequence calls are 0–0.4% at every size, at sensitivity close to training on a large input.
- **The preset reference should include fragments.** Presets built from complete transcripts only lose sensitivity on fragments: 0.806 vs 0.880 trained on the input, at 5,000 sequences. Building from complete transcripts plus fragments ("mix") closes the gap (0.873) without hurting complete transcripts (0.956 vs 0.958).
- **Presets need to be clade-matched.** Mouse ↔ human transfers well. Human → Arabidopsis or Drosophila loses a lot: fragments 0.72 vs 0.87 (mix, same species).

## Design (PyTD exp/short-orf-taper, 6b024ca; local only)
- `pytransdecoder build-preset -t ref.fa --run-dir <decoy-rule run> --name N -o DIR` writes:
  - `preset.json`: genetic code, GC %, GC-table auto-accept length, strata, q levels;
  - `hexamer.scores`;
  - `training_orfs.cds`, which fixes mu;
  - `kozak.tsv`;
  - `thresholds.tsv`: for each stratum and q level (0.001–0.5), the minimum frame contrast accepted in the reference's q-values.
- `predict --preset DIR --decoy-q Q`:
  - uses the preset's training ORFs and coding model (no training on the input) and its start model;
  - accepts eligible ORFs whose contrast is at least the stratum threshold at the largest stored level ≤ Q;
  - a stratum absent from the reference takes the nearest ORF-length bin with the same transcript-length bin and type;
  - long-ORF auto-accept uses the preset's GC-table length.

## Test (held out)
- **Split.** Coding transcripts of human, mouse, A. thaliana and D. melanogaster are split by an MD5 hash of the ID. Random and lnc sequences are hashed by their base ID.
  - **Reference A:** complete transcripts, plus the A-half lncRNAs for human and mouse. "mix" adds the A-half fragments (`setup_presets.py` and the inline script in the session; refs in `work/bench_v3/presets/ref/`).
  - **Tests from half B:**
    - coding + random 1:1, complete and fragments, at 50, 200, 500 and 1,000 (3 replicates each) and 5,000 (1 replicate);
    - coding-only sets of 50, 200 and 1,000 (3 replicates).
- **Methods compared:**
  - trained on the input, with the release rule (cliff_c200) or the decoy rule (decoy_q0.05_longgc);
  - same-species preset (complete-only, or mix);
  - other-species preset: mouse for human, human for mouse, Arabidopsis and Drosophila.
- **Runs.** 912 fresh runs, all at 6b024ca. Scores: `data/benchmark/presets/preset_scores.tsv` (score_negatives.py). Preset metadata and thresholds: `data/benchmark/presets/preset_meta/`.

## Results (median over organisms × replicates; stop sensitivity / % random with a prediction)
| complete + random | trained on input, release | trained on input, decoy | preset same | preset same mix | preset other |
|---|---|---|---|---|---|
| 50 | 1.000 / 68 | 1.000 / 82 | 0.980 / 0 | 0.960 / 0 | 0.940 / 0 |
| 200 | 0.990 / 64 | 0.990 / 53 | 0.960 / 0 | 0.960 / 0 | 0.915 / 0 |
| 500 | 0.984 / 15 | 0.986 / 7.6 | 0.952 / 0.4 | 0.948 / 0 | 0.900 / 0.4 |
| 1,000 | 0.967 / 0.9 | 0.979 / 0.8 | 0.956 / 0.4 | 0.958 / 0.2 | 0.913 / 0.3 |
| 5,000 | 0.963 / 0.5 | 0.967 / 0.6 | 0.958 / 0.3 | 0.956 / 0.2 | 0.914 / 0.3 |

| fragments + random | trained on input, release | trained on input, decoy | preset same | preset same mix | preset other |
|---|---|---|---|---|---|
| 50 | 0.880 / 52 | 0.920 / 62 | 0.840 / 0 | 0.860 / 0 | 0.740 / 0 |
| 200 | 0.875 / 45 | 0.920 / 57 | 0.815 / 0 | 0.875 / 0 | 0.740 / 0 |
| 500 | 0.884 / 24 | 0.898 / 25 | 0.804 / 0.4 | 0.866 / 0.2 | 0.718 / 0 |
| 1,000 | 0.825 / 3.8 | 0.868 / 3.7 | 0.796 / 0.3 | 0.854 / 0.2 | 0.712 / 0.2 |
| 5,000 | 0.816 / 0.1 | 0.880 / 0.3 | 0.806 / 0.2 | 0.873 / 0.3 | 0.723 / 0.2 |

Coding-only sets (stop sensitivity / exact sensitivity / stop precision):

| set | trained on input, decoy | preset same mix |
|---|---|---|
| complete, 50 | 0.980 / 0.810 / 1.000 | 0.960 / 0.880 / 1.000 |
| fragments, 50 | 0.910 / 0.730 / 0.958 | 0.840 / 0.740 / 1.000 |
| fragments, 1,000 | 0.869 / 0.770 / 0.997 | 0.859 / 0.761 / 0.995 |

## Implications
- Build bundled presets from complete transcripts plus fragments, one per clade.
- Document `build-preset` for users with a related-species transcriptome.
- Warn below about 1,000 transcripts. The exact criterion could be based on how many usable training ORFs there are.
- The earlier subsets (`work/bench_v3/subsets`, drawn from the full sets, not held out) are superseded by this held-out test.
