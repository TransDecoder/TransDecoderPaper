# TransDecoder paper

Working repository for the TransDecoder manuscript: algorithm description, application examples, and accuracy benchmarking.

- `notes/` — technical notes reconstructed from the code, plus inventories of existing assets
- `manuscript/` — LaTeX source: git submodule [TransDecoderPaper_LaTeX](https://github.com/TransDecoder/TransDecoderPaper_LaTeX) (Overleaf-synced). Clone with `git clone --recurse-submodules`.
- `notes/EXPERIMENT_LOG.md` — internal log of every experiment, idea and design choice tried (adopted / not adopted / open), with key numbers, lessons and pointers to the evidence records
- `evidence/` — one record per analysis supporting the manuscript's claims: data, software versions, exact commands, outputs, results, checks (see `evidence/README.md`)
- `analyses/`, `genome_annotation/`, `benchmark_v2/` — analysis scripts and committed result tables
- `benchmark/` — benchmarking pipeline, scripts and notebooks (to come)
- `figures/` — figure-generating code and outputs (to come)
