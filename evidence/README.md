# Evidence records

Each claim in the manuscript traces to one record here. A record states what was done, with what
data and software, using which exact commands, and where the outputs are. Together with the scripts
in this repository, the records are meant to serve as the paper's supplementary code and methods
detail.

## Record format
One Markdown file per analysis, named `E<NN>_<topic>.md`, with these sections:
1. **Claims supported**: the claims, quoted or paraphrased, and the manuscript section(s) that make them.
2. **Data**: accessions and URLs, download date, local paths, checksums or read counts where useful.
3. **Software**: tool versions, and repository commits for PyTransDecoder, TransDecoder and this repo.
4. **Commands**: the exact commands, in order. Use scripts in this repo where they exist, and record
   one-off commands verbatim.
5. **Outputs**: committed result tables (`*/results/`) and large local outputs (paths under the
   project work area, which are not in git).
6. **Results**: the numbers the manuscript uses.
7. **Checks**: reproduction runs, sanity checks and controls.
8. **Caveats**.

Rules:
- Every number in the manuscript must appear in a committed result table or in a record.
- Manuscript source lines carry `% evidence: evidence/E<NN>_<topic>.md`.
- Paths: `$P` = project root (`/broad/hptmp/bhaas/TRANSDECODER_PAPER` on the analysis machine);
  `$T` = this repository (`$P/TransDecoderPaper`).
- Notebooks (`*.ipynb`) are used where a figure or interactive exploration helps. They are committed
  with outputs, next to the record that cites them.

## Index

| Record | Topic | Manuscript | Status |
|---|---|---|---|
| [E01](E01_gm12878_riboseq.md) | GM12878 Ribo-seq support of predicted ORFs | Results: GM12878 evidence; Methods: Ribo-seq support | complete |
| [E02](E02_gm12878_ms_peptides.md) | GM12878 MS peptide support of predicted ORFs | Results: GM12878 evidence; Methods: MS support | complete |
| E03 | GM12878 three-assembly annotation vs MANE/GENCODE (final defaults) | Results: application | to write (notes/06) |
| E04 | Minimum ORF length sweep (benchmark lncRNA/random negatives; GM12878 post-hoc) | Results: short ORFs | to write (notes/05, notes/06) |
| E05 | Expressed-gene sensitivity (GM12878) | Results: application | to write (notes/06) |
| [E06](E06_gm12878_homology_series.md) | Homology restriction series (clade-excluded Swiss-Prot); Alu-in-UTY_PANTR retention artifact | Results: homology | complete |
| E07 | Cross-species benchmark (bench_v3): data, truth sets, tools, scoring | Methods/Results: benchmark | to write |
| E08 | Perl/Python concordance and runtime | Results | to write (regression/) |
| E09 | Start-codon refinement and Kozak | Methods/Results | to write (analyses/start_*, figures/kozak) |
| E10 | Long-ORF retention rule and GC | Methods/Results | to write (analyses/long_orf_threshold) |
| E11 | Trinity genome-guided assembly of GM12878 | Methods: application data | to write |

## Large files (decision, Brian 2026-10-21)
- Until the paper's analyses are final, large outputs stay in the project work area (`$P/work`, `$P/data`). Each record lists their paths.
- Once the analyses are final, the important large files go into this repo via Git LFS: benchmark inputs and truth sets, the GM12878 assemblies and genome-projected predictions, and the per-ORF tables.
- A subset chosen as the paper's supplementary material is also deposited separately on Zenodo, for a DOI.
