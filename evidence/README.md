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
| [E03](E03_gm12878_three_assemblies.md) | GM12878 three-assembly annotation vs MANE/GENCODE (final defaults) | Results: application | complete |
| [E04](E04_min_orf_length.md) | Minimum ORF length sweep (benchmark lncRNA/random negatives; GM12878 post-hoc) | Results: short ORFs | complete |
| [E05](E05_expressed_sensitivity.md) | Expressed-gene sensitivity (GM12878) | Results: application | complete |
| [E06](E06_gm12878_homology_series.md) | Homology restriction series (clade-excluded Swiss-Prot); Alu-in-UTY_PANTR retention artifact | Results: homology | complete |
| [E07](E07_benchmark_v3.md) | Cross-species benchmark (bench_v3): data, truth sets, tools, scoring | Methods/Results: benchmark | complete except CodAn scores (71/80 runs); homology at final defaults added 2026-10-22 |
| [E08](E08_perl_python_concordance.md) | Perl/Python concordance and runtime | Methods: implementations | complete (rerun at final commits 2026-10-21 and after Kozak-in-bits 2026-10-22) |
| [E09](E09_start_refinement_kozak.md) | Start-codon refinement and Kozak | Methods/Results | complete (gaps listed) |
| [E10](E10_long_orf_rule_gc.md) | Long-ORF retention rule and GC | Methods/Results | complete (gaps listed) |
| [E11](E11_gm12878_trinity_gg.md) | Trinity genome-guided assembly of GM12878 | Methods: application data | complete |
| [E12](E12_long_orf_failures.md) | True CDSs that fail the coding criteria; long-ORF rule rescues; causes (composition, low complexity, unconserved models); GM12878 MANE | Methods: long-ORF rule (E10 gap); Results: coding-test failures | complete |
| [E13](E13_start_param_sensitivity.md) | Start-refinement parameter sensitivity: alpha x Kozak weight grid, leave-one-organism-out, GM12878 sub-grid | Methods: start-codon refinement | complete |
| [E14](E14_references.md) | Literature references added via paperclip/Crossref: TGD and genome, *Drosophila* chr4/Y codon bias, Kozak, assembly and alignment tools | Methods, Results: coding-test failures | complete (Singh2014 not read in full) |

## Large files (decision, Brian 2026-10-21)
- Until the paper's analyses are final, large outputs stay in the project work area (`$P/work`, `$P/data`). Each record lists their paths.
- Once the analyses are final, the important large files go into this repo via Git LFS: benchmark inputs and truth sets, the GM12878 assemblies and genome-projected predictions, and the per-ORF tables.
- A subset chosen as the paper's supplementary material is also deposited separately on Zenodo, for a DOI.
