# Reproduction detail A: GM12878 application (E01–E03, E05, E06, E11) and simulated Trinity assemblies (E30–E32)

> **Paths into `$P/work`.** This section was written on 2026-11-03, before `$P/work` was deleted. Files that were kept are now under `$P/data/work_keep/` at the same path relative to `work/`. For example, `$P/work/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout` is now `$P/data/work_keep/bench_v3/runs/human/complete/pytd_v7_homology/pfam.domtblout`. Everything else under `$P/work` has to be regenerated with the steps given here. The tables in Part 2 record what was kept and why; see `../REPRODUCE.md` for the full list.

`P=/broad/hptmp/bhaas/TRANSDECODER_PAPER`, `T=$P/TransDecoderPaper`. All scripts hardcode `P` (edit the
`P=` line, or keep the same layout). Environments are `$P/env/<name>`, specs in `environment/conda_envs/<name>.yml`
(`td`, `rnaseq_asm`, `trinity215`, `polyester`, `riboseq`, `sratools`, `proteo`, `bedtools`).

## 0. Code versions and how to plug in a pip-installed PyTD 0.4.0

| label in runs | commit | same output as 0.4.0? | used by |
|---|---|---|---|
| `pytd_v7`, `final_v7_*` | PyTD 3a793ff (`work/trees/PY_v7`) | yes for ab initio and retention-only (E40: byte-identical to d3892b5) | E25 GM12878 release runs, E30-E32 |
| `hseries_v7h`, `final_v7h_*` | PyTD a5b24c3 (`work/trees/PY_hstart`) | yes (same `pytransdecoder/` tree as d3892b5) | E06 release series |
| `rel` (release_pin) | PyTD d3892b5 (`work/trees/PY_rel_d3892b5`) | yes (pin) | E40 reruns |
| `final_abinitio`, `hseries` | PyTD 48c870c (`work/trees/PY_48c870c`) | **no** (development; fixed short-ORF test) | E01/E02/E03/E05/E06 development tables (kept for comparison only) |
| `perl_td_v600` | TransDecoder v6.0.0 66d4712 (`work/trees/TD_base`) | n/a (comparator) | E30-E32 |

Worktrees are `git worktree`s of `$P/repos/{PyTransDecoder,TransDecoder}` (tags v0.2.0 ... v0.4.0 exist there).
All option names the scripts use (`--retain-blastp-hits`, `--retain-pfam-hits`, `--no-homology-start-assist`,
`--start-contrast-alpha`, `--start-kozak-weight`, `--max-qvalue`, `--no-refine-starts`, `--precise`, `--sensitive`,
`genome-map`) exist in v0.4.0.

Every TD call is `PYTHONPATH=<tree> $P/env/td/bin/python -m pytransdecoder.cli ...`, preceded by
`(cd <tree> && git rev-parse --short HEAD) > pytd_version.txt` under `set -e`. So a bare `pip install
pytransdecoder==0.4.0` does **not** drop in: the `git rev-parse` fails. Two options:
1. Simplest: `git clone https://github.com/TransDecoder/PyTransDecoder && git -C PyTransDecoder checkout v0.4.0`, and
   point the env vars at the clone: `PYTD_V7=` (run_tool.sh, gm12878_final_v7.sh), `PYTD_REL=` (run_gm_v1.sh). Then
   `pytd_version.txt` holds the v0.4.0 hash.
2. pip install into `env/td` and edit the scripts: drop `PYTHONPATH=...` and replace the rev-parse line with
   `python -c 'import importlib.metadata as m; print(m.version("pytransdecoder"))' > pytd_version.txt`.
   `homology_series_v7h*.sh` and `08_homology_series_v7h.sh` also **assert `pytd_version.txt == a5b24c3`**, and
   `homology_series_v7h*.sh` hardcode `PYTHONPATH=$P/work/trees/PY_hstart`. Both must be edited.

---------------------------------------------------------------------------------------------------------------
## A. GM12878 upstream data (shared by E01-E06)

### E11 Trinity genome-guided assembly of GM12878 and upstream data  (status: in paper, Methods)
- **Feeds**: Methods "GM12878 transcriptome assemblies"; `06_gm12878_inputs.tsv` (notebook 06); input to E03/E05/E06/E01/E02.
- **Inputs**:
  - ENCODE ENCSR000AED rep 1: `ENCFF001REK` (R1, md5 0055876c86d1909b7b7ff0e5497bd92b) and `ENCFF001REJ` (R2, md5 b484a76ccf390e100d60f373a9c54dd6).
    URL form `https://www.encodeproject.org/files/<ENCFF>/@@download/<ENCFF>.fastq.gz` -> `$P/data/gm12878/reads/short/`.
    **These were deleted on 2026-10-29** (notes/DISK_CLEANUP_2026-10-29.md).
  - ENCSR838WFC rep 2 PacBio FLNC `ENCFF694DIE.fastq.gz` (md5 38c1eaea4283c60b6d8f0b3025d540cd) -> `$P/data/gm12878/reads/long/` (present).
  - GENCODE v50 `GRCh38.primary_assembly.genome.fa.gz` (da1a1125...) and `gencode.v50.primary_assembly.annotation.gtf.gz`
    (289b91e5...), from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/ -> `$P/data/gm12878/ref/` (present).
- **Software**: STAR 2.7.11b, StringTie 3.0.3, minimap2 2.31-r1302, samtools 1.24, gffread 0.12.9 (`env/rnaseq_asm`);
  Trinity 2.15.2 + salmon 1.10.3 pinned (`env/trinity215`; `micromamba create -p env/trinity215 -c conda-forge -c bioconda --override-channels trinity=2.15.2 salmon=1.10.3`);
  LRAA v0.44.3 (`us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:0.44.3` via apptainer, `$P/data/gm12878/longread_lraa/image/lraa-core_0.44.3.sif`).
  **Must** replace the serial bioconda ParaFly:
  `B=$P/env/trinity215/bin/trinity-plugins/BIN; mv $B/ParaFly $B/ParaFly.bioconda_serial; cp genome_annotation/scripts/data_prep/ParaFly_xargs.sh $B/ParaFly; chmod +x $B/ParaFly`.
  Unset `PERL5LIB` (trinity_gg.sh does this).
- **Steps** (`genome_annotation/scripts/data_prep/`, see its README; `D=$P/data/gm12878`):
  1. `download.sh <outdir> <url> <md5>` for each FASTQ.
  2. `star_index.sh` (`--sjdbGTFfile gencode.v50... --sjdbOverhang 100`), `star_align_stringtie.sh` (2-pass, `--outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS`;
     then `stringtie ... --rf -p 4 -l STRG -o stringtie.gtf`, `gffread -w stringtie.transcripts.fa`).
  3. `longread_minimap2.sh`, `pull_lraa.sh`, `lraa_denovo.sh`, `lraa_post.sh` (-> `$D/longread_lraa/lraa_denovo/LRAA.{gtf,transcripts.fa}`).
  4. `CPU=16 bash trinity_gg.sh`: `Trinity --genome_guided_bam $D/star_stringtie/GM12878_ENCSR000AED_rep1.STAR.bam --genome_guided_max_intron 100000 --SS_lib_type RF --max_memory 50G --CPU 16`,
     then `minimap2 -t 16 -ax splice:hq -uf --secondary=no` -> `$D/trinity/Trinity.minimap2.bam`.
- **Outputs**: no committed table of its own; counts appear in `data/gm12878/final_v7/06_gm12878_inputs.tsv`.
  Kept in `$P/data/gm12878/`: `trinity/Trinity.fasta` (428,960 tx), `Trinity.minimap2.bam`, `star_stringtie/stringtie.{gtf,transcripts.fa}`, `longread_lraa/lraa_denovo/LRAA.*`.
- **Cost**: Trinity GG ~11 h of partition assembly at CPU 16 (plus ~4 h partitioning), `--max_memory 50G`; STAR index ~31 GB on disk; minimap2 164 s.
  Rebuilding the assemblies is **not verified byte-identical**; the kept `Trinity.fasta` / StringTie / LRAA files are the reference copies.

---------------------------------------------------------------------------------------------------------------
## B. GM12878 PyTD runs and evaluation (E03, E05) and their evidence (E01, E02)

Paper numbers are from the **release** runs: `final_v7_*` (E25, 3a793ff), byte-identical to `release_pin/gm/<asm>/rel`
(d3892b5, E40). The development records E01/E02/E03/E05 (48c870c, `final_abinitio`) describe the methods; their
development tables (`data/gm12878/gm12878_*.tsv`, `eval_summaries/*_final_abinitio*`) are only used as comparison
columns in notebooks 06/07.

### E03 GM12878 three-assembly annotation vs MANE/GENCODE  (status: in paper; release values)
- **Feeds**: Table `tab:gm12878`, Fig. 6b, Results "annotating three GM12878 assemblies", Methods; `gm_runs.tsv` also feeds Fig. 7 / E13-E40 grid text.
- **Inputs**: section A assemblies; MANE Select v1.5 `$P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz`; GENCODE v50 GTF (`$D/ref/`);
  GENCODE Ribo-seq ORF catalogues phase 1 (`$P/data/riboseq_orfs/phase1/Ribo-seq_ORFs.all.bed`, GENCODE 35) and phase 2 (`phase2/Ribo-seq_ORFs.comprehensive.bed`, GENCODE 45).
  Homology mode: `$P/data/homology/{sprot_minus_9606.dmnd,Pfam-A.hmm}` (Swiss-Prot 2026_03, Pfam-A 38.2).
- **Software**: PyTD (section 0); DIAMOND 2.2.8, HMMER 3.4 (`env/td`). **Rerun with 0.4.0 reproduces the paper: yes.**
- **Steps** (release path):
  1. `PYTD_V7=<v0.4.0 clone>` and run what `$P/work/gm12878_td/launch_v7.sh` ran (**not in the repo**; content):
     `genome_annotation/scripts/gm12878_final_v7.sh <asm> <fasta> <structures> abinitio 4 [precise|sensitive]` for lraa, stringtie, trinity (default),
     lraa+stringtie with precise and sensitive, then `gm12878_final_v7.sh trinity $D/trinity/Trinity.fasta $D/trinity/Trinity.minimap2.bam homology 8`.
     (fasta/structures: `$D/longread_lraa/lraa_denovo/LRAA.{transcripts.fa,gtf}`, `$D/star_stringtie/stringtie.{transcripts.fa,gtf}`, `$D/trinity/Trinity.{fasta,minimap2.bam}`.)
     Each run: `longorfs -S`; `predict` [+ `--retain-blastp-hits/--retain-pfam-hits`]; `genome-map`; `eval_genome_cds.py` -> `$P/work/gm12878_td/eval/<asm>_final_v7_<mode>{.summary,.ref_status,.pred_status}.tsv`;
     `analyses/riboseq/match_riboseq.py` -> `$P/work/gm12878_td/ribo/`.
  2. Equivalent pinned rerun (E40): `PYTD_REL=<clone> bash pipeline/release_pin/run_gm_v1.sh <asm> <cfg>` (`cfg` = rel, precise, sensitive, m50_noq, m100_noq, norefine, a<A>_w<W>, retonly)
     -> `$P/work/release_pin/gm/<asm>/<cfg>/`; then `python pipeline/release_pin/12_gm_tables_v1.py` -> `data/gm12878/release_pin/{gm_runs,no_overlap_by_length}.tsv`.
     Configs other than rel/m*_noq symlink rel's LongOrfs (run rel first). **`retonly` reads hits from `$P/work/gm12878_td/<asm>/final_v7_homology/{blastp.outfmt6,pfam.domtblout}`.**
  3. `bash pipeline/derive/06_gm12878_application.sh` -> `data/gm12878/final_v7/06_*.tsv` (needs `$P/work/gm12878_td/eval/*final_v7_abinitio*`,
     and for the "earlier" rows also the 48c870c `final_abinitio` runs and `eval/status/`; it sets `PYTHONPATH=$P/work/trees/PY_v7`).
  4. Eval summaries were copied to `data/benchmark/final_v7/eval.*_final_v7_abinitio.summary.tsv`.
- **Outputs**: `data/benchmark/final_v7/eval.*summary.tsv`, `data/gm12878/final_v7/06_*`, `data/gm12878/release_pin/gm_runs.tsv`; consumed by `notebooks/06_gm12878_application.Rmd` (Fig. 6, Table) and 08.
  Development: `data/gm12878/eval_summaries/*_final_{abinitio,homology}.summary.tsv`, `figures/fig6_gm12878/data/pred_class_by_length.tsv` (old matplotlib fig, superseded by notebook 06).
- **Cost**: release_pin rel runtimes 480 s (LRAA), 783 s (StringTie), 1,898 s (Trinity), single process. Trinity homology adds DIAMOND (633,606 peptides) and hmmsearch `--cpu 8` (hours).

### E05 Expressed-gene sensitivity  (status: in paper; release values)
- **Feeds**: Fig. 6a, Results "restricted to expressed MANE genes" (74-81%, 11,076 genes), Methods.
- **Inputs**: STAR BAM (section A; **deleted**, rebuild) -> StringTie `-e` quantification `$P/data/gm12878/expr/gencode_v50.ref_quant.gtf` (**present**, 950 MB); E03 `ref_status.tsv`.
- **Software**: StringTie 3.0.3, Python (`env/td`). 0.4.0: yes.
- **Steps**:
  1. `genome_annotation/scripts/data_prep/stringtie_expr.sh` (reconstruction of the one-off command:
     `stringtie <STAR.bam> --rf -e -G ref/gencode.v50.primary_assembly.annotation.gtf -p 8 -A expr/gencode_v50.gene_abund.tsv -o expr/gencode_v50.ref_quant.gtf`). Skip if `expr/` is kept.
  2. Release: inside `pipeline/derive/06_gm12878_application.sh`: `expressed_sensitivity.py $P/data/gm12878/expr/gencode_v50.ref_quant.gtf lraa=$E/lraa_final_v7_abinitio.ref_status.tsv stringtie=... trinity=...`.
- **Outputs**: `data/gm12878/final_v7/06_gm12878_v7_expressed_sensitivity.tsv` (notebook 06). Development: `data/gm12878/gm12878_expressed_sensitivity.tsv`.

### E01 GM12878 Ribo-seq support  (status: in paper; release values)
- **Feeds**: Fig. 6c,e; Results "Translation and protein evidence"; Methods "Ribosome profiling support"; rule-support numbers (E40).
- **Inputs**: Cenik et al. 2015 GSE65912, SRR1803147 + SRR1803148 from `https://sra-pub-run-odp.s3.amazonaws.com/sra/<run>/<run>`;
  rRNA/tRNA set = U13369.1 + GENCODE v50 rRNA/rRNA_pseudogene/Mt_rRNA/Mt_tRNA (575 seqs); STAR index of section A.
  Kept: `$P/data/gm12878/riboseq_cenik2015/footprints.genome.bam(.bai)`, `rrna/` bowtie2 index, `logs/`, `star/*.Log.final.out`. SRA/FASTQ/per-run BAMs deleted.
- **Software**: sra-tools 3.4.1 (`env/sratools`), cutadapt 5.2, bowtie2 2.5.5, samtools 1.24, pysam 0.24.1 (`env/riboseq`), STAR 2.7.11b. 0.4.0: yes.
- **Steps** (only if the footprint BAM is lost): `genome_annotation/scripts/gm12878_riboseq_download.sh`, `gm12878_riboseq_rrna_index.sh`, `gm12878_riboseq_process.sh`
  (`cutadapt -a AGATCGGAAGAGCACACGTCT -O 3 --discard-untrimmed --trim-n -m 25 -M 34`; `bowtie2 --very-sensitive -L 15 --un-gz`; `STAR --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMultimapNmax 1`).
  Scoring (release): **no committed script** for the `final_v7_scoring` run (gap 1). Reconstruction from E01 + E25:
  ```bash
  cd $P; G=work/gm12878_td; E=$G/eval
  env/riboseq/bin/python $T/analyses/riboseq/riboseq_orf_support.py --bam data/gm12878/riboseq_cenik2015/footprints.genome.bam \
    --mane data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
    --pred lraa=$G/lraa/final_v7_abinitio/transcripts.fa.transdecoder.genome.gff3,$E/lraa_final_v7_abinitio.pred_status.tsv \
    [--pred lraa_precise=..., lraa_sensitive=..., stringtie..., stringtie_precise..., stringtie_sensitive..., trinity=...] \
    -o $G/final_v7_scoring/ribo/gm12878_v7
  cp $G/final_v7_scoring/ribo/gm12878_v7.summary.tsv $T/data/benchmark/final_v7/gm12878_v7_riboseq_support.tsv   # verified identical today
  bash $T/pipeline/derive/07_gm12878_evidence_v7.sh   # read processing, by-expression (expr_matched_support.awk), MS tables
  ```
  Pinned alternative: `bash pipeline/release_pin/13_gm_support_v1.sh` (rel + m50_noq from `work/release_pin/gm`) -> `data/gm12878/release_pin/riboseq_summary.tsv` (its rel rows reproduce E25 exactly).
- **Outputs**: `data/benchmark/final_v7/gm12878_v7_riboseq_{support,calibration}.tsv`, `data/gm12878/final_v7/gm12878_v7_riboseq_{read_processing,stringtie_by_expression}.tsv`,
  `data/gm12878/release_pin/{riboseq_summary,rule_support}.tsv` -> notebooks 06 (Fig. 6c,e) and 07. Development: `data/gm12878/gm12878_riboseq_*.tsv`.
- **Cost**: not recorded for scoring; download ~25 GB.

### E02 GM12878 MS peptide support  (status: in paper; release values. C3orf86P claim: only in `--sensitive`, not in text)
- **Feeds**: Fig. 6d; Results/Methods "Mass-spectrometry support"; Discussion.
- **Inputs**: Sinitcyn et al. 2023, MassIVE MSV000086944 (PXD024364; PRIDE has no record) `search/GM12878.zip` (5.52 GB, md5 f062982347bb27a12b42c52823275b08; `peptides.txt` md5 8bf9a444...);
  refdb GENCODE v50 `pc_translations` + UniProt UP000005640 with isoforms. All under `$P/data/gm12878/proteome_sinitcyn2023/` (present, 6.4 GB).
- **Software**: `env/proteo` (Python 3.11, pyahocorasick, pandas). 0.4.0: yes.
- **Steps** (release; same gap as E01): `env/proteo/bin/python -I analyses/proteomics/ms_peptide_support.py --peptides $D/search/GM12878/peptides.txt --ref <gencode> --ref <uniprot> --pred <label>=<final_v7 pep>,<pred_status> ... -o $P/work/gm12878_td/final_v7_scoring/proteomics`;
  copy `ms_peptide_support.summary.tsv` -> `data/benchmark/final_v7/gm12878_v7_ms_peptide_support.tsv`; then `pipeline/derive/07_gm12878_evidence_v7.sh` (runs `novel_peptide_loci.py`). Pinned: `13_gm_support_v1.sh` -> `release_pin/ms_summary.tsv`.
- **Outputs**: `data/benchmark/final_v7/gm12878_v7_ms_peptide_support.tsv`, `data/gm12878/final_v7/gm12878_v7_ms_*`, `release_pin/ms_summary.tsv` -> notebooks 06/07.
- **Cost**: 175 s, 2.2 GB RAM.

### E06 Homology restriction series (clade-excluded Swiss-Prot)  (status: in paper, release series `v7h`; development series not adopted)
- **Feeds**: Fig. 7c,d; Supplementary Table S7 (`tab:hseries`); Results "Homology evidence and evolutionary distance", "Homology to transposable elements".
- **Inputs**: Swiss-Prot 2026_03 `$P/data/homology/uniprot_sprot.fasta.gz`; NCBI `$P/data/taxonomy/nodes.dmp` (2026-10-21); Pfam-A 38.2 `$P/data/homology/Pfam-A.hmm`;
  clade DBs `$P/data/homology/sprot_minus_clade_{9443,40674,7742,33208}.{fasta,dmnd}` + `sprot_minus_9606.dmnd` (all present);
  UCSC hg38 `rmsk.txt.gz` `$P/data/repeats/ucsc_hg38_rmsk/` (TE screen). Upstream: E03 `final_v7_abinitio` LongOrfs (prep `cmp`s against them).
- **Software**: PyTD a5b24c3 (= d3892b5 tree; 0.4.0: yes), DIAMOND (`--max-target-seqs 1 --evalue 1e-5 --outfmt 6`, default sensitivity), HMMER `hmmsearch --cut_ga --noali`, bedtools 2.31.1 (`env/bedtools`).
- **Steps**:
  1. DBs (already in data): `python3 -I genome_annotation/scripts/make_sprot_clade_excluded.py uniprot_sprot.fasta.gz nodes.dmp <taxid> > sprot_minus_clade_<taxid>.fasta; diamond makedb ...`.
  2. `homology_series_v7h_prep.sh <asm> [threads]` (fresh LongOrfs, `cmp` with `final_v7_abinitio`, DIAMOND x5 clades; **copies Pfam** from `$P/work/gm12878_td/<asm>/pfam.domtblout`
     (LRAA/StringTie = 12-chunk searches, `homology_series_driver.sh` step 2: `hmmsearch --cpu 1 --cut_ga --noali` per chunk, `xargs -P 12`) or `trinity/final_v7_homology/pfam.domtblout`).
     Without these files, run hmmsearch on `_longorfs/.../longest_orfs.pep` (E-values differ, predictions do not).
  3. `homology_series_v7h_driver_v2.sh` (36 configs via `homology_series_v7h.sh <asm> <taxid|none> <pfam|nopfam>`, 3 at a time) -> `$P/work/gm12878_td/<asm>/hseries_v7h/<cfg>`, `eval/hseries_v7h/`.
     (driver_v1 ran prep + 9 Trinity configs; v2 is the complete list minus those. For a fresh run, run all 36.)
  4. `bash pipeline/derive/08_homology_series_v7h.sh` (asserts a5b24c3), `python pipeline/derive/08_homology_series_v7h_table.py`, `bash pipeline/derive/08_te_screen_check_v7h.sh`.
- **Outputs**: `data/gm12878/gm12878_homology_series_v7h.{tsv,tex}`, `eval_summaries/hseries_v7h/`, `gm12878_homology_primate_retention_v7h.txt`, `gm12878_te_screen_v7h.tsv` -> notebook 08 (Fig. 7c,d), supplement.tex S7.
  Development (not adopted): `gm12878_homology_series.tsv`, `eval_summaries/hseries/` (48c870c, `homology_series_driver.sh`).
- **Cost**: 36 Predict+genome-map+eval runs (3 parallel; Trinity ~30 min each); hmmsearch is the expensive part (Trinity 633,606 peptides).

---------------------------------------------------------------------------------------------------------------
## C. Simulated Trinity pipeline (E30-E32), `pipeline/simtrinity/` and `pipeline/simtrinity_athal/`

Shared software: salmon 1.10.3 + Trinity 2.15.2 (ParaFly drop-in) `env/trinity215`; Polyester 1.38.0, R 4.3.3, Biostrings 2.70.1 `env/polyester`;
STAR 2.7.11b, minimap2 2.31, samtools 1.24, gffread `env/rnaseq_asm`; predictors via unmodified `benchmark_v2/scripts/run_tool.sh`:
`pytd_v7` (PyTD 3a793ff; **0.4.0 reproduces: yes**, ab initio path byte-identical, but the input assembly must be the kept one), `perl_td_v600`
(TransDecoder v6.0.0 66d4712, `TD_V600=$P/work/trees/TD_base`), Prodigal 2.6.3, GeneMarkS-T 5.1 (licensed, `$HOME/software/genemarks_t/dist`, key `~/.gm_key`),
CPAT 3.0.5 (human model, `$P/tools/models/cpat`), CPC2 1.0.1, RNAsamba 0.2.5 (`$P/tools/models/rnasamba`), Borf 1.2. Longest-ORF baseline = `pytd_v7/longest_orf.preds.tsv`.

### E30 Simulated GM12878 Trinity assembly benchmark  (status: in paper, Results text; d1 is the full-depth level in Supp. Figs S5/S6)
- **Feeds**: results_benchmark.tex simulated-assembly paragraph; results_coding_failures.tex (start rules on simulated assemblies); notebook 02 section "Simulated GM12878 Trinity assembly" (`figures/figS_simtrinity.*`, not included in supplement.tex).
- **Inputs**: GENCODE v50 (section A); `gffread -w` -> `$P/data/simtrinity/ref/gencode.v50.transcripts.fa` (present); GM12878 short reads (salmon only; **deleted**);
  STAR index of section A (**deleted**); MANE v1.5. Kept: `$P/data/simtrinity/{ref,salmon/quant,sim/chunk_*.{fa,counts.tsv},labels,inputs,trinity/Trinity.fasta,Trinity.minimap2.bam}`.
  Simulated FASTQ, STAR BAM, salmon index, Trinity working dirs **deleted**.
- **Steps**:
  1. `python pipeline/simtrinity/01_gencode_tables.py <gencode.gtf> $P/data/simtrinity/ref/gencode.v50.transcripts.fa $P/data/simtrinity/ref`
  2. `bash 02_salmon.sh` (`-k 31 --keepDuplicates`; `quant -l ISR --seqBias --gcBias -p 6 --validateMappings`).
  3. `python 03_make_counts.py $P/data/simtrinity/salmon/quant ref/kept.tsv ref/kept.fa $P/data/simtrinity/sim 5`
  4. `bash 06_simulate.sh` (Polyester, seed 20261025+k, 2 internal patches; R1 = Polyester mate 2).
  5. `bash 07_star.sh` (index `$P/data/gm12878/star_stringtie/star_index`, `--limitBAMsortRAM 40000000000`, 6 threads; peak 37 GB -> was run on methods02 via `$P/work/simtrinity/run_star_methods02.sh`, not in repo).
  6. `bash 08_trinity_gg.sh` (`--min_kmer_cov 2 --max_memory 50G --CPU 8`; then minimap2 `-uf`).
  7. Label + E29 statistic + runs + score for d1: **no single committed script** (gap 2). Commands = `12_depth_label_and_run.sh` with d1 paths:
     `label_assembly.py --bam $P/data/simtrinity/trinity/Trinity.minimap2.bam --fasta .../Trinity.fasta --gencode ref/gencode_v50.tx.tsv.gz --kept ref/kept.tsv --salmon salmon/quant/quant.sf --out $P/data/simtrinity/labels/sim_trinity`;
     `pipeline/derive/03_utr_gm12878_truth.py ... labels/sim_trinity_e29`; inputs `$P/data/simtrinity/inputs/fragments_simtrinity.{fa,truth.tsv}`;
     runs as `$P/work/simtrinity/launch_simtrinity_runs.sh` (not in repo: `run_tool.sh <tool> fragments_simtrinity.fa $P/work/simtrinity/runs/<tool> human` for pytd_v7 perl_td_v600 gmst prodigal cpat cpc2 rnasamba borf, `THREADS=1`, 4 parallel);
     `score_simtrinity.py <fa> $P/data/simtrinity/labels/sim_trinity $P/work/simtrinity/runs $T/data/simtrinity` (filtered preds land in `$P/work/simtrinity/filtered_preds`).
  8. Real-assembly comparison: `label_assembly.py` on `$P/data/gm12878/trinity/` -> `$P/work/simtrinity/realtrinity_label/real_trinity.*` (work only); then
     `09_summarize.py labels/sim_trinity work/simtrinity/realtrinity_label/real_trinity sim/counts.tsv <sim_e29.pairs.tsv.gz> <real_e29.pairs.tsv.gz> $T/data/simtrinity`, `10_incompatible_breakdown.py ... > incompatible_breakdown.tsv`.
  9. `scripts/knit_notebooks.sh 02_benchmark.Rmd`.
- **Outputs**: `data/simtrinity/*.tsv` (scores_overall/by_type/by_expr, negatives, excluded_preds, pred_types, label_counts*, type_by_expr, complete_by_source_expr, e29_five_prime, strand_check_*, selection_counts, sim_params) -> notebook 02.
- **Cost**: Polyester 1.2-1.4 h per chunk (5 parallel); STAR 37 GB peak; Trinity GG 8.5 h at CPU 8; pytd_v7 913 s.

### E31 Read-depth ablation (human)  (status: supplement, Figs S5, S6)
- **Feeds**: Supplementary Figs S5 (`figS5_simtrinity_depth`) and S6 (`figS6_simtrinity_depth_sens`); Results text.
- **Inputs**: E30 STAR BAM `$P/data/simtrinity/star/sim_GM12878.STAR.bam` (**deleted**; rebuild with E30 steps 1-5). Kept per level: `$P/data/simtrinity/depth/d{2..32}/{Trinity.fasta,Trinity.minimap2.bam,labels,inputs,pair_counts.tsv}`.
- **Steps** (driver `$P/work/simtrinity/depth/run_lane_v1.sh`, not in repo; at most 2 Trinity runs at a time):
  `bash pipeline/simtrinity/11_depth_subsample.sh dN <frac>` (seed 20261026; nested); `CPU=4 bash 08b_trinity_gg_level.sh dN`; `bash 12_depth_label_and_run.sh dN`;
  then `$P/env/td/bin/python 13_depth_collect.py $T/data/simtrinity/depth d1 d2 d4 d8 d16 d32` (or `13b_depth_collect.py`, byte-identical); knit notebook 02.
- **Outputs**: `data/simtrinity/depth/*.tsv` -> notebook 02 -> `figures/figS_simtrinity_depth{,_sens}` (copied to manuscript as figS5/S6).
- **Cost**: per level one Trinity GG run (CPU 4) + 8 tools; d32 had an NFS flock stall (E31 caveat).

### E32 Simulated Arabidopsis Trinity assembly + depth series  (status: supplement Fig. S7 + Results text)
- **Inputs**: TAIR12 `$P/data/annot/athal/{Athaliana_TAIR12_genomic_Chr_softmasked.fa.gz,TAIR12_genome_annotation_TAIR.gff3.gz}` (present);
  ENA ERR4140056 (PRJEB38272, Col-0 leaf rep 3; md5 R1 bbac8950..., R2 394b6671...) from the ENA FTP -> `$P/data/simtrinity_athal/reads/` (**deleted**; no committed download script, gap 3).
  Kept: `$P/data/simtrinity_athal/{ref,salmon/quant,sim/chunk_*,depth/dN/{Trinity.fasta,Trinity.minimap2.bam,labels,inputs}}`; STAR BAM/index, salmon index, sim FASTQ deleted.
- **Steps** (E32 block, `A=pipeline/simtrinity_athal`, `S=pipeline/simtrinity`):
  `01_tair12_tables.py <gff3.gz> <fa.gz> $P/data/simtrinity_athal/ref`; `02_salmon_athal.sh`; `03_star_index_athal.sh`; `$S/03_make_counts.py ... sim 4` (note: E32's `$R` = `$P/data/simtrinity_athal/ref`);
  `06_simulate_athal.sh` (seed 20261027+k); `07_star_athal.sh` (~5 GB, local); `CPU=8 bash $A/run_level_athal.sh d1:1`; `run_level_athal.sh d2:0.5 ... d32:0.03125` (CPU 4, 2 lanes);
  `$S/13b_depth_collect.py --data-root $P/data/simtrinity_athal --work-root $P/work/simtrinity_athal $T/data/simtrinity_athal/depth d1 d2 d4 d8 d16 d32`;
  `$A/14_negatives_by_class.py $T/data/simtrinity_athal/depth/negatives_by_class.tsv d1 ... d32`. Tools with organism `athal` (TD `-G Universal`).
- **Outputs**: `data/simtrinity_athal/{selection_counts,sim_params,strand_check_genome}.tsv`, `data/simtrinity_athal/depth/*.tsv` -> notebook 02 -> `figures/figS_simtrinity_athal` (manuscript figS7).
- **Cost**: not recorded per step; STAR ~5 GB.

---------------------------------------------------------------------------------------------------------------
## PART 2. Things in $P/work that are inputs or expensive to regenerate

Sizes from `du -sh` on 2026-11-03. "Data copy" = an equivalent exists in `$P/data` or the repo.

| path | size | what | needed by | data copy? | recommendation |
|---|---|---|---|---|---|
| `work/trees/` (PY_v7, PY_hstart, PY_rel_d3892b5, PY_48c870c, TD_base, TD_rel_41fc5ee, TD_v571...) | 1.4G | git worktrees of `$P/repos/*` | every script default (`PYTD_V7`, `TD_V600`, `PY_hstart`) | commits in `$P/repos` (if kept) and GitHub | safe to delete **if $P/repos is kept** (then `git worktree prune`); otherwise recreate from GitHub at the commits in section 0 |
| `work/gm12878_td/trinity/final_v7_homology/{blastp.outfmt6,pfam.domtblout}` (12 MB + 55 MB; dir 392M) | 392M | DIAMOND + Pfam hits on Trinity LongOrfs (633,606 peptides) | E03 Trinity homology, E40 `retonly` (hard path), E06 v7h prep (copies Pfam) | none | **move the two hit files to `$P/data/gm12878/homology_hits/trinity/`** (hmmsearch = hours) |
| `work/gm12878_td/{lraa,stringtie}/pfam.domtblout` (24 MB, 25 MB) + `pfam_chunks/` (34M, 51M) | ~135M | Pfam hits on LRAA/StringTie LongOrfs | E06 v7h prep (hard path) | none | **move `pfam.domtblout` to `$P/data/gm12878/homology_hits/<asm>/`**; chunks deletable |
| `work/gm12878_td/<asm>/hseries_v7h/blastp.minus_*.outfmt6` (15 files, 1.4-12 MB) + `pfam.domtblout` | ~120M of 1.3G/2.0G/4.5G | clade-excluded DIAMOND hits (verified identical to dev `hseries/`) | E06 Fig. 7c,d / Table S7 | none | **move hit files** to the same data dir; run dirs regenerable (~36 Predict runs) |
| `work/gm12878_td/eval/` (incl. `status/`, `hseries_v7h/`) | 572M | `ref_status`/`pred_status`/summaries for all GM runs | E01/E02 scoring, E05, 06/07/08 derive scripts, fig6 old panel b | summaries only (repo) | regenerable by rerunning E03/E06 (~1 h total); keep if cheap; otherwise delete |
| `work/gm12878_td/{lraa,stringtie,trinity}/final_v7_abinitio*` | 103M/153M/348M | release ab initio predictions (= release_pin `rel`) | E01/E02/E05/E06 prep `cmp` | release_pin/gm duplicates | regenerable (8-32 min each); delete one of the two copies at most |
| `work/gm12878_td/{*/final_abinitio,*/hseries,eval/status,...}` (48c870c) | ~5G of 21G | development runs | only the "earlier" comparison rows of `06_gm12878_tables.py` and `figures/fig6_gm12878/derive_class_by_length.py` | committed tables | safe to delete (needs PyTD 48c870c to regenerate; tables are committed) |
| `work/gm12878_td/final_v7_scoring/` | 104M | per-ORF Ribo-seq/MS tables of the release runs | `07_gm12878_evidence_v7.sh`, `08_rule_support_v7.py` | summaries in repo | regenerable in minutes from the footprint BAM / peptides.txt; delete OK |
| `work/gm12878_td/{ribo,ribo_gm12878,proteomics,proteomics_repro,te_screen_v7h,comparators}` | 16M/4M/45M/45M/92M/688M | dev scoring outputs, TE screen intermediates, comparator runs | E01/E02 dev, E06 TE screen | summaries in repo | delete OK (regenerable) |
| `work/gm12878_td/launch_v7.sh`, `overnight.sh`, `final_defaults.sh`, `homology_*.sh`, `logs/` | <3M | launchers; `launch_v7.sh` is **not in the repo** | E03 release runs | others are in `genome_annotation/scripts` | **commit `launch_v7.sh`** (or paste into REPRODUCE.md); logs optional |
| `work/release_pin/gm/` | 6.1G | pinned GM runs (55 configs) | `12_gm_tables_v1.py`, `13_gm_support_v1.sh` | tables in `data/gm12878/release_pin/` | regenerable (~1 day at 4 processes with grid); delete OK |
| `work/release_pin/gm_scoring/` | 108M | protein-match/Ribo/MS scoring of pinned runs | `release_pin/*.tsv` | yes (repo tables) | delete OK |
| `work/simtrinity/realtrinity_label/` | 20M | labels of the real GM12878 Trinity GG vs simulated transcript set | E30 `09_summarize.py` (real column, e29 statistic) | none | **move to `$P/data/simtrinity/labels/real_trinity/`** (cheap but small) |
| `work/simtrinity/runs/` | 612M | d1 predictor outputs (`preds.tsv`) incl. GeneMarkS-T (licensed) and TD v6 | E30 scores | scores in repo | keep `*/preds.tsv` + `runtime_sec.txt` (small) in `$P/data/simtrinity/runs_d1/`, rest delete |
| `work/simtrinity/filtered_preds/` | 16M | d1 scored predictions | `13_depth_collect.py` d1 (hard path) | none | move with runs (collector reads it) |
| `work/simtrinity/depth/dN/{runs,scores}` | 1.6G | per-level predictor runs and scores | `13_depth_collect.py` (reads `work/simtrinity/depth/dN/scores`) | collected tables in repo | keep `scores/` (small) or delete; regenerable from kept assemblies (~hours per level) |
| `work/simtrinity_athal/depth/dN/{runs,scores}` | 1.3G | same, Arabidopsis | `13b_depth_collect.py --work-root` | tables in repo | same as above |
| `work/simtrinity/{run_star_methods02.sh,launch_simtrinity_runs.sh,depth/run_lane_v1.sh}` | <1M | launchers **not in repo** | E30/E31 | none | **commit into `pipeline/simtrinity/`** |
| `work/simtrinity/{strandcheck,polytest,diag,scores_e39,strand_sample.tsv,tx_strand.tsv}` | ~45M | strand checks, Polyester test, E39 TD2 scores | E30 checks, E39 | `strand_check_*.tsv` in repo | delete OK |

Already outside `$P/work` and kept (no action): Trinity/StringTie/LRAA assemblies, Ribo-seq footprint BAM, MS search, StringTie `-e` quantification, homology DBs, simulated `chunk_*.fa`/counts, all simulated assemblies and labels.
Already deleted 2026-10-29 (regenerable but costly): GM12878 short-read FASTQ, STAR index + BAM, Ribo-seq FASTQ/BAMs, simulated FASTQs and STAR BAMs, salmon indexes, ERR4140056 reads.

---------------------------------------------------------------------------------------------------------------
## PART 3. Gaps

1. **GM12878 release Ribo-seq/MS scoring (`work/gm12878_td/final_v7_scoring/`) has no committed script.** E25 gives only "E01/E02 arguments on the v7 files". The pinned
   `pipeline/release_pin/13_gm_support_v1.sh` reproduces the default rows exactly but covers rel + m50_noq only (no precise/sensitive labels, which
   `data/benchmark/final_v7/gm12878_v7_{riboseq,ms_peptide}_support.tsv` and `07_gm12878_evidence_v7.sh` use). Suggest committing a `final_v7_scoring` script.
2. **Launchers only in `$P/work`:** `gm12878_td/launch_v7.sh` (E03/E25 release runs), `simtrinity/launch_simtrinity_runs.sh` (E30 d1 runs), `simtrinity/run_star_methods02.sh`,
   `simtrinity/depth/run_lane_v1.sh` (E31). E30 d1 labelling/E29/scoring and `09_summarize.py`/`10_incompatible_breakdown.py` invocations are not scripted (usage lines only).
   The real-Trinity labelling command (`realtrinity_label`) is not recorded.
3. **ERR4140056 download** (E32) is not scripted; only "ENA FTP" + md5 prefixes in E32. Full md5s should be copied from ENA into REPRODUCE.md.
4. **Raw inputs needed to restart from reads were deleted** (DISK_CLEANUP_2026-10-29): GM12878 short reads, STAR index/BAM, simulated FASTQs/BAMs, ERR4140056.
   Reassembly (Trinity GG real and simulated, 6+6 depth levels) is **not verified byte-identical**, so exact paper numbers require the kept assemblies in `$P/data`, not a rebuild.
5. **E05 expression** was a one-off command; `stringtie_expr.sh` is a reconstruction (output kept in `$P/data/gm12878/expr/`).
6. **Hardcoded version checks** in `homology_series_v7h{,_prep}.sh` and `pipeline/derive/08_homology_series_v7h.sh` (`a5b24c3`) and `PYTHONPATH=$P/work/trees/PY_hstart`; `06_gm12878_application.sh` uses `PY_v7`;
   all `git rev-parse` calls require a git checkout (section 0).
7. **Development tables** (`data/gm12878/gm12878_*` without `v7`, `eval_summaries/*final_abinitio*`, `hseries/`, `figures/fig6_gm12878/data/`) need PyTD 48c870c
   (`work/trees/PY_48c870c`) to regenerate; 0.4.0 will not reproduce them. They are comparison-only.
8. **Per-ORF outputs never committed** (only in work): E01 `orf_support.tsv`, E02 `model_hits`/`novel_specific_hits`, E03 `ref_status`/`pred_status`, predictions GFF3/pep.
   All regenerable from section B commands; nothing in the paper reads them directly.
9. **Cost not recorded** for E32 steps, E06 hmmsearch wall time, E01 scoring, and E31 per-level Trinity runtime.
10. `$P/repos` and `$P/tools` (CodAn, CPAT/RNAsamba models) are neither data nor env: confirm they survive, or the worktrees and comparators cannot be recreated locally. GeneMarkS-T lives in `$HOME` (licensed).
