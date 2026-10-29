# Disk cleanup 2026-10-29 (/broad/hptmp quota)

`P=/broad/hptmp/bhaas/TRANSDECODER_PAPER`, `T=$P/TransDecoderPaper`. Scope: `$P` only (not `$P/env`, not git-tracked
files in `T`). Authorised by Brian (quota over its soft limit, grace period about 3 days). Deletions ran niced on methods01,
with at most 2 concurrent rm lanes. No process on methods01 or methods02 was using `$P` (checked with `ps` and `/proc/<pid>/cwd`).

## Quota (`quota -s`, helium-nfs:/ifs/data/broad/hptmp)

| | space | soft | hard | grace | files |
|---|---|---|---|---|---|
| before | 3142G* | 3072G | 10240G | 4days | 2591k |
| after  | 1973G  | 3072G | 10240G | -     | 1167k |

The quota dropped by 1,169 GB, but `du` (disk blocks) counts about 587 GB freed. The quota accounting is larger than `du`
on this Isilon volume (logical size and/or protection overhead). The per-row sizes below are `du -sk`, before minus after.

## What the paper reads (checked before and after)

* The notebooks (`notebooks/*.Rmd`) read only `T/data/...` (via `read_tab`/`data_path` in `notebooks/R/common.R`). After
  the cleanup, all 67 literal `read_tab` paths exist. `$P/work/...` appears in the notebooks only in prose.
* 243 concrete `$P/data|work/...` paths are referenced in `pipeline/`, `evidence/`, `regression/`, `notebooks/`,
  `genome_annotation/` and `analyses/`. After the cleanup 31 of them do not exist:
  * 12 were deleted below on purpose (all regenerable): GM12878 short-read FASTQs, STAR BAM and index; the killed de novo
    Trinity run; Arabidopsis real reads, simulated FASTQ, STAR BAM and index; GM12878 simulated FASTQ and STAR BAM.
  * The other 19 are file-name prefixes (`labels/sim_trinity`, `sprot_minus_`, `logs/final_` ...), `dN` placeholders,
    temp dirs (`sorttmp`, `tmp/sort`, `*_tmp`) or de novo Trinity outputs that were never produced. None is caused by
    this cleanup.
* Untouched: `work/bench_v3` (its run dirs, LongOrfs dirs and homology hits are reused by pytd_v7h runs, the Kozak figure
  (E09), the stringency sweeps and step-4 error anatomy), `work/gm12878_td` (TD predictions, eval, ribo, proteomics),
  `work/start_classifier` (step-3 and athal_depth searches), `data/bench_v3`, `data/homology`, `data/annot`, and
  `data/gm12878/proteome_sinitcyn2023` (6.4 GB, MS search results).

## Deleted

| path (pattern) | freed (du) | what it was | kept | regenerate |
|---|---|---|---|---|
| `$P/_trash/*` | ~37 GB | killed or failed E33 step-4/5 runs, simtrinity test dirs, superseded internals | nothing | n/a (already discarded) |
| `data/start_classifier/step5/db/uniref50_2026_03.dmnd`, `db/uniref50_{rep_taxonmap,seqinfo,xml_rep}*.gz`, `taxonomy/` (NCBI taxdump), empty `dl/`, `tmp/` | 15.4 GB | UniRef50 2026_03 DIAMOND db and its build intermediates (E33 step 5; finished, negative) | `db/own_*.tsv`, `db/sprot_seqinfo.tsv.gz`, `queries/`, `searches/`, `hits/`, `features/`, `pilot/`, `timing/`, `variants/`, all logs (0.55 GB). The result tables are in `T/data/start_classifier/step5` | `pipeline/start_classifier/step5_01_uniref50_build_v3.sh` (downloads the UniRef50 XML/FASTA and the taxdump; note that UniRef50 2026_03 may no longer be current on the FTP) |
| `data/gm12878/trinity_denovo_killed/` | 63.9 GB | abandoned de novo Trinity run (E11: "not used") | nothing | `genome_annotation/scripts/data_prep/trinity_denovo.sh` (projected to take weeks) |
| `data/gm12878/trinity/trinity_gg_out/*` except `Trinity-GG.fasta`, `Trinity-GG.fasta.gene_trans_map`, `Trinity.timing`, `FailedCommands` | 61.1 GB | genome-guided Trinity working files: `Dir_*` partitions, norm_200 BAM, per-strand SAM, frag_coverage.wig, read_coords, cmds | assembly (`Trinity.fasta` symlink still resolves), `Trinity.minimap2.bam(.bai)`, flagstat, minimap2.log | `genome_annotation/scripts/data_prep/trinity_gg.sh` (needs the STAR BAM below) |
| `data/gm12878/star_stringtie/{GM12878_ENCSR000AED_rep1.STAR.bam(.bai), star_index/, star_align_tmp/, *._STARpass1, *._STARgenome}` | 52.8 GB (index 30.8, BAM and tmp ~22) | STAR GENCODE v50 index and the GM12878 short-read alignment | `stringtie.gtf`, `stringtie.transcripts.fa`, `Log.final.out` (read by `pipeline/derive/06`), `SJ.out.tab`, flagstat, logs | `genome_annotation/scripts/data_prep/star_index.sh`, then `star_align_stringtie.sh` (StringTie output is kept, so it does not need to be rerun) |
| `data/gm12878/reads/short/ENCFF001REK.fastq.gz`, `ENCFF001REJ.fastq.gz` | 17.3 GB | ENCODE ENCSR000AED rep 1 PE101 reads | long reads `reads/long/ENCFF694DIE.fastq.gz` (1.3 GB) | `genome_annotation/scripts/data_prep/download.sh`; URLs and md5s are in `data/gm12878/README.md` §1a |
| `data/gm12878/riboseq_cenik2015/{sra/, fastq/, trim/, tmp/, star/SRR180314{7,8}.Aligned.sortedByCoord.out.bam, star/*._STARtmp}` | 24.9 GB | Cenik 2015 Ribo-seq downloads, trimmed and rRNA-depleted reads, per-run STAR BAMs | merged `footprints.genome.bam(.bai)` (used for ORF support), `logs/*` and `star/*.Log.final.out` (read by `pipeline/derive/07_gm12878_evidence_v7.sh`), `rrna/` index | `riboseq_cenik2015/download.sh`, then `genome_annotation/scripts/gm12878_riboseq_process.sh` (needs the STAR index above) |
| `data/gm12878/longread_lraa/{lraa_denovo/GM12878_ENCSR838WFC_rep2.LRAA.chunked_work, tmp}` | 4.0 GB | LRAA chunk working dir | `LRAA.gtf` and the other LRAA outputs, `ENCFF694DIE.minimap2.bam` | `genome_annotation/scripts/data_prep/lraa_denovo.sh` |
| `data/simtrinity/trinity/trinity_gg_out/*` (same keep list as for GM12878) | 39.6 GB | E30 genome-guided Trinity working files | `Trinity-GG.fasta`, `Trinity.minimap2.bam`, flagstat, labels and inputs | `pipeline/simtrinity/08_trinity_gg.sh` |
| `data/simtrinity/depth/d{2,4,8,16,32}/trinity_gg_out/*` (same keep list) and `star.sub.bam(.bai)` | 159.8 GB | E31 per-depth Trinity working files and subsampled STAR BAMs | per level: `Trinity-GG.fasta` (`Trinity.fasta` symlink), `Trinity.minimap2.bam`, `labels/`, `inputs/`, `pair_counts.tsv`, minimap2.log | `pipeline/simtrinity/11_depth_subsample.sh <level> <frac>` (samtools `--subsample-seed 20261026`), then `08b_trinity_gg_level.sh` and `12_depth_label_and_run.sh` |
| `data/simtrinity/star/{sim_GM12878.STAR.bam(.bai), star_tmp, _STARpass1, _STARgenome}` | 3.0 GB | STAR alignment of the simulated reads | `SJ.out.tab`, `Log.final.out`, flagstat | `pipeline/simtrinity/07_star.sh` (needs the GM12878 STAR index and the simulated FASTQs) |
| `data/simtrinity/sim/{fastq/, polyester_chunk_*/}` | 4.8 GB | Polyester simulated PE reads | `chunk_*.fa`, `chunk_*.counts.tsv`, `counts.tsv`, `sim_params.tsv`, logs | `pipeline/simtrinity/06_simulate.sh` (`04_polyester_chunk.R` + `05_polyester_to_fastq.py`; chunk k seed = 20261025 + k) |
| `data/simtrinity/salmon/index/` | 1.5 GB | salmon index | `salmon/quant/quant.sf` (used for labels) | `pipeline/simtrinity/02_salmon.sh` (its quant step also needs the GM12878 short reads) |
| `data/simtrinity_athal/depth/d{1,2,4,8,16,32}/trinity_gg_out/*` (same keep list) and `star.sub.bam(.bai)` (d1: symlink) | 30.7 GB | E32 per-depth Trinity working files and subsampled BAMs | as for E31 per level | `pipeline/simtrinity_athal/11_depth_subsample_athal.sh` (seed 20261026), `08_trinity_gg_athal.sh`, `12_label_and_run_athal.sh` / `run_level_athal.sh` |
| `data/simtrinity_athal/reads/ERR4140056_{1,2}.fastq.gz` | 4.8 GB | real Arabidopsis Col-0 leaf reads (ENA ERR4140056), used only for salmon | `salmon/quant/quant.sf` | re-download ERR4140056 from ENA, then `pipeline/simtrinity_athal/02_salmon_athal.sh` |
| `data/simtrinity_athal/sim/{fastq/, polyester_chunk_*/}` | 2.3 GB | Polyester simulated reads | `chunk_*.fa`, counts, `sim_params.tsv`, `fld.tsv`, logs | `pipeline/simtrinity_athal/06_simulate_athal.sh` (chunk k seed = 20261027 + k) |
| `data/simtrinity_athal/star/{sim_athal.STAR.bam(.bai), star_tmp, _STARpass1, _STARgenome}`, `ref/star_index/`, `salmon/index/` | 3.3 GB | STAR BAM, TAIR12 STAR index, salmon index | STAR logs, `SJ.out.tab`, flagstat, `quant.sf` | `03_star_index_athal.sh`, `07_star_athal.sh`, `02_salmon_athal.sh` |
| `work/start_param_grid/{runs,runs_bits,gm12878}/**`: `*.transdecoder_dir/`, `*.transdecoder.{cds,pep,bed}` | 31.2 GB | E13 start-parameter grid (1,264 runs), all already tabulated in `T/data/start_param_sensitivity/grid_scores.tsv` and related tables | per run: `preds.tsv`, `score.tsv`, `runtime_sec.txt`, `predict.log`, `*.transdecoder.gff3`; GM12878 runs also keep `eval.*.tsv`, `*.genome.{gff3,bed}`, `genome_map.log`. `summarise.py` and `open_end_mix.py` can still be rerun | `T/analyses/start_param_sensitivity/run_grid.py`, `run_gm12878.sh` |
| `work/bench_v2/runs/**`: `*.transdecoder_dir/` (except `pytd_norefine`), `*.transdecoder.{cds,pep,bed}`, CodAn `codan_out/`, GeneMarkS-T `seq_bin_N` temp files | 29.5 GB | superseded benchmark v2 (RefSeq truth; E07: not to be quoted) | per run: `preds.tsv`, `*.transdecoder.gff3`, tool GFFs, logs, homology hits; `pytd_norefine/input.fa.transdecoder_dir` (read by `analyses/start_refinement_stability/run_strategy_grid.sh`, E09); `summary_partial.tsv` | `benchmark_v2/scripts/run_tool.sh <tool> <input> <outdir> <org>` on `data/bench_v2/inputs` |

Total by `du`: about 587 GB (data/gm12878 224.0, data/simtrinity 208.7, data/simtrinity_athal 41.4, `_trash` ~37,
work/start_param_grid 31.2, work/bench_v2 29.5, step5 UniRef50 15.4).

Cost of regeneration: rerunning anything that starts from STAR (the simulated Trinity levels, GM12878 Trinity-GG, Ribo-seq)
first needs the GM12878 STAR index rebuilt (`star_index.sh`). The Polyester reads are seeded (seeds above), and the
depth subsamples are hash-seeded and nested. That a rerun reproduces the same assemblies byte for byte has not been
verified, so the kept `Trinity-GG.fasta`, labels and scores are the reference copies.
