# GM12878 data preparation (as run; absolute paths refer to the analysis machine, `$P=/broad/hptmp/bhaas/TRANSDECODER_PAPER`)

| Script | Purpose |
|---|---|
| `annot_download.sh` | MANE / GENCODE reference annotation downloads (`$P/data/annot`) |
| `download.sh` | GM12878 short- and long-read RNA-seq downloads |
| `star_index.sh`, `star_align_stringtie.sh` | STAR index (GENCODE v50 junctions), alignment, StringTie assembly |
| `longread_minimap2.sh`, `pull_lraa.sh`, `lraa_denovo.sh`, `lraa_post.sh` | long-read alignment and LRAA isoform reconstruction |
| `trinity_gg.sh` | Trinity genome-guided assembly (+ minimap2 transcript alignment) |
| `trinity_denovo.sh` | Trinity de novo (not used in the final comparison) |

Ribo-seq processing: `../gm12878_riboseq_*.sh`. TransDecoder runs: `../gm12878_final_defaults.sh`, `../homology_series*.sh`.
