#!/bin/bash
# One fresh PyTD run on one GM12878 assembly with the pinned release code (E40): predict, genome-map,
# eval_genome_cds.py, match_riboseq.py. Procedure as genome_annotation/scripts/gm12878_final_v7.sh.
# usage: run_gm_v1.sh <lraa|stringtie|trinity> <cfg>
# cfg: rel | precise | sensitive | a<A>_w<W> | norefine | m<N>_noq | retonly
# Configs other than rel and m<N>_noq reuse rel's LongOrfs output (same commit, same settings), symlinked.
set -euo pipefail
asm=$1; cfg=$2
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; T=$P/TransDecoderPaper
PYTD=${PYTD_REL:-$P/work/trees/PY_rel_d3892b5}; PY=$P/env/td/bin/python
export PATH=$P/env/td/bin:$PATH OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=$PYTD
case $asm in
  stringtie) fa=$D/star_stringtie/stringtie.transcripts.fa; st=$D/star_stringtie/stringtie.gtf;;
  lraa) fa=$D/longread_lraa/lraa_denovo/LRAA.transcripts.fa; st=$D/longread_lraa/lraa_denovo/LRAA.gtf;;
  trinity) fa=$D/trinity/Trinity.fasta; st=$D/trinity/Trinity.minimap2.bam;;
esac
out=$P/work/release_pin/gm/$asm/$cfg
[ -s $out/DONE ] && { echo "done already: $out"; exit 0; }
[ -e $out ] && { mkdir -p $P/_trash/release_pin_partial; mv $out $P/_trash/release_pin_partial/gm.$asm.$cfg.$(date +%s); }
mkdir -p $out; cd $out
ln -s $(readlink -f $fa) transcripts.fa
(cd $PYTD && git rev-parse --short HEAD) > pytd_version.txt
td() { $PY -m pytransdecoder.cli "$@"; }
t0=$(date +%s)
pr=()
case $cfg in
  rel) td longorfs -t transcripts.fa -S > longorfs.log 2>&1;;
  m*_noq) n=${cfg#m}; td longorfs -t transcripts.fa -S -m ${n%_noq} > longorfs.log 2>&1; pr=(--max-qvalue 1000000);;
  *) src=$P/work/release_pin/gm/$asm/rel/transcripts.fa.transdecoder_dir
     [ -s $src/longest_orfs.pep ] || { echo "no rel LongOrfs" >&2; exit 1; }
     mkdir transcripts.fa.transdecoder_dir
     for f in longest_orfs.cds longest_orfs.gff3 longest_orfs.pep base_freqs.dat longorfs.params.json; do ln -s $src/$f transcripts.fa.transdecoder_dir/$f; done
     case $cfg in
       precise|sensitive) pr=(--$cfg);;
       norefine) pr=(--no-refine-starts);;
       a*_w*) a=${cfg#a}; a=${a%_w*}; w=${cfg#*_w}; pr=(--start-contrast-alpha $a --start-kozak-weight $w);;
       retonly) h=$P/work/gm12878_td/$asm/final_v7_homology
          cmp transcripts.fa.transdecoder_dir/longest_orfs.pep $h/transcripts.fa.transdecoder_dir/longest_orfs.pep
          pr=(--retain-blastp-hits $h/blastp.outfmt6 --retain-pfam-hits $h/pfam.domtblout --no-homology-start-assist);;
       *) echo "bad cfg $cfg" >&2; exit 1;;
     esac;;
esac
td predict -t transcripts.fa ${pr[@]+"${pr[@]}"} > predict.log 2>&1
td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $st -t transcripts.fa \
   -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
$PY $T/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
   --structures $st --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
   --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${asm}_$cfg -o eval > eval.log 2>&1
PYTHONPATH=$T/genome_annotation/scripts:$PYTHONPATH $PY $T/analyses/riboseq/match_riboseq.py \
   --pred transcripts.fa.transdecoder.genome.gff3 --phase1 $P/data/riboseq_orfs/phase1/Ribo-seq_ORFs.all.bed \
   --phase2 $P/data/riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.bed > ribo.tsv
echo $(( $(date +%s) - t0 )) > runtime_sec.txt
date > DONE
