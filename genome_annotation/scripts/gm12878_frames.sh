#!/bin/bash
# GM12878, frame-contrast experiment (E26): PyTD exp/contrast-frames (PY_frames worktree), --contrast-frames <mode>:
# PyTD -> genome-map -> MANE/GENCODE eval -> Ribo-seq match, for one assembly. Fresh output directory
# (an existing one is refused).
# usage: gm12878_frames.sh <name> <transcripts.fa> <structures gtf|bam> <forward|all> [threads] [default|precise|sensitive]
set -euo pipefail
name=$1; fa=$(readlink -f $2); structs=$(readlink -f $3); frames=$4; mode=abinitio; threads=${5:-8}; qmode=${6:-default}
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH; export PYTHONPATH=${PYTD_FRAMES:-$P/work/trees/PY_frames}
tag=frames_$frames; [ $qmode != default ] && tag=${tag}_$qmode
out=$P/work/gm12878_td/$name/final_$tag
[ -e $out ] && { echo "$out exists; refusing to reuse" >&2; exit 1; }
mkdir -p $out; cd $out
(cd $PYTHONPATH && git rev-parse --short HEAD) > pytd_version.txt
qa=(); [ $qmode != default ] && qa=(--$qmode)
ln -sf $fa transcripts.fa
td() { $PY -m pytransdecoder.cli "$@"; }
[ -s transcripts.fa.transdecoder_dir/longest_orfs.pep ] || td longorfs -t transcripts.fa -S > longorfs.log 2>&1
extra=()
if [ $mode = homology ]; then
  pep=transcripts.fa.transdecoder_dir/longest_orfs.pep
  [ -s blastp.outfmt6 ] || diamond blastp --query $pep --db $P/data/homology/sprot_minus_9606.dmnd \
      --max-target-seqs 1 --evalue 1e-5 --outfmt 6 --threads $threads --out blastp.outfmt6 > diamond.log 2>&1
  [ -s pfam.domtblout ] || hmmsearch --cpu $threads --cut_ga --noali -o /dev/null --domtblout pfam.domtblout \
      $P/data/homology/Pfam-A.hmm $pep > hmmsearch.log 2>&1
  extra=(--retain-blastp-hits blastp.outfmt6 --retain-pfam-hits pfam.domtblout)
fi
td predict -t transcripts.fa --contrast-frames $frames ${extra[@]+"${extra[@]}"} ${qa[@]+"${qa[@]}"} > predict.log 2>&1
td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $structs -t transcripts.fa \
   -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
$PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
   --structures $structs --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
   --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${name}_final_$tag \
   -o $P/work/gm12878_td/eval/${name}_final_$tag > eval.log 2>&1
PYTHONPATH=$P/TransDecoderPaper/genome_annotation/scripts:$PYTHONPATH $PY $P/TransDecoderPaper/analyses/riboseq/match_riboseq.py \
   --pred transcripts.fa.transdecoder.genome.gff3 --phase1 $P/data/riboseq_orfs/phase1/Ribo-seq_ORFs.all.bed \
   --phase2 $P/data/riboseq_orfs/phase2/Ribo-seq_ORFs.comprehensive.bed > $P/work/gm12878_td/ribo/${name}_final_$tag.ribo.tsv
echo "DONE $name $tag"
