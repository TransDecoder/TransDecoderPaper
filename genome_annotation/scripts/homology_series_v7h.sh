#!/bin/bash
# E06 release rerun: one configuration of the GM12878 homology series with PyTD a5b24c3 (PY_hstart;
# homology-assisted start refinement on, the default when hits are supplied). Fresh run directory with a
# copy of the per-assembly fresh LongOrfs (homology_series_v7h_prep.sh), then Predict, genome-map and
# eval_genome_cds.py as homology_series.sh (E06) and gm12878_final_v7.sh (E25).
# usage: homology_series_v7h.sh <lraa|stringtie|trinity> <clade taxid|none> <pfam|nopfam>
#   none nopfam = ab initio (no hits); none pfam = Pfam only.
set -euo pipefail
name=$1; clade=$2; pf=$3
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; D=$P/data/gm12878; W=$P/work/gm12878_td
PY=$P/env/td/bin/python; export PATH=$P/env/td/bin:$PATH; export PYTHONPATH=$P/work/trees/PY_hstart
declare -A S=([lraa]=$D/longread_lraa/lraa_denovo/LRAA.gtf [stringtie]=$D/star_stringtie/stringtie.gtf
              [trinity]=$D/trinity/Trinity.minimap2.bam)
structs=$(readlink -f ${S[$name]})
tag=blast_${clade}_$pf; [ $clade = none ] && tag=$pf; [ $clade = none ] && [ $pf = nopfam ] && tag=abinitio
base=$W/$name/hseries_v7h; out=$base/$tag
[ -e $out ] && { echo "exists: $out" >&2; exit 1; }
mkdir -p $out; cd $out
(cd $PYTHONPATH && git rev-parse --short HEAD) > pytd_version.txt
[ "$(cat pytd_version.txt)" = a5b24c3 ] || { echo "PyTD is not a5b24c3" >&2; exit 1; }
ln -s $(readlink -f $base/_longorfs/transcripts.fa) transcripts.fa
rsync -a --exclude '__checkpoints*' $base/_longorfs/transcripts.fa.transdecoder_dir/ transcripts.fa.transdecoder_dir/
md5sum -c --quiet $base/longest_orfs.pep.md5
extra=()
[ $clade != none ] && extra+=(--retain-blastp-hits $base/blastp.minus_$clade.outfmt6)
[ $pf = pfam ] && extra+=(--retain-pfam-hits $base/pfam.domtblout)
td() { $PY -m pytransdecoder.cli "$@"; }
td predict -t transcripts.fa ${extra[@]+"${extra[@]}"} > predict.log 2>&1
td genome-map --orfs transcripts.fa.transdecoder.gff3 --structures $structs -t transcripts.fa \
   -o transcripts.fa.transdecoder.genome.gff3 > genome_map.log 2>&1
mkdir -p $W/eval/hseries_v7h
$PY $P/TransDecoderPaper/genome_annotation/scripts/eval_genome_cds.py --pred transcripts.fa.transdecoder.genome.gff3 \
   --structures $structs --mane $P/data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
   --gencode $D/ref/gencode.v50.primary_assembly.annotation.gtf --label ${name}_$tag \
   -o $W/eval/hseries_v7h/${name}_$tag > eval.log 2>&1
echo "DONE $name $tag"
