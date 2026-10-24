#!/bin/bash
# fresh PyTD (release branch) run on a bundled-preset reference, then build-preset
set -euo pipefail
clade=$1; gcode=$2; org=$3
R=/broad/hptmp/bhaas/TRANSDECODER_PAPER; B=$R/work/bench_v3/presets/bundled
export PYTHONPATH=$R/work/trees/PY_taper; PYB=$R/env/td/bin/python
out=$B/runs/$clade
if [ -e $out ]; then mkdir -p $B/runs/_trash; mv $out $B/runs/_trash/$clade.$(date +%s); fi
mkdir -p $out; cd $out; ln -s $B/ref/$clade.fa $clade.fa
(cd $R/work/trees/PY_taper && git rev-parse --short HEAD) > pytd_version.txt
$PYB -m pytransdecoder.cli longorfs -t $clade.fa -G $gcode -S > longorfs.log 2>&1
$PYB -m pytransdecoder.cli predict -t $clade.fa -G $gcode > predict.log 2>&1
$PYB -m pytransdecoder.cli build-preset -t $clade.fa --name $clade -o $B/presets/$clade \
  --source "$org benchmark annotation (TransDecoder paper 2026): complete coding transcripts + their fragments$( [ $org = human ] && echo ' + GENCODE lncRNAs')" > build.log 2>&1
echo done
