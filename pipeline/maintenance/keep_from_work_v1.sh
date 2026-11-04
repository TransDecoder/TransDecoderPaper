#!/bin/bash
# Copy the parts of $P/work that are inputs, expensive to regenerate, or the only record, into
# $P/data/work_keep/<same path relative to $P/work>, before $P/work is deleted (BH 2026-11-03).
# What and why: REPRODUCE.md, section "What was kept from work/". Copies only; never deletes.
# usage: keep_from_work_v1.sh [--dry-run]
set -u   # no -e: a missing optional path must not stop the list; rsync failure is checked below
P=/broad/hptmp/bhaas/TRANSDECODER_PAPER; W=$P/work; K=$P/data/work_keep
dry=(); [ "${1:-}" = --dry-run ] && dry=(--dry-run)
mkdir -p $K
cd $W
list=$(mktemp)
{
  # 1. homology hits (DIAMOND vs sprot_minus_<taxid>, hmmsearch Pfam-A 38.2) and the peptides they were run on
  find bench_v3/runs -path '*/pytd_v7_homology/*' \( -name blastp.outfmt6 -o -name pfam.domtblout -o -name longest_orfs.pep \)
  find gm12878_td -maxdepth 4 \( -name 'blastp*.outfmt6' -o -name pfam.domtblout \) -not -path '*/pfam_chunks/*'
  find start_classifier -path '*/searches/*' -type f
  [ -d start_classifier/homology ] && find start_classifier/homology -type f
  # 2. standardized predictions + run times of every benchmark run (incl. CodAn 1,629 CPU-h, licensed GeneMarkS-T,
  #    unreleased-code rules such as cliff_c200, historical releases)
  find bench_v3/runs bench_v3/rules bench_v3/history bench_v3/rerun_old16 bench_v3/frames bench_v3/stringency bench_v3/minlen \
       bench_v3/presets bench_v3/first_atg_v7 simtrinity simtrinity_athal \
       \( -name 'preds*.tsv' -o -name 'longest_orf.preds.tsv' -o -name runtime_sec.txt -o -name '*version*.txt' \) 2>/dev/null
  # 3. inputs that cannot be regenerated or exist only here
  find bench_v3/presets/ref bench_v3/presets/tests -type f
  find history/sourceforge history/sf_runs bench_v3/scores -type f
  find lnc -maxdepth 1 \( -name '*.pick.*' -o -name '*.py' \)
  find thresh_big/raw -type f
  find utr_prior/annot -type f; find utr_prior/gm12878 -maxdepth 1 -type f \( -name '*truth*' -o -name '*pairs*' \)
  find simtrinity/realtrinity_label simtrinity/filtered_preds -type f
  find simtrinity/depth simtrinity_athal/depth -path '*/scores/*' -type f
  find analysis/pwm_audit analysis/kozak start_grid -type f 2>/dev/null; ls analysis/fork_findings.md
  # 4. release-pin per-run small files (tables can be re-derived without rerunning)
  find release_pin \( -name 'preds*.tsv' -o -name 'score*.tsv' -o -name 'eval*' -o -name 'ribo*.tsv' -o -name runtime_sec.txt \
       -o -name pytd_version.txt -o -name '*.kozak.tsv' -o -name '*.log' -o -name 'preset.json' \) -type f -size -50M
  # 5. every script, note and launcher (several launchers were never committed)
  find . \( -name '*.sh' -o -name '*.py' -o -name '*.R' -o -name '*.awk' -o -name '*.md' -o -name 'jobs*.txt' \) -type f -size -5M
} 2>/dev/null | sed 's#^\./##' | sort -u > $list
echo "$(wc -l < $list) files"
rsync -a ${dry[@]+"${dry[@]}"} --files-from=$list $W/ $K/ || { echo "rsync failed" >&2; exit 1; }
du -sh $K
rm -f $list
