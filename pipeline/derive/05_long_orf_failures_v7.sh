#!/bin/bash
# Coding-test failures of true CDSs under the release rule (pytd_v7 runs; evidence E12, E25).
# Full outputs: $P/work/long_orf_failures_v7/ ; compact tables are then copied to
# data/long_orf_failures/v7/ by 05_fig5_tables_v7.py (last step). 2 processes; about 1-3 h on a shared host.
set -euo pipefail
P=${TD_PAPER_ROOT:-/broad/hptmp/bhaas/TRANSDECODER_PAPER}
PY=$P/env/td/bin/python
D=$P/TransDecoderPaper/pipeline/derive
OUT=$P/work/long_orf_failures_v7
mkdir -p $OUT
$PY $D/05_long_orf_failures_v7.py bench --input complete --jobs 2 \
    --orgs tthermophila human mouse creinhardtii athal dmel spombe pfalciparum --out $OUT
$PY $D/05_long_orf_failures_v7.py bench --input complete_fullutr --jobs 2 --out $OUT
$PY $D/05_long_orf_failures_v7.py altrules --jobs 2 --out $OUT
$PY $D/05_long_orf_failures_v7.py summarize --out $OUT
$PY $D/05_length_vs_homology_v7.py $OUT
$PY $D/05_fig5_tables_v7.py $OUT
