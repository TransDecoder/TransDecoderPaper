#!/bin/bash
# Homology-assisted PyTD at the final defaults (PY_dev 48c870c): DIAMOND vs Swiss-Prot minus the
# organism + Pfam (chunked hmmsearch), all inputs. Replaces pytd_homology (run at 27bfb5e).
cd /broad/hptmp/bhaas/TRANSDECODER_PAPER
JOBS=3 THREADS=4 TransDecoderPaper/benchmark_v2/scripts/run_all.sh data/bench_v3/inputs work/bench_v3/runs pytd_dev_v2_homology
echo ALL_DONE17
