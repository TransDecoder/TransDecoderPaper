#!/usr/bin/env python3
"""E34 check: a fresh PyTD a5b24c3 Predict --sensitive (02_fresh_runs_v1.sh sensitive) gives the same predictions
as the offline q 0.2 selection + PyTD refinement of 01_trace_missed_v1.py on the release run's intermediates.
usage: python 07_check_sensitive_v1.py <big_dir> <dataset> ...
"""
import re
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
big = sys.argv[1]
for ds in sys.argv[2:]:
    fresh = set()
    for line in open(f"{P}/work/missed_cds/sensitive/{ds}/runs/pytd_v7/preds.tsv"):
        x = line.rstrip("\n").split("\t")
        if x[0] != "transcript_id":
            fresh.add((x[0], int(x[1]), int(x[2]), x[3]))
    off = set()
    for line in open(f"{big}/{ds}/td/q0.2.revised.gff3"):
        x = line.rstrip("\n").split("\t")
        if len(x) == 9 and x[2] == "mRNA":
            off.add((x[0], int(x[3]), int(x[4]), x[6]))
    print(f"{ds}\tfresh_sensitive={len(fresh)}\toffline_q0.2={len(off)}\tonly_fresh={len(fresh - off)}\tonly_offline={len(off - fresh)}")
