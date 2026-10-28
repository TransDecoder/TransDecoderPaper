#!/usr/bin/env python3
"""Release predictions on GENCODE lncRNA negatives that are at least the long-ORF auto-accept length.

ORFs at least L*(GC) long are accepted on length alone, so these calls are not controlled by the q-value
test, whether or not they would also pass it. The count is an upper bound on calls accepted by length
only. Evidence: E25 (section "lncRNA calls at or above L*").

Method:
  - For human and mouse, complete.wLnc and fragments.wLnc, take the release run
    <runs>/<org>/<variant>/pytd_v7/.
  - Read L* from its predict.log ("min auto-accept length = <nt> nt").
  - Count the predictions in preds.tsv on lncRNA negatives (transcript IDs starting with "lnc_"), and
    those whose length, end - start + 1 (nt, stop codon included), is >= L*.

usage: 04_lnc_calls_by_autoaccept_length.py [P]
writes $T/data/benchmark/final_v7/lnc_calls_by_autoaccept_length.tsv
"""
import re
import sys
from pathlib import Path

P = Path(sys.argv[1] if len(sys.argv) > 1 else "/broad/hptmp/bhaas/TRANSDECODER_PAPER")
RUNS = P / "work/bench_v3/runs"
OUT = P / "TransDecoderPaper/data/benchmark/final_v7/lnc_calls_by_autoaccept_length.tsv"

rows = []
for org in ("human", "mouse"):
    for variant in ("complete.wLnc", "fragments.wLnc"):
        run = RUNS / org / variant / "pytd_v7"
        m = re.search(r"min auto-accept length = (\d+) nt", (run / "predict.log").read_text())
        if not m:
            sys.exit(f"no auto-accept length in {run}/predict.log")
        lstar = int(m.group(1))
        n = n_long = 0
        with open(run / "preds.tsv") as fh:
            next(fh)
            for line in fh:
                tid, start, end = line.split("\t")[:3]
                if tid.startswith("lnc_"):
                    n += 1
                    n_long += int(end) - int(start) + 1 >= lstar
        rows.append((org, variant, lstar, n, n_long, n_long / n))

with open(OUT, "w") as fh:
    fh.write("org\tvariant\tL_star_nt\tlnc_calls\tlnc_calls_ge_L_star\tfrac_ge_L_star\n")
    for r in rows:
        fh.write("\t".join(map(str, r[:5])) + f"\t{r[5]:.4f}\n")
print(OUT.read_text(), end="")
