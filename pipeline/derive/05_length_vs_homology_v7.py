#!/usr/bin/env python3
"""True CDSs that fail the coding test but are retained by the long-ORF rule: how many also have a
homology hit (Swiss-Prot minus the organism's own entries, or Pfam), i.e. would also be retained by the
homology fallback, and how many only by length.
Copy of analyses/long_orf_failures/length_vs_homology_fallback.py for the release-rule tables
(05_long_orf_failures_v7.py); only the default directory differs. Homology flags now cover ORFs >= 50 aa.

usage: 05_length_vs_homology_v7.py [results_dir]   (default: $P/work/long_orf_failures_v7)
Reads <dir>/per_cds/<org>.per_cds.tsv.gz (complete inputs); writes <dir>/length_vs_homology_fallback.tsv
"""
import glob
import sys
from pathlib import Path

import pandas as pd

res = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/long_orf_failures_v7")


def flag(s):
    return s.astype(str).str.lower().isin(["true", "1", "1.0"])


rows = []
for f in sorted(glob.glob(str(res / "per_cds" / "*.per_cds.tsv.gz"))):
    if "fullutr" in f:
        continue
    d = pd.read_csv(f, sep="\t", low_memory=False)
    fail = d[d.pass_coding.notna() & ~flag(d.pass_coding)]
    resc = fail[flag(fail.long_rule)]
    hit = flag(resc.has_sprot) | flag(resc.has_pfam)
    rows.append(dict(org=Path(f).name.split(".")[0], fail_coding=len(fail), rescued_by_length=len(resc),
                     with_homology_data=int(flag(resc.homology_orf).sum()), also_homology_hit=int(hit.sum()),
                     length_only=int((flag(resc.homology_orf) & ~hit).sum())))
t = pd.DataFrame(rows)
tot = t.drop(columns="org").sum()
t = pd.concat([t, pd.DataFrame([dict(org="ALL", **tot.to_dict())])], ignore_index=True)
t.to_csv(res / "length_vs_homology_fallback.tsv", sep="\t", index=False)
print(t.to_string(index=False))
