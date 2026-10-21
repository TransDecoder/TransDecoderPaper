#!/usr/bin/env python3
"""True CDSs that fail the coding test but are retained by the long-ORF rule: how many also have a
homology hit (Swiss-Prot minus the organism's own entries, or Pfam), i.e. would also be retained by the
homology fallback, and how many only by length.

usage: length_vs_homology_fallback.py [results_dir]   (default: results/)
Reads results/per_cds/<org>.per_cds.tsv.gz (complete inputs); writes results/length_vs_homology_fallback.tsv
"""
import glob
import sys
from pathlib import Path

import pandas as pd

res = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results"


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
