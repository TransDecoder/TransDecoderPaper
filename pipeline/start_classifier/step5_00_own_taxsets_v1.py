#!/usr/bin/env python3
"""E33 step 5: for the 8 benchmark organisms, list each taxid and its NCBI descendants.
Usage: step5_00_own_taxsets_v1.py nodes.dmp > own_taxsets.tsv   (columns: taxid, organism_taxid)
These are the taxa that `diamond --taxon-exclude <org>` removes."""
import sys
ORGS = ['9606', '10090', '3702', '7227', '284812', '3055', '36329', '312017']
ch = {}
for line in open(sys.argv[1]):
    f = line.split('\t|\t')
    ch.setdefault(f[1], []).append(f[0])
print('taxid\torg_taxid')
for o in ORGS:
    st = [o]
    seen = set()
    while st:
        x = st.pop()
        if x in seen:
            continue
        seen.add(x)
        print(f'{x}\t{o}')
        st.extend(c for c in ch.get(x, []) if c != x)
