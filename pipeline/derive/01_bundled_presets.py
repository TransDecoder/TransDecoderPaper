#!/usr/bin/env python3
"""Summarise the seven bundled presets of the release (Table tab:presets; evidence E21, E24).

usage: 01_bundled_presets.py [P]      (P = /broad/hptmp/bhaas/TRANSDECODER_PAPER)
Reads the preset directories shipped with PyTransDecoder main 3a793ff (work/trees/PY_v7/pytransdecoder/presets)
and TransDecoder Perl devel 41fc5ee (work/trees/TD_v7/PerlLib/presets). Writes
data/benchmark/presets/bundled_presets.tsv: one row per preset with preset.json fields, counts of
q levels, strata and training ORFs, and whether every file is byte-identical between the two packages.
Run with $P/env/td/bin/python.
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

P = Path(sys.argv[1] if len(sys.argv) > 1 else "/broad/hptmp/bhaas/TRANSDECODER_PAPER")
PY = P / "work/trees/PY_v7/pytransdecoder/presets"
PL = P / "work/trees/TD_v7/PerlLib/presets"
OUT = P / "TransDecoderPaper/data/benchmark/presets/bundled_presets.tsv"
ORDER = ["mammal", "insect", "plant", "fungus", "apicomplexan", "ciliate", "green_alga"]


def md5(f):
    return hashlib.md5(f.read_bytes()).hexdigest()


def gz_lines(f):
    with gzip.open(f, "rt") as fh:
        return [l for l in fh]


cols = ["name", "source", "genetic_code", "gc_pct", "auto_accept_nt", "n_q_levels", "q_min", "q_max",
        "n_strata", "n_threshold_rows", "n_training_orfs", "n_files", "perl_py_identical"]
rows = []
for name in ORDER:
    d = PY / name
    meta = json.loads((d / "preset.json").read_text())
    thr = [l.rstrip("\n").split("\t") for l in gz_lines(d / "thresholds.tsv.gz") if not l.startswith("#")]
    strata = {r[0] for r in thr}
    n_train = sum(1 for l in gz_lines(d / "training_orfs.cds.gz") if l.startswith(">"))
    files = sorted(p.name for p in d.iterdir())
    same = all((PL / name / f).exists() and md5(d / f) == md5(PL / name / f) for f in files) and \
        sorted(p.name for p in (PL / name).iterdir()) == files
    q = meta["q_levels"]
    aa = meta.get("auto_accept_nt")
    rows.append([name, meta.get("source", ""), meta["genetic_code"], meta["gc_pct"],
                 "none" if aa in (None, 0) else aa, len(q), min(q), max(q), len(strata), len(thr),
                 n_train, len(files), same])

with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
print(OUT)
