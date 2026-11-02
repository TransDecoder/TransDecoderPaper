#!/usr/bin/env python3
"""GM12878 tables from the pinned-commit runs (E40; work/release_pin/gm/<asm>/<cfg>).
Writes to data/gm12878/release_pin/:
  gm_runs.tsv        one row per assembly x config: all eval_genome_cds.py summary metrics, prediction-side counts of
                     complete predictions (as E13: exact / compatible not exact / start share), and the
                     catalogued Ribo-seq ORF matches of match_riboseq.py
  no_overlap_by_length.tsv   as pipeline/derive/08_no_overlap_by_length_v7.sh, rel only
usage: 12_gm_tables_v1.py
"""
import collections
import csv
import os
import re

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
GM = f"{P}/work/release_pin/gm"
OUT = f"{P}/TransDecoderPaper/data/gm12878/release_pin"
os.makedirs(OUT, exist_ok=True)
rows = []
keys = None
for asm in ("lraa", "stringtie", "trinity"):
    for cfg in sorted(os.listdir(f"{GM}/{asm}")):
        d = f"{GM}/{asm}/{cfg}"
        if not os.path.exists(f"{d}/DONE"):
            continue
        kv = dict(line.rstrip("\n").split("\t") for line in open(f"{d}/eval.summary.tsv"))
        kv.pop("label", None)
        m = re.match(r"^a([0-9.]+)_w([0-9.]+)$", cfg)
        r = dict(assembly=asm, config=cfg, alpha=m.group(1) if m else "NA", weight=m.group(2) if m else "NA")
        r.update(kv)
        pc = collections.Counter()
        with open(f"{d}/eval.pred_status.tsv") as fh:
            for x in csv.DictReader(fh, delimiter="\t"):
                if x["gencode_compatible"] != "1":
                    continue
                if x["orf_type"] == "complete":
                    pc["e" if x["gencode_exact"] == "1" else "n"] += 1
                else:
                    pc["p"] += 1
        r.update(pred_complete_gencode_exact=pc["e"], pred_complete_compatible_not_exact=pc["n"],
                 pred_partial_gencode_compatible=pc["p"],
                 pred_complete_start_precision=f"{pc['e'] / max(1, pc['e'] + pc['n']):.4f}")
        n = sm = em = 0
        with open(f"{d}/ribo.tsv") as fh:
            for x in csv.DictReader(fh, delimiter="\t"):
                n += 1
                sm += x["stop_match"] == "1"
                em += x["exact_match"] == "1"
        r.update(ribo_catalog_models=n, ribo_catalog_stop_match=sm, ribo_catalog_exact_match=em)
        rows.append(r)
cols = list(rows[0])
for r in rows:
    cols += [k for k in r if k not in cols]
with open(f"{OUT}/gm_runs.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
with open(f"{OUT}/no_overlap_by_length.tsv", "w") as fh:
    fh.write("assembly\tconfig\tlength_bin\tn\tno_gencode_cds_overlap\n")
    for asm in ("lraa", "stringtie", "trinity"):
        for cfg in ("rel", "m50_noq"):
            f = f"{GM}/{asm}/{cfg}/eval.pred_status.tsv"
            if not os.path.exists(f):
                continue
            n, o = collections.Counter(), collections.Counter()
            for x in csv.DictReader(open(f), delimiter="\t"):
                L = int(x["aa"])
                b = "<100" if L < 100 else "100-149" if L < 150 else "150-199" if L < 200 else "200-299" if L < 300 else ">=300"
                n[b] += 1
                o[b] += int(x["no_gencode_cds_overlap"])
            for b in ("<100", "100-149", "150-199", "200-299", ">=300"):
                fh.write(f"{asm}\t{cfg}\t{b}\t{n[b]}\t{o[b]}\n")
print(len(rows), "rows")
