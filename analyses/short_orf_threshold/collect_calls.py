#!/usr/bin/env python3
"""Per-call table (score, length, label) for PyTD runs, to explore score thresholds for short ORFs.

usage: collect_calls.py <inputs_dir> <runs_dir> <tool> <variant> [<variant> ...] > calls.tsv
label: tp (stop end within 2 nt of the truth on a coding transcript), wrong (other call on a
coding transcript), lnc (on a GENCODE lncRNA 'lnc_*'), random (on 'rand_*').
Also the 6 frame scores of the ORF from longest_orfs.cds.scores.
"""
import os, re, sys
inputs, runs, tool, variants = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
print("org\tvariant\tmodel_id\torf_type\taa\tscore\tbest_other_frame\tlabel")
for org in sorted(os.listdir(runs)):
    for v in variants:
        rd = os.path.join(runs, org, v, tool)
        g = os.path.join(rd, "input.fa.transdecoder.gff3")
        if not os.path.exists(g):
            continue
        truth = {}
        with open(os.path.join(inputs, org, v + ".truth.tsv")) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for l in fh:
                r = dict(zip(hdr, l.rstrip("\n").split("\t")))
                truth[r["transcript_id"]] = (int(r["cds_start"]), int(r["cds_end"]), r["strand"])
        fs = {}
        for l in open(os.path.join(rd, "input.fa.transdecoder_dir", "longest_orfs.cds.scores")):
            if l[0] != "#":
                x = l.split("\t"); fs[x[0]] = [float(y) for y in x[3:9]]
        cds = {}
        for l in open(g):
            f = l.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "CDS":
                cds[re.search(r"Parent=([^;]+)", f[8]).group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
        for l in open(g):
            f = l.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "mRNA":
                continue
            mid = re.search(r"ID=([^;]+)", f[8]).group(1)
            typ = re.search(r"ORF type:(\w+)", f[8]).group(1)
            sc = float(re.search(r";score=([-\d.]+)", f[8]).group(1))
            tid, a, b, s = cds[mid]
            aa = (b - a + 1) // 3
            fr = fs.get(mid, [sc] + [0] * 5)
            if tid.startswith("lnc_"): lab = "lnc"
            elif tid.startswith("rand_"): lab = "random"
            else:
                t = truth.get(tid)
                p3 = b if s == "+" else a
                t3 = None if t is None else (t[1] if t[2] == "+" else t[0])
                lab = "tp" if t and t[2] == s and abs(p3 - t3) <= 2 else "wrong"
            print(f"{org}\t{v}\t{mid}\t{typ}\t{aa}\t{sc}\t{max(fr[1], fr[2]):.2f}\t{lab}")
