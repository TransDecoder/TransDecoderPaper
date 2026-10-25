#!/usr/bin/env python3
"""Per-ORF start-refinement features of one finished PyTransDecoder run (evidence E29).

usage: PYTHONPATH=<PyTD tree> python 03_utr_features.py <run_dir> <genetic_code> <out.tsv.gz> [truth.tsv]

Re-derives, from the run's own files (best_candidates.gff3 = selected ORFs before refinement,
hexamer.scores, top_training_orfs.cds, start_refinement.kozak.tsv), the quantities the release
start refinement compares for every selected ORF that is open at its 5' end:
  d        implied 5' UTR length (nt) if the ORF were trimmed to its first in-frame ATG
           (0-based transcript index of that ATG, in ORF orientation)
  seg      length of the segment between the open end and the ATG
  coding   frame contrast of the segment / mu  (bases' worth of coding signal)
  kozak    Kozak score of the ATG (bits; 0 without a model)
  trim     1 if the release rule trims (alpha*seg + w*kozak > coding)
and checks 'trim' against the run's best_candidates.gff3.revised_starts.gff3.
With a truth table (benchmark inputs), adds the reference CDS ends (transcript orientation),
the reference type and whether the ORF's 3' end matches the reference 3' end (+-2 nt).
"""
import gzip
import re
import sys
from pathlib import Path

from pytransdecoder.predict import TransDecoderPredict, DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT
from pytransdecoder.core import kozak as kozak_mod
from pytransdecoder.core.sequence import reverse_complement

run = Path(sys.argv[1]).resolve()
gcode = sys.argv[2]
out = sys.argv[3]
truth_file = sys.argv[4] if len(sys.argv) > 4 else None
ALPHA, W = DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT

fa = run / "input.fa" if (run / "input.fa").exists() else run / "transcripts.fa"
tdp = TransDecoderPredict(transcripts_file=fa, output_dir=run, genetic_code=gcode)
wd = tdp.workdir
hexfile = wd / "hexamer.scores"
train = wd / "__checkpoints_TDpredict" / "top_training_orfs.cds"
scores = tdp._load_hexamer_scores(hexfile)
mu = tdp._coding_contrast(train, hexfile)
kfile = wd / "start_refinement.kozak.tsv"
kmodel = kozak_mod.read_model(kfile) if kfile.exists() else None
tx = tdp._load_transcripts()

revised = {}
rfile = wd / "best_candidates.gff3.revised_starts.gff3"
if rfile.exists():
    for line in open(rfile):
        p = line.rstrip("\n").split("\t")
        if len(p) == 9 and p[2] == "mRNA":
            revised[re.search(r"ID=([^;]+)", p[8]).group(1)] = (int(p[3]), int(p[4]), "start_revised=true" in p[8])

truth = {}
if truth_file:
    with open(truth_file) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            truth[r["transcript_id"]] = r

cols = ["transcript_id", "orf_id", "strand", "tx_len", "orf5", "orf3", "open5", "d", "seg", "coding",
        "kozak", "utr_score", "trim", "trim_run", "mu", "truth_type", "truth5", "truth3", "stop_match"]
n_mis = 0
with gzip.open(out, "wt") as o:
    o.write("\t".join(cols) + "\n")
    for line in open(wd / "best_candidates.gff3"):
        p = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(p) != 9 or p[2] != "mRNA":
            continue
        tid, s, e, strand = p[0], int(p[3]), int(p[4]), p[6]
        oid = re.search(r"ID=([^;]+)", p[8]).group(1)
        seq = tx[tid]
        L = len(seq)
        if strand == "-":
            oriented = reverse_complement(seq)
            o5 = L - e  # 0-based
        else:
            oriented = seq
            o5 = s - 1
        olen = e - s + 1
        orf = oriented[o5:o5 + olen]
        open5 = int(orf[:3] != "ATG")
        d = seg = coding = kz = us = trim = ""
        if open5:
            atg = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"), None)
            if atg is not None:
                c = tdp._frame_contrast(orf[:atg], scores) / mu
                k = kozak_mod.score(kmodel, oriented, o5 + atg) if kmodel else 0.0
                u = ALPHA * atg + W * k
                d, seg, coding, kz, us, trim = o5 + atg, atg, f"{c:.2f}", f"{k:.3f}", f"{u:.2f}", int(u > c)
        rv = revised.get(oid)
        trim_run = "" if rv is None else int(rv[2])
        if trim != "" and trim_run != "" and trim != trim_run:
            n_mis += 1
        t = truth.get(tid)
        tt = t5 = t3 = sm = ""
        if t:
            ts, te = int(t["cds_start"]), int(t["cds_end"])
            if t["strand"] == "-":
                t5, t3 = L - te + 1, L - ts + 1
            else:
                t5, t3 = ts, te
            tt = t["type"]
            sm = int(t["strand"] == strand and abs((o5 + olen) - t3) <= 2)
        o.write("\t".join(map(str, [tid, oid, strand, L, o5 + 1, o5 + olen, open5, d, seg, coding, kz, us,
                                    trim, trim_run, f"{mu:.5f}", tt, t5, t3, sm])) + "\n")
print(f"{run}: mu={mu:.4f} kozak={'yes' if kmodel else 'no'} trim mismatches vs run: {n_mis}", file=sys.stderr)
