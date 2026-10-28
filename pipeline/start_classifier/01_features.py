#!/usr/bin/env python3
"""Labelled start-decision candidates of one finished PyTransDecoder run (evidence E33, step 1).

usage: PYTHONPATH=<PyTD tree> python 01_features.py <run_dir> <genetic_code> <truth.tsv> <dataset> <out.tsv.gz>

For every ORF selected before start refinement (workdir best_candidates.gff3; refinement changes only
start codons, so this is the --no_refine_starts selection) that is open at its 5' end and has an
in-frame ATG, the script recomputes with PyTD's own code (imported from the tree on PYTHONPATH) the
quantities of the release rule (trim when alpha*|u| + w*K > C(u)/mu), reads the run's own decision
(best_candidates.gff3.revised_starts.gff3) and checks the two agree.

Only candidates whose ORF 3' end matches a truth CDS (same strand, +-2 nt, as score_predictions.py)
are labelled:
  label = truncated   truth type 5prime_partial / internal (true CDS begins before the transcript)
  label = utr         truth type complete / 3prime_partial (true start is an ATG inside the ORF)
  true_atg_rank       for utr: 1 if the true start is the first in-frame ATG, 2.. if a later one
                      (0 = truth start upstream of the first ATG, i.e. inconsistent; reported)
Unmatched candidates are written with label 'nostop' (kept for the input-level priors).

Decoy null (column dec_coding): for every candidate, the transcript is dinucleotide-shuffled exactly
as PyTD's decoys (core/decoys.py, Xoshiro128ss seed '42:<tid>:0'); in the ORF orientation a random
stop-free in-frame stretch of length |u| is taken from the shuffle (random frame, random start among
all valid positions; seeded) and scored as C/mu like u. Empty if the shuffle has no such stretch.
"""
import gzip
import random
import re
import statistics
import sys
from pathlib import Path

from pytransdecoder.predict import TransDecoderPredict, DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT
from pytransdecoder.core import kozak as kozak_mod
from pytransdecoder.core import decoys
from pytransdecoder.core.prng import Xoshiro128ss
from pytransdecoder.core.sequence import reverse_complement

run, gcode, truth_file, dataset, out = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
ALPHA, W = DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT

fa = run / "input.fa"
tdp = TransDecoderPredict(transcripts_file=fa, output_dir=run, genetic_code=gcode)
wd = tdp.workdir
hexfile = wd / "hexamer.scores"
train = wd / "__checkpoints_TDpredict" / "top_training_orfs.cds"
scores = tdp._load_hexamer_scores(hexfile)
mu = tdp._coding_contrast(train, hexfile)
kfile = wd / "start_refinement.kozak.tsv"
kmodel = kozak_mod.read_model(kfile) if kfile.exists() else None
tx = tdp._load_transcripts()
stops = set(tdp.stop_codons)

revised = {}
for line in open(wd / "best_candidates.gff3.revised_starts.gff3"):
    p = line.rstrip("\n").split("\t")
    if len(p) == 9 and p[2] == "mRNA":
        revised[re.search(r"ID=([^;]+)", p[8]).group(1)] = "start_revised=true" in p[8]

truth = {}
with open(truth_file) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(hdr, line.rstrip("\n").split("\t")))
        truth[r["transcript_id"]] = r

orfs = []
for line in open(wd / "best_candidates.gff3"):
    p = line.rstrip("\n").split("\t")
    if line.startswith("#") or len(p) != 9 or p[2] != "mRNA":
        continue
    orfs.append((p[0], int(p[3]), int(p[4]), p[6], re.search(r"ID=([^;]+)", p[8]).group(1)))
med_orf_len = statistics.median(e - s + 1 for _, s, e, _, _ in orfs)


def decoy_segment(tid, oriented_strand, seglen, rng):
    """Random stop-free in-frame stretch of length seglen from PyTD's first decoy shuffle of tid."""
    shuf = decoys.dinucleotide_shuffle(tx[tid], Xoshiro128ss(f"42:{tid}:0"))
    if oriented_strand == "-":
        shuf = reverse_complement(shuf)
    cands = []
    for f in range(3):
        run_start = f
        for i in range(f, len(shuf) - 2, 3):
            if shuf[i:i + 3] in stops:
                if i - run_start >= seglen:
                    cands.append((run_start, i))
                run_start = i + 3
        last = f + 3 * ((len(shuf) - f) // 3)
        if last - run_start >= seglen:
            cands.append((run_start, last))
    if not cands:
        return None
    # pick a valid start uniformly over all positions
    weights = [(e - s - seglen) // 3 + 1 for s, e in cands]
    k = rng.choices(range(len(cands)), weights=weights)[0]
    s, e = cands[k]
    start = s + 3 * rng.randrange(weights[k])
    return shuf[start:start + seglen]


rng = random.Random(f"E33:{dataset}")
cols = ["dataset", "transcript_id", "orf_id", "strand", "tx_len", "orf5", "orf_len", "seg", "coding",
        "coding_pb", "kozak", "kozak_off", "utr_score", "trim", "trim_run", "internal", "utr3",
        "orf_len_rel", "seg_gc", "mu", "truth_type", "label", "true_atg_rank", "true_start_rel",
        "dec_coding"]
n_mis = n_open = n_noatg = 0
with gzip.open(out, "wt") as o:
    o.write("\t".join(cols) + "\n")
    for tid, s, e, strand, oid in orfs:
        seq = tx[tid]
        L = len(seq)
        if strand == "-":
            oriented, o5 = reverse_complement(seq), L - e
        else:
            oriented, o5 = seq, s - 1
        olen = e - s + 1
        orf = oriented[o5:o5 + olen]
        if orf[:3] == "ATG":
            continue
        n_open += 1
        atg = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"), None)
        if atg is None:
            n_noatg += 1
            continue
        u = orf[:atg]
        c = tdp._frame_contrast(u, scores) / mu
        k = kozak_mod.score(kmodel, oriented, o5 + atg) if kmodel else 0.0
        koff = sum(1 for b in kozak_mod.window(oriented, o5 + atg) if b is None)
        us = ALPHA * atg + W * k
        trim = int(us > c)
        trim_run = int(revised[oid])
        n_mis += trim != trim_run
        internal = int(orf[-3:] not in stops)
        utr3 = L - (o5 + olen)
        gc = (u.count("G") + u.count("C")) / len(u)
        t = truth.get(tid)
        tt = lab = rank = rel = ""
        if t and t["strand"] == strand:
            ts, te = int(t["cds_start"]), int(t["cds_end"])
            t5, t3 = (L - te + 1, L - ts + 1) if strand == "-" else (ts, te)
            if abs((o5 + olen) - t3) <= 2:
                tt = t["type"]
                rel = (t5 - 1) - o5  # true start relative to the ORF's open 5' end (nt)
                if tt in ("5prime_partial", "internal"):
                    lab = "truncated"
                else:
                    lab = "utr"
                    atgs = [i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"]
                    near = [j + 1 for j, i in enumerate(atgs) if abs(i - rel) <= 2]
                    rank = near[0] if near else 0
        if not lab:
            lab = "nostop"
        dseg = decoy_segment(tid, strand, atg, rng) if lab != "nostop" else None
        dc = "" if dseg is None else f"{tdp._frame_contrast(dseg, scores) / mu:.3f}"
        o.write("\t".join(map(str, [dataset, tid, oid, strand, L, o5, olen, atg, f"{c:.3f}", f"{c / atg:.4f}",
                                    f"{k:.3f}", koff, f"{us:.3f}", trim, trim_run, internal, utr3,
                                    f"{olen / med_orf_len:.4f}", f"{gc:.4f}", f"{mu:.5f}", tt, lab, rank,
                                    rel, dc])) + "\n")
print(f"{dataset}\t{run}\tmu={mu:.4f}\tkozak={'yes' if kmodel else 'no'}\topen5={n_open}\tno_atg={n_noatg}"
      f"\ttrim_mismatch_vs_run={n_mis}", file=sys.stderr)
