#!/usr/bin/env python3
"""E35: start-decision quantities of every 5'-open candidate of one finished PyTD run (read only).

usage: PYTHONPATH=<PyTD a5b24c3 tree> python 01_candidates_v1.py <run_dir> <genetic_code> <dataset> <out.tsv.gz>

For every ORF selected before start refinement (workdir best_candidates.gff3) that does not begin with
ATG and has a first in-frame ATG, recompute with PyTD's own code exactly as _revise_start_sites:
  seg     = |u|, nt between the open 5' end and the first in-frame ATG
  coding  = frame contrast of u / training median per-base contrast (mu)
  kozak   = Kozak score K (bits) of the ATG with the run's model (start_refinement.kozak.tsv; 0 if none)
  trim_ab = release ab initio rule: 0.4*seg + 17*K > coding
Homology runs (start_refinement.homology.tsv present): the run's own BLAST qstart/sstart and Pfam
env_from per ORF are passed to PyTD's homology_start_call -> hcall (trim / keep / conflict / none);
final release decision = hcall if trim/keep, else trim_ab.
trim_run = the run's own decision (start_revised=true in best_candidates.gff3.revised_starts.gff3).
Floats are written with repr() (full precision) so variant decisions are exact.
Stderr: one summary line incl. the number of recomputed-vs-run mismatches (must be 0).
"""
import gzip
import re
import sys
from pathlib import Path

from pytransdecoder.predict import TransDecoderPredict, DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT
from pytransdecoder.core import kozak as kozak_mod
from pytransdecoder.core.homology_starts import homology_start_call
from pytransdecoder.core.sequence import reverse_complement

run, gcode, dataset, out = Path(sys.argv[1]).resolve(), sys.argv[2], sys.argv[3], sys.argv[4]
ALPHA, W = DEFAULT_START_CONTRAST_ALPHA, DEFAULT_START_KOZAK_WEIGHT
assert (ALPHA, W) == (0.4, 17), (ALPHA, W)

tdp = TransDecoderPredict(transcripts_file=run / "input.fa", output_dir=run, genetic_code=gcode)
wd = tdp.workdir
hexfile = wd / "hexamer.scores"
scores = tdp._load_hexamer_scores(hexfile)
mu = tdp._coding_contrast(wd / "__checkpoints_TDpredict" / "top_training_orfs.cds", hexfile)
kfile = wd / "start_refinement.kozak.tsv"
kmodel = kozak_mod.read_model(kfile) if kfile.exists() and not (wd / "start_refinement.kozak.SKIPPED").exists() else None
tx = tdp._load_transcripts()

revised = {}
for line in open(wd / "best_candidates.gff3.revised_starts.gff3"):
    p = line.rstrip("\n").split("\t")
    if len(p) == 9 and p[2] == "mRNA":
        revised[re.search(r"ID=([^;]+)", p[8]).group(1)] = "start_revised=true" in p[8]

htab = None
hfile = wd / "start_refinement.homology.tsv"
if hfile.exists():
    htab = {}
    with open(hfile) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            htab[r["orf_id"]] = r

cols = ["dataset", "transcript_id", "orf_id", "strand", "tx_len", "orf5", "orf_len", "seg", "coding", "kozak",
        "trim_ab", "hcall", "trim_rel", "trim_run"]
n = n_mis = n_hmis = 0
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
            oriented, o5 = reverse_complement(seq), L - e
        else:
            oriented, o5 = seq, s - 1
        olen = e - s + 1
        orf = oriented[o5:o5 + olen]
        if orf[:3] == "ATG":
            continue
        atg = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"), None)
        if atg is None:
            continue
        n += 1
        c = tdp._frame_contrast(orf[:atg], scores) / mu
        k = kozak_mod.score(kmodel, oriented, o5 + atg) if kmodel else 0.0
        trim_ab = c < ALPHA * atg + W * k
        hcall = "none"
        trim_rel = trim_ab
        if htab is not None:
            h = htab[oid]
            if int(h["atg_offset_nt"]) != atg:
                n_hmis += 1
            blast = None if h["blast_qstart"] == "." else (int(h["blast_qstart"]), int(h["blast_sstart"]))
            env = None if h["pfam_env_from"] == "." else int(h["pfam_env_from"])
            _, call = homology_start_call(atg // 3 + 1, blast, env)
            hcall = call or "none"
            if call in ("trim", "keep"):
                trim_rel = call == "trim"
            if (h["decision"] == "trim") != trim_rel or (h["ab_initio_decision"] == "trim") != trim_ab:
                n_hmis += 1
        n_mis += trim_rel != revised[oid]
        o.write("\t".join(map(str, [dataset, tid, oid, strand, L, o5, olen, atg, repr(c), repr(k), int(trim_ab), hcall,
                                    int(trim_rel), int(revised[oid])])) + "\n")
print(f"{dataset}\t{run}\tmu={mu!r}\tkozak={'yes' if kmodel else 'no'}\thomology={'yes' if htab is not None else 'no'}"
      f"\tcandidates={n}\tmismatch_vs_run={n_mis}\tmismatch_vs_homology_table={n_hmis}", file=sys.stderr)
