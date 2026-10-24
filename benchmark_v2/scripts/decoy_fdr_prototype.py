#!/usr/bin/env python3
"""Prototype: target-decoy FDR for candidate ORFs, with decoys from shuffled transcripts (E20).

Targets are the LongOrfs candidates of a PyTD run (longest_orfs.cds.scores). Decoys are the LongOrfs
candidates of shuffled copies of the same transcripts (same ORF finder, genetic code and strandedness),
scored here with the run's own trained hexamer model exactly as Predict scores candidates.
Statistic: frame contrast = in-frame score - max(two shifted frames); eligible ORFs have in-frame
score > 0 and contrast > 0 (as now). Within each ORF-length bin (log-spaced, ratio BIN_RATIO):
  FDR(c) = decoys with contrast >= c / targets with contrast >= c;  q(c) = min FDR over cutoffs <= c.
Bins without decoys give q = 0 (length alone suffices).
Strata (env STRATA, default 'len'): 'len' = ORF-length bin only; 'len_tx' = also transcript-length bin
(log2, from 500 nt), since a short transcript caps the ORF length; 'len_tx_type' = also complete vs
partial ORF (open at either transcript end). Decoys keep each transcript's length, so they carry the
same caps. env NSHUF: number of shuffled copies per transcript in the decoy run (decoy counts are
divided by it).

Evaluation (post-filter on the run's selected ORFs, one per transcript): for each q cutoff, stop
sensitivity on coding transcripts and the fraction of negatives (no truth row) keeping a prediction.

usage: decoy_fdr_prototype.py <run_dir> <shuffled_run_dir> <truth.tsv> <genetic_code> <label> > out.tsv
"""
import bisect
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_predictions import TOL, ends, read_truth  # noqa: E402

BIN_RATIO = 1.1
CUTOFFS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def markov_score(seq, scores, stops, order=5):
    """As TransDecoderPredict._score_cds_via_markov."""
    n, total = len(seq), 0.0
    if n < order + 1:
        return 0.0
    for i in range(n):
        frame = i % 3
        if i == n - 3 and frame == 0 and seq[i:i + 3] in stops:
            break
        total += scores.get(f"{seq[i - min(i, order):i + 1]}-{frame}", 0.0)
    return total


def length_bin(codons):
    return int(math.log(codons / 50.0) / math.log(BIN_RATIO)) if codons >= 50 else -1


STRATA = os.environ.get("STRATA", "len")
NSHUF = float(os.environ.get("NSHUF", "1"))


def stratum(codons, tx_len, orf_type):
    k = (length_bin(codons),)
    if STRATA in ("len_tx", "len_tx_type"):
        k += (max(0, int(math.log2(max(tx_len, 1) / 500.0)) + 1) if tx_len >= 500 else 0,)
    if STRATA == "len_tx_type":
        k += (orf_type == "complete",)
    return k


def cds_info(wd, fasta):
    """ORF id -> (type, transcript length) from longest_orfs.cds headers and the transcripts."""
    tlen = {k: len(v) for k, v in read_fasta(fasta)}
    info = {}
    with open(os.path.join(wd, "longest_orfs.cds")) as fh:
        for line in fh:
            if line.startswith(">"):
                x = line[1:].split()
                info[x[0]] = (x[1].split(":")[1], tlen.get(x[2].rsplit(":", 1)[0], 0))
    return info


def main():
    run, shuf, truth_file, gcode, label = sys.argv[1:6]
    wd = os.path.join(run, "input.fa.transdecoder_dir")
    sys.path.insert(0, os.environ.get("PYTD", "/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/trees/PY_taper"))
    from pytransdecoder.core.translator import Translator
    stops = frozenset(Translator(gcode).get_stop_codons())
    hex_scores = {}
    with open(os.path.join(wd, "hexamer.scores")) as fh:
        for line in fh:
            if not line.startswith("#"):
                p = line.rstrip("\n").split("\t")
                hex_scores[p[0]] = float(p[3])
    tinfo = cds_info(wd, os.path.join(run, "input.fa"))
    dinfo = cds_info(os.path.join(shuf, "input.fa.transdecoder_dir"), os.path.join(shuf, "input.fa"))
    targets = {}  # orf id -> (stratum, contrast) for eligible ORFs
    tbins, dbins = {}, {}
    with open(os.path.join(wd, "longest_orfs.cds.scores")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            s = [float(x) for x in f[3:6]]
            c = s[0] - max(s[1], s[2])
            if s[0] > 0 and c > 0:
                b = stratum(int(f[2]) // 3, tinfo[f[0]][1], tinfo[f[0]][0])
                targets[f[0]] = (b, c)
                tbins.setdefault(b, []).append(c)
    n_decoy = 0
    for did, seq in read_fasta(os.path.join(shuf, "input.fa.transdecoder_dir", "longest_orfs.cds")):
        seq = seq.upper()
        s0 = markov_score(seq, hex_scores, stops)
        c = s0 - max(markov_score(seq[1:], hex_scores, stops), markov_score(seq[2:], hex_scores, stops))
        n_decoy += 1
        if s0 > 0 and c > 0:
            dbins.setdefault(stratum(len(seq) // 3, dinfo[did][1], dinfo[did][0]), []).append(c)
    # q-values per bin
    qmap = {}
    for b, tc in tbins.items():
        tc = sorted(tc)
        dc = sorted(dbins.get(b, []))
        fdr = [(len(dc) - bisect.bisect_left(dc, c)) / NSHUF / (len(tc) - i) for i, c in enumerate(tc)]
        # q(c) = smallest FDR of any cutoff that accepts c, i.e. min over cutoffs <= c
        q, run_min = [0.0] * len(tc), math.inf
        for i in range(len(tc)):
            run_min = min(run_min, fdr[i])
            q[i] = run_min
        qmap[b] = (tc, q)

    def qvalue(orf):
        if orf not in targets:
            return 1.0
        b, c = targets[orf]
        tc, q = qmap[b]
        return q[bisect.bisect_left(tc, c)]

    truth = read_truth(truth_file)
    sel = []
    with open(os.path.join(run, "input.fa.transdecoder.gff3")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "CDS":
                sel.append((f[0], int(f[3]), int(f[4]), f[6], re.search(r"Parent=([^;]+)", f[8]).group(1)))
    n_neg = sum(1 for k, _ in read_fasta(os.path.join(run, "input.fa")) if k not in truth)
    rows = [(tid, qvalue(orf), bool(truth.get(tid)) and st == truth[tid][2]
             and abs(ends(s, e, st)[1] - ends(*truth[tid][:3])[1]) <= TOL, tid not in truth)
            for tid, s, e, st, orf in sel]
    print(f"# {label}: strata {STRATA}, {len(targets)} eligible targets, {n_decoy} decoys", file=sys.stderr)
    print("label\tq_cutoff\tsens_stop\tneg_rate\tn_pred")
    for a in CUTOFFS:
        keep = [r for r in rows if r[1] <= a]
        print(f"{label}\t{a:g}\t{sum(r[2] for r in keep) / len(truth):.4f}\t"
              f"{sum(r[3] for r in keep) / n_neg if n_neg else 0:.4f}\t{len(keep)}")
    print(f"{label}\tbase\t{sum(r[2] for r in rows) / len(truth):.4f}\t"
          f"{sum(r[3] for r in rows) / n_neg if n_neg else 0:.4f}\t{len(rows)}")


if __name__ == "__main__":
    main()
