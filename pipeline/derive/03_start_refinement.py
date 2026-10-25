#!/usr/bin/env python3
"""Compact tables for notebook 03 (start-codon refinement; evidence E09, E13, E25).

usage: PYTHONPATH=$P/work/trees/PY_v7 $P/env/td/bin/python 03_start_refinement.py [P]
       (P = /broad/hptmp/bhaas/TRANSDECODER_PAPER; PY_v7 = PyTransDecoder 3a793ff, for core/kozak.py)

Writes data/start_refinement/:
  kozak_models.tsv    Kozak (start-context) counts and log2-odds, 8 organisms x 8 positions x 4 bases,
                      for two models trained on the stranded 'complete' input:
                        self   = the model the release run trained (pytd_v7,
                                 input.fa.transdecoder_dir/start_refinement.kozak.counts.tsv)
                        truth  = annotated CDS starts of transcripts with an annotated 5' UTR vs
                                 in-frame ATGs inside those CDSs (same training code, kozak.py),
                                 as in figures/kozak/make_kozak_figure.py
  training_mu.tsv     per organism and input (complete, fragments), from the pytd_v7 predict.log:
                      the median per-base frame contrast of the training ORFs (mu) and the number
                      of Kozak training starts
  later_atg.tsv       annotated starts (complete input, + strand, annotated 5' UTR, ATG start,
                      20 nt upstream and 10 downstream in the transcript) classified as in
                      pwm_audit/audit_signal.py Q4: closed/open (an in-frame stop upstream within
                      the transcript or not) x first/later (an in-frame ATG upstream of the
                      annotated start with no stop in between)
  seed_start_changes.tsv  v6.0.0-era PWM seeds 2-5 vs seed 1 on the PASA and StringTie sample data
                      ($P/work/start_seed_runs; E09): ORFs, same stop with a different start, and
                      ORFs whose stop is not shared
"""
import re
import sys
from pathlib import Path

from pytransdecoder.core import kozak

P = Path(sys.argv[1] if len(sys.argv) > 1 else "/broad/hptmp/bhaas/TRANSDECODER_PAPER")
INPUTS = P / "data/bench_v3/inputs"
RUNS = P / "work/bench_v3/runs"
SEEDS = P / "work/start_seed_runs"
OUT = P / "TransDecoderPaper/data/start_refinement"
ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
COMP = str.maketrans("ACGT", "TGCA")
STOPS = {1: {"TAA", "TAG", "TGA"}, 6: {"TGA"}}
LABELS = kozak.position_labels()


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def truth_rows(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(hdr, line.rstrip("\n").split("\t")))


def read_counts(path):
    pos, neg = kozak.new_counts(), kozak.new_counts()
    for line in open(path):
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        k = LABELS.index(x[0])
        for i, b in enumerate(kozak.BASES):
            pos[k][b] = int(x[1 + i])
            neg[k][b] = int(x[5 + i])
    return pos, neg


def truth_counts(seqs, truth):
    """Same as make_kozak_figure.py train() on truth_sites(): annotated 5' UTR, ATG starts only."""
    pos, neg, n = kozak.new_counts(), kozak.new_counts(), 0
    for r in truth:
        if r["utr5_source"] != "annotated":
            continue
        seq = seqs[r["transcript_id"]]
        start, end = int(r["cds_start"]) - 1, int(r["cds_end"])
        if r["strand"] == "-":
            seq = seq.translate(COMP)[::-1]
            start, end = len(seq) - end, len(seq) - start
        if seq[start:start + 3] != "ATG":
            continue
        n += 1
        kozak.add_site(pos, seq, start)
        for i in range(start + 3, end - 3, 3):
            if seq[i:i + 3] == "ATG":
                kozak.add_site(neg, seq, i)
    return pos, neg, n


def later_atg(seqs, truth):
    """audit_signal.py Q4 (pwmlib.truth_sets): + strand, annotated 5' UTR, ATG, full 31-nt window."""
    c = dict(n=0, closed_first=0, closed_later=0, open_first=0, open_later=0)
    for r in truth:
        if r["utr5_source"] != "annotated" or r["strand"] != "+":
            continue
        s = seqs[r["transcript_id"]]
        cs = int(r["cds_start"]) - 1
        if s[cs:cs + 3] != "ATG":
            continue
        b, e = cs - 20, cs + 13
        if b < 0 or e > len(s) or any(ch not in "GATC" for ch in s[b:e]):
            continue
        stops = STOPS[int(r["transl_table"])]
        i, ups = cs - 3, 0
        while i >= 0:
            codon = s[i:i + 3]
            if codon in stops:
                break
            if codon == "ATG":
                ups += 1
            i -= 3
        key = ("closed" if i >= 0 else "open") + ("_first" if ups == 0 else "_later")
        c["n"] += 1
        c[key] += 1
    return c


def gff_cds(path):
    spans = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        parent = re.search(r"Parent=([^;]+)", f[8]).group(1)
        lo, hi = int(f[3]), int(f[4])
        s = spans.setdefault(parent, [f[0], f[6], lo, hi])
        s[2], s[3] = min(s[2], lo), max(s[3], hi)
    return {tuple(v) for v in spans.values()}


def stop_key(o):
    return (o[0], o[1], o[3] if o[1] == "+" else o[2])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mu_pat = re.compile(r"median frame contrast ([0-9.]+) per base")
    kz_pat = re.compile(r"Kozak model from (\d+) training starts")
    with open(OUT / "kozak_models.tsv", "w") as km, open(OUT / "training_mu.tsv", "w") as tm, \
            open(OUT / "later_atg.tsv", "w") as la:
        km.write("org\tmodel\tn_starts\tposition\tbase\tstart_count\tinternal_count\tlog2odds\n")
        tm.write("org\tinput\tmu_per_base\tkozak_training_starts\n")
        la.write("org\tn\tclosed_first\tclosed_later\topen_first\topen_later\n")
        for org in ORGS:
            for inp in ("complete", "fragments"):
                log = (RUNS / org / inp / "pytd_v7" / "predict.log").read_text()
                tm.write(f"{org}\t{inp}\t{mu_pat.search(log).group(1)}\t{kz_pat.search(log).group(1)}\n")
            seqs = fasta(INPUTS / org / "complete.fa")
            truth = list(truth_rows(INPUTS / org / "complete.truth.tsv"))
            d = RUNS / org / "complete" / "pytd_v7" / "input.fa.transdecoder_dir"
            pos, neg = read_counts(d / "start_refinement.kozak.counts.tsv")
            n_self = int(re.search(r"(\d+) training starts", (d / "start_refinement.kozak.tsv").read_text()).group(1))
            tpos, tneg, n_truth = truth_counts(seqs, truth)
            for model, (p, ng, n) in (("self", (pos, neg, n_self)), ("truth", (tpos, tneg, n_truth))):
                lo = kozak.log_odds_model(p, ng)
                for k, lab in enumerate(LABELS):
                    for b in kozak.BASES:
                        km.write(f"{org}\t{model}\t{n}\t{lab}\t{b}\t{p[k][b]}\t{ng[k][b]}\t{lo[k][b]:.4f}\n")
            c = later_atg(seqs, truth)
            la.write("\t".join([org] + [str(c[k]) for k in
                                        ("n", "closed_first", "closed_later", "open_first", "open_later")]) + "\n")
            print(org, "done", file=sys.stderr)
    with open(OUT / "seed_start_changes.tsv", "w") as sc:
        sc.write("dataset\tseed\tn_orfs\tn_orfs_seed1\tsame_stop_diff_start\tstop_not_shared\n")
        for ds in ("pasa", "stringtie"):
            a = gff_cds(SEEDS / ds / "seed1" / f"{ds}.fa.transdecoder.gff3")
            for seed in range(2, 6):
                b = gff_cds(SEEDS / ds / f"seed{seed}" / f"{ds}.fa.transdecoder.gff3")
                ao = {stop_key(o) for o in a - b}
                bo = {stop_key(o) for o in b - a}
                diff = len(ao & bo)
                sc.write(f"{ds}\t{seed}\t{len(b)}\t{len(a)}\t{diff}\t{len(ao ^ bo)}\n")


if __name__ == "__main__":
    main()
