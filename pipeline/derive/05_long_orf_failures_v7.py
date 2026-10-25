#!/usr/bin/env python3
"""Coding-test failures of true CDSs under the RELEASE rule (PyTD 3a793ff, `pytd_v7` runs).

Copy of analyses/long_orf_failures/long_orf_failures.py (E12), changed for the release (E25):

  rule 3 (coding)  S1 > 0, C = S1 - max(S2, S3) > 0, and target-decoy q <= 0.05
                   (q from the run's input.fa.transdecoder_dir/decoy_qvalues.tsv; replaces the fixed
                   short-ORF test C >= 200 * mu below 300 nt of the development rule)
  rule 2 (length)  ORF length >= L*(GC) (off for TGA-only codes, e.g. T. thermophila)

  fail_reason      S1<=0 | shifted_frame>=S1 | decoy_q>0.05
  runs             work/bench_v3/runs/<org>/<input>/pytd_v7 (LongOrfs at 50 aa, as pytd_dev_v2)
  homology flags   pytd_v7_homology where finished, else pytd_dev_v2_homology. Both search the same
                   LongOrfs ORFs (>= 50 aa) against the same Swiss-Prot (organism removed) and Pfam;
                   their hits were checked identical for athal, dmel and human.
  truth-trained    alternative Markov models (two folds by gene; all / rand500 / long500 / GC3-matched
  models           500) now pass an ORF only with S1 > 0, C > 0 and q <= 0.05, where q is computed as
                   PyTD does, from decoys scored with that model: the decoy ORFs of the run are
                   regenerated exactly as PyTD (dinucleotide shuffles, seed 42, same ORF finder) and
                   every LongOrfs ORF and decoy ORF is rescored with each model. With the run's own
                   model this reproduces decoy_qvalues.tsv (agreement written to the meta table).
  altrules         A (S1 > 0 dropped) recomputes q with C > 0 as the only eligibility condition.
  GM12878          not part of this copy (not used in the coding-failure section).

usage:
  05_long_orf_failures_v7.py bench [--orgs ...] [--input complete|complete_fullutr] [--jobs 2]
  05_long_orf_failures_v7.py altrules
  05_long_orf_failures_v7.py summarize
  05_long_orf_failures_v7.py all [--jobs 2]
Outputs: $P/work/long_orf_failures_v7/ (override with --out). Run with $P/env/td/bin/python.
"""
import argparse
import gzip
import os
import random
import re
import sys
import zlib
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOF = HERE.parents[1] / "analyses/long_orf_failures"  # markov.py, features.py (unchanged)
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(LOF))
P = Path(os.environ.get("TD_PAPER_ROOT", "/broad/hptmp/bhaas/TRANSDECODER_PAPER"))
T = P / "TransDecoderPaper"
PYTD = Path(os.environ.get("PYTD_TREE", P / "work/trees/PY_v7"))
sys.path.insert(0, str(PYTD))
sys.path.insert(0, str(T / "genome_annotation/scripts"))

from markov import MarkovModel, median_contrast, read_bg  # noqa: E402
import features as F  # noqa: E402

ORGS = ["athal", "creinhardtii", "dmel", "human", "mouse", "pfalciparum", "spombe", "tthermophila"]
RUN = "pytd_v7"
DECOY_Q = 0.05          # release default (predict.py DEFAULT_DECOY_Q)
DECOY_SEED = 42
DECOY_TARGET_NT = 20_000_000
DECOY_MAX_SHUFFLES = 50
LEN_BINS = [(0, 300), (300, 600), (600, 1000), (1000, 2000), (2000, 5000), (5000, 10 ** 9)]
GC_TABLE = [(25, 519), (30, 552), (35, 573), (40, 615), (45, 663), (50, 726), (55, 804), (60, 927),
            (65, 1083), (70, 1347), (75, 1746), (80, 2436)]
TRAIN_N = 500


def len_bin(n):
    for lo, hi in LEN_BINS:
        if lo <= n < hi:
            return f"{lo}-{hi - 1}" if hi < 10 ** 9 else f">={lo}"
    return "NA"


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf).upper()


def read_fasta_headers(path):
    out = {}
    for line in open(path):
        if line.startswith(">"):
            p = line[1:].split()
            m = re.search(r"type:(\S+)", line)
            out[p[0]] = m.group(1) if m else "NA"
    return out


def read_orf_gff3(path):
    """transcript -> list of (orf_id, start, end, strand) from the mRNA lines of longest_orfs.gff3."""
    out = defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] == "mRNA":
            oid = re.search(r"ID=([^;]+)", f[8]).group(1)
            out[f[0]].append((oid, int(f[3]), int(f[4]), f[6]))
    return out


def read_scores(path):
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        out[p[0]] = [float(x) for x in p[3:9]]
    return out


def parse_predict_log(path):
    """L* (inf when off), the training ORFs' median per-base frame contrast (start refinement line),
    transcriptome GC %, and the number of decoy shuffles per transcript."""
    txt = open(path).read()
    m = re.search(r"min auto-accept length = (\S+) nt", txt)
    lstar = float(m.group(1)) if m else float("inf")
    if "no length-based auto-accept" in txt:
        lstar = float("inf")
    m = re.search(r"median frame contrast ([0-9.]+) per base", txt)
    mu = float(m.group(1)) if m else None
    m = re.search(r"Transcriptome GC content: ([0-9.]+)%", txt)
    gc = float(m.group(1)) if m else float("nan")
    m = re.search(r"Decoys: (\d+) shuffle", txt)
    nshuf = int(m.group(1)) if m else None
    return lstar, mu, gc, nshuf


def read_qvalues(path):
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        out[p[0]] = float(p[3])
    return out


def read_preds_tsv(path):
    out = defaultdict(list)
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if p[0] == "transcript_id":
            continue
        out[p[0]].append((int(p[1]), int(p[2]), p[3]))
    return out


def stop_end(start, end, strand):
    return end if strand == "+" else start


def match_stop(orfs, s, e, strand, tol=2):
    """ORF (from (id, start, end, strand) list) with the same strand and stop within tol nt."""
    target = stop_end(s, e, strand)
    best = None
    for o in orfs:
        if o[3] != strand:
            continue
        d = abs(stop_end(o[1], o[2], o[3]) - target)
        if d <= tol and (best is None or d < best[0]):
            best = (d, o)
    return best[1] if best else None


def passes_rule3(sc, q):
    """Release coding test: S1 > 0, frame contrast > 0 and target-decoy q <= DECOY_Q.
    q is None for ORFs without a q-value (not eligible), which then fail."""
    s1, c = sc[0], sc[0] - max(sc[1], sc[2])
    q_ok = q is not None and q <= DECOY_Q
    return s1 > 0 and c > 0 and q_ok, s1 > 0, c > 0, q_ok


# ------------------------------------------------------------------ decoys and batch rescoring

def load_pytd_decoys():
    from pytransdecoder.core import decoys
    from pytransdecoder.core.prng import Xoshiro128ss
    from pytransdecoder.core.orf_finder import ORFFinder
    return decoys, Xoshiro128ss, ORFFinder


def decoy_orfs(run, cache, n_shuf_log=None):
    """Decoy ORFs of a run, regenerated as predict.py _decoy_qvalues does (same shuffles, seed and
    ORF finder settings). Returns (n_shuf, [(codons, tx_len, complete, cds_seq)]); cached as a
    gzipped FASTA in `cache`."""
    import json
    decoys, Xoshiro128ss, ORFFinder = load_pytd_decoys()
    wd = run / "input.fa.transdecoder_dir"
    transcripts = {k: v for k, v in read_fasta(run / "input.fa")}
    total_nt = sum(len(x) for x in transcripts.values())
    n_shuf = max(1, min(DECOY_MAX_SHUFFLES, -(-DECOY_TARGET_NT // max(total_nt, 1))))
    assert n_shuf_log is None or n_shuf == n_shuf_log, (n_shuf, n_shuf_log)
    out = []
    if cache.exists():
        with gzip.open(cache, "rt") as fh:
            hdr = None
            for line in fh:
                if line.startswith(">"):
                    hdr = line[1:].split()
                else:
                    out.append((int(hdr[1]), int(hdr[2]), hdr[3] == "1", line.strip()))
        return n_shuf, out
    params = json.loads((wd / "longorfs.params.json").read_text())
    finder = ORFFinder(min_protein_length=params["min_protein_length"],
                       allow_5prime_partial=not params["complete_orfs_only"],
                       allow_3prime_partial=not params["complete_orfs_only"],
                       forward_strand=True, reverse_strand=not params["strand_specific"],
                       genetic_code=params["genetic_code"],
                       complete_orfs_only=params["complete_orfs_only"])
    for tid in sorted(transcripts):
        seq = transcripts[tid]
        for k in range(n_shuf):
            shuf = decoys.dinucleotide_shuffle(seq, Xoshiro128ss(f"{DECOY_SEED}:{tid}:{k}"))
            for orf in finder.find_all_orfs(shuf, f"decoy{k}_{tid}"):
                cds = orf.sequence.upper()
                out.append((len(cds) // 3, len(seq), orf.orf_type == "complete", cds))
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(".tmp")
    with gzip.open(tmp, "wt") as fh:
        for i, (cod, tl, comp, cds) in enumerate(out):
            fh.write(f">d{i} {cod} {tl} {int(comp)}\n{cds}\n")
    tmp.rename(cache)
    return n_shuf, out


def target_orfs(wd, seqs):
    """LongOrfs ORFs: id -> (codons, tx_len, complete, cds_seq), as predict.py builds its targets."""
    out = {}
    for line_id, cds in read_fasta(wd / "longest_orfs.cds"):
        out[line_id] = [len(cds) // 3, None, None, cds]
    for line in open(wd / "longest_orfs.cds"):
        if line.startswith(">"):
            x = line[1:].split()
            tid = x[2].rsplit(":", 1)[0]
            out[x[0]][1] = len(seqs.get(tid, ""))
            out[x[0]][2] = x[1].split(":", 1)[1] == "complete"
    return {k: tuple(v) for k, v in out.items()}


def _tables(model):
    """Flattened 6-mer table (3*4096) and the short k-mer tables of a markov.MarkovModel."""
    return model.tables[6].reshape(-1), {k: model.tables[k] for k in range(1, 6)}


def batch_frame_scores(seqs, models, stops, chunk_nt=4_000_000):
    """Forward frame scores S1, S2, S3 (offsets 0, 1, 2) of every sequence under every model, as
    markov.MarkovModel.score (= PyTD _score_cds_via_markov), rounded to 2 decimals as PyTD writes
    them. Vectorised over sequences: one pass builds the 6-mer indices of a chunk of sequences and
    every model is applied to it. Returns {model_name: array (n, 3)}."""
    import numpy as np
    from markov import encode
    tabs = {name: _tables(m) for name, m in models.items()}
    res = {name: np.zeros((len(seqs), 3)) for name in models}
    i0 = 0
    while i0 < len(seqs):
        i1, tot = i0, 0
        while i1 < len(seqs) and (tot < chunk_nt or i1 == i0):
            tot += len(seqs[i1])
            i1 += 1
        chunk = seqs[i0:i1]
        lens = np.array([len(x) for x in chunk], dtype=np.int64)
        starts = np.concatenate([[0], np.cumsum(lens)[:-1]])
        codes = encode("".join(chunk)).astype(np.int64)
        n = len(codes)
        sid = np.repeat(np.arange(len(chunk)), lens)
        pos = np.arange(n) - starts[sid]
        idx = np.zeros(n, dtype=np.int64)
        valid = pos >= 5
        for t in range(6):                       # k-mer codes[j-5 .. j]
            sh = 5 - t
            c = np.empty(n, dtype=np.int64)
            c[sh:] = codes[:n - sh]
            c[:sh] = 4
            valid &= c < 4
            idx = idx * 4 + np.where(c < 4, c, 0)
        term = np.array([x[-3:] in stops for x in chunk])
        masks = []
        for o in range(3):
            Lo = lens - o
            n_o = np.where((Lo >= 3) & ((Lo - 3) % 3 == 0) & term, Lo - 3, Lo)
            ok_seq = (Lo >= 6) & (n_o > 0)
            mk = valid & (pos >= o + 5) & (pos < o + n_o[sid]) & ok_seq[sid]
            fr = (pos - o) % 3
            masks.append((mk, fr[mk] * 4096 + idx[mk], sid[mk], n_o, ok_seq))
        # short k-mers (k = i+1 <= 5) at the start of each offset
        short = []
        for o in range(3):
            n_o, ok_seq = masks[o][3], masks[o][4]
            for i in range(5):
                sel = ok_seq & (i < n_o) & (o + i < lens)
                kidx = np.zeros(len(chunk), dtype=np.int64)
                kval = sel.copy()
                for t in range(i + 1):
                    p = np.minimum(starts + o + t, n - 1)
                    c = codes[p]
                    kval &= c < 4
                    kidx = kidx * 4 + np.where(c < 4, c, 0)
                short.append((o, i, kval, kidx))
        for name in models:
            t6, tk = tabs[name]
            out = res[name]
            for o in range(3):
                mk, flat, sids, _, _ = masks[o]
                out[i0:i1, o] = np.bincount(sids, weights=t6[flat], minlength=len(chunk))
            for o, i, kval, kidx in short:
                add = np.where(kval, tk[i + 1][i % 3][kidx], 0.0)
                out[i0:i1, o] += add
        i0 = i1
    # Python round(): the correctly rounded value, = float(f"{x:.2f}") as written to
    # longest_orfs.cds.scores and as predict.py rounds decoy scores (np.round can differ by 1 ulp)
    return {k: np.array([round(float(x), 2) for x in v.ravel()]).reshape(v.shape) for k, v in res.items()}


def model_qvalues(t_meta, t_scores, d_meta, d_scores, n_shuf, require_s1=True):
    """q-values (PyTD core/decoys.py) of target ORFs from their and the decoys' frame scores.
    t_meta/d_meta: lists of (codons, tx_len, complete); scores: arrays (n, 3).
    Eligible: S1 > 0 (if require_s1) and contrast > 0. Returns array of q (nan if not eligible)."""
    import numpy as np
    decoys, _, _ = load_pytd_decoys()

    def elig(sc):
        c = sc[:, 0] - np.maximum(sc[:, 1], sc[:, 2])   # not rounded, as predict.py
        e = c > 0
        if require_s1:
            e &= sc[:, 0] > 0
        return e, c
    te, tc = elig(t_scores)
    de, dc = elig(d_scores)
    targets = {i: (decoys.stratum(m[0], m[1], m[2], "len_tx_type"), float(tc[i]))
               for i, m in enumerate(t_meta) if te[i]}
    dl = [(decoys.stratum(m[0], m[1], m[2], "len_tx_type"), float(dc[i]))
          for i, m in enumerate(d_meta) if de[i]]
    q = decoys.qvalues(targets, dl, n_shuf)
    out = np.full(len(t_meta), np.nan)
    for i, v in q.items():
        out[i] = float(f"{v:.5g}")
    return out


def orgs_table():
    rows = {}
    lines = open(T / "benchmark_v2/organisms.tsv").read().splitlines()
    hdr = lines[0].split("\t")
    for line in lines[1:]:
        r = dict(zip(hdr, line.split("\t")))
        rows[r["org"]] = r
    return rows


def annotation_path(org):
    for line in open(T / "benchmark_v2/truth_sources.tsv").read().splitlines()[1:]:
        p = line.split("\t")
        if p[0] == org:
            return P / "data/annot" / p[2]
    return None


def homology_source(org, inp="complete"):
    for t in ("pytd_v7_homology", "pytd_dev_v2_homology"):
        d = P / f"work/bench_v3/runs/{org}/{inp}/{t}"
        if (d / "preds.tsv").exists() and (d / "blastp.outfmt6").exists() and (d / "pfam.domtblout").exists():
            return t
    return None


def homology_lookup(org, tids, inp="complete"):
    """(transcript, stop) -> Pfam / Swiss-Prot hits of the LongOrfs ORF (>= 50 aa) from the finished
    homology run (pytd_v7_homology, else pytd_dev_v2_homology: same ORFs and searches); matched by
    coordinates, not IDs."""
    d = P / f"work/bench_v3/runs/{org}/{inp}/{homology_source(org, inp)}"
    orfs = read_orf_gff3(d / "input.fa.transdecoder_dir/longest_orfs.gff3")
    pfam = defaultdict(list)
    for line in open(d / "pfam.domtblout"):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) > 12 and p[4].startswith("PF"):
            if p[3] not in pfam[p[0]]:
                pfam[p[0]].append(p[3])
    sprot = {}
    for line in open(d / "blastp.outfmt6"):
        p = line.rstrip("\n").split("\t")
        if p[0] not in sprot:
            sprot[p[0]] = (p[1], float(p[2]), p[10])
    return orfs, pfam, sprot


def train_oracle_models(truth_rows, stops, bg, seed=1):
    """Markov models trained on annotated CDSs, two folds by gene: for fold f, models trained on
    genes of fold 1-f. Returns {fold: {name: (model, mu)}} and GC3 quintile breakpoints."""
    by_fold = {0: [], 1: []}
    for r in truth_rows:
        by_fold[r["fold"]].append(r)
    gc3s = sorted(r["gc3"] for r in truth_rows)
    qs = [gc3s[int(len(gc3s) * q / 5)] for q in range(1, 5)]
    models = {}
    for f in (0, 1):
        train = by_fold[1 - f]
        rng = random.Random(seed + f)
        seqs = [r["cds_seq"] for r in train]
        sets = {
            "all": seqs,
            "rand500": rng.sample(seqs, min(TRAIN_N, len(seqs))),
            "long500": [s for s in sorted((s for s in seqs if len(s) <= 5000), key=len, reverse=True)][:TRAIN_N],
        }
        for q in range(5):
            pool = [r["cds_seq"] for r in train if gc3_quintile(r["gc3"], qs) == q]
            sets[f"gc3q{q}"] = rng.sample(pool, min(TRAIN_N, len(pool)))
        models[f] = {}
        for name, ss in sets.items():
            m = MarkovModel.train(ss, bg, stops)
            mu_set = ss if len(ss) <= 2000 else random.Random(seed).sample(ss, 2000)
            models[f][name] = (m, median_contrast(m, mu_set, stops))
    return models, qs


def gc3_quintile(v, qs):
    q = 0
    while q < 4 and v >= qs[q]:
        q += 1
    return q


def process_org(org, out_dir, inp="complete", with_models=True):
    orgs = orgs_table()
    gcode = orgs[org]["td_genetic_code"]
    from pytransdecoder.core.translator import Translator
    tr = Translator(gcode)
    stops = set(tr.get_stop_codons())
    run = P / f"work/bench_v3/runs/{org}/{inp}/{RUN}"
    wd = run / "input.fa.transdecoder_dir"
    lstar, mu_log, gc_log, nshuf_log = parse_predict_log(run / "predict.log")
    qvals = read_qvalues(wd / "decoy_qvalues.tsv")
    seqs = dict(read_fasta(run / "input.fa"))
    orfs = read_orf_gff3(wd / "longest_orfs.gff3")
    orf_types = read_fasta_headers(wd / "longest_orfs.cds")
    scores = read_scores(wd / "longest_orfs.cds.scores")
    train = dict(read_fasta(wd / "__checkpoints_TDpredict/top_training_orfs.cds"))
    preds = read_preds_tsv(run / "preds.tsv")
    model = MarkovModel.from_scores_file(wd / "hexamer.scores")
    bg = read_bg(wd / "base_freqs.dat")
    mu = median_contrast(model, list(train.values()), stops)   # per-base; descriptive only
    annot = F.read_annotation(annotation_path(org))
    h_orfs, h_pfam, h_sprot = homology_lookup(org, seqs, inp)

    # training-set codon usage and GC3
    train_codons = defaultdict(int)
    train_gc3 = []
    for s in train.values():
        for c, n in F.codon_counts(s, stops).items():
            train_codons[c] += n
        train_gc3.append(F.gc3(s))

    truth_rows = []
    for line in open(P / f"data/bench_v3/inputs/{org}/{inp}.truth.tsv").read().splitlines()[1:]:
        p = line.split("\t")
        tid, cs, ce, strand, tlen, gid = p[0], int(p[1]), int(p[2]), p[3], int(p[4]), p[5]
        utr_src = (p[8], p[9]) if len(p) > 9 else ("", "")
        tseq = seqs[tid]
        cds = tseq[cs - 1:ce] if strand == "+" else F_revcomp(tseq[cs - 1:ce])
        truth_rows.append(dict(org=org, transcript_id=tid, gene_id=gid, cds_start=cs, cds_end=ce,
                               strand=strand, tx_len=tlen, cds_seq=cds, gc3=F.gc3(cds),
                               fold=zlib.crc32(gid.encode()) % 2, utr_src=utr_src))
    all_codons = defaultdict(int)
    for r in truth_rows:
        for c, n in F.codon_counts(r["cds_seq"], stops).items():
            all_codons[c] += n
    sense = sorted(set(all_codons) | set(train_codons))
    tot_t = sum(train_codons.values()) + len(sense)
    tot_a = sum(all_codons.values()) + len(sense)
    import math
    codon_llr = {c: math.log((train_codons[c] + 1) / tot_t) - math.log((all_codons[c] + 1) / tot_a)
                 for c in sense}

    oracle, gc3_qs = train_oracle_models(truth_rows, stops, bg)
    check = {}
    model_q = {}      # model key -> {orf_id: (S1, C, q)} for truth-matched ORFs
    if with_models:
        import numpy as np
        tgt = target_orfs(wd, seqs)
        tids = sorted(tgt)
        n_shuf, dec = decoy_orfs(run, out_dir / "decoys" / f"{org}.{inp}.decoys.fa.gz", nshuf_log)
        models = {"td": model}
        for f in (0, 1):
            for name, (m, _) in oracle[f].items():
                models[f"{name}.f{f}"] = m
        allseq = [tgt[k][3] for k in tids] + [d[3] for d in dec]
        sc = batch_frame_scores(allseq, models, stops)
        nt = len(tids)
        t_meta = [tgt[k][:3] for k in tids]
        d_meta = [d[:3] for d in dec]
        want = set()
        for r in truth_rows:
            o = match_stop(orfs.get(r["transcript_id"], []), r["cds_start"], r["cds_end"], r["strand"])
            if o:
                want.add(o[0])
        for key, v in sc.items():
            q = model_qvalues(t_meta, v[:nt], d_meta, v[nt:], n_shuf)
            model_q[key] = {k: (v[i, 0], round(v[i, 0] - max(v[i, 1], v[i, 2]), 2),
                                None if q[i] != q[i] else q[i]) for i, k in enumerate(tids) if k in want}
            if key == "td":
                # reproduction checks: frame scores vs longest_orfs.cds.scores, q vs decoy_qvalues.tsv
                ds = [abs(v[i, j] - scores[k][j]) for i, k in enumerate(tids) for j in range(3)]
                check["check_td_scores_maxdiff"] = max(ds)
                check["check_td_scores_ndiff_gt0.011"] = sum(x > 0.011 for x in ds)
                same = sum(1 for i, k in enumerate(tids)
                           if (k in qvals) == (q[i] == q[i]) and (k not in qvals or abs(qvals[k] - q[i]) <= 1e-9))
                check["check_td_q_identical"] = same
                check["check_td_q_n"] = nt
                check["n_decoy_orfs"] = len(dec)
                check["n_decoy_shuffles"] = n_shuf

    # training purity: TD training ORFs that are annotated CDSs (same stop)
    truth_by_tid = {r["transcript_id"]: r for r in truth_rows}
    train_ids = set(train)
    pure = 0
    train_hdr_re = re.compile(r"(\S+):(\d+)-(\d+)\(([+-])\)")
    train_coords = {}
    for line in open(wd / "__checkpoints_TDpredict/top_training_orfs.cds"):
        if line.startswith(">"):
            m = train_hdr_re.search(line)
            train_coords[line[1:].split()[0]] = (m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
    for oid, (t, a, b, st) in train_coords.items():
        r = truth_by_tid.get(t)
        if r and st == r["strand"] and abs(stop_end(min(a, b), max(a, b), st) -
                                           stop_end(r["cds_start"], r["cds_end"], r["strand"])) <= 2:
            pure += 1

    rows = []
    for r in truth_rows:
        tid = r["transcript_id"]
        tseq = seqs[tid]
        cds = r["cds_seq"]
        name, desc, seqid = annot.get(tid, ("", "", ""))
        row = dict(org=org, transcript_id=tid, gene_id=r["gene_id"], gene_name=name,
                   description=desc, seqid=seqid, organelle=int(F.is_organelle(seqid)),
                   cds_len=len(cds), len_bin=len_bin(len(cds)), tx_len=r["tx_len"],
                   utr5_len=r["cds_start"] - 1, utr3_len=r["tx_len"] - r["cds_end"],
                   utr5_source=r["utr_src"][0], utr3_source=r["utr_src"][1].strip(),
                   gc=F.gc_content(cds), gc3=r["gc3"], gc3_quintile=gc3_quintile(r["gc3"], gc3_qs),
                   tx_gc=F.gc_content(tseq), fold=r["fold"], L_star=lstar, mu=mu)
        cc = F.codon_counts(cds, stops)
        ncod = sum(cc.values())
        row["codon_llr_train"] = sum(n * codon_llr.get(c, 0) for c, n in cc.items()) / ncod if ncod else float("nan")
        row.update(F.protein_features(tr.translate(cds, frame=0)))
        # truth CDS itself under the TD model (no 5' extension)
        sc_cds = model.frame_scores(cds, stops)
        row["S1_cds"] = sc_cds[0]
        row["C_cds"] = sc_cds[0] - max(sc_cds[1], sc_cds[2])
        # the bare CDS is not a LongOrfs candidate and has no q-value: score part of the test only
        row["pass_scores_cds"] = int(sc_cds[0] > 0 and row["C_cds"] > 0)

        o = match_stop(orfs.get(tid, []), r["cds_start"], r["cds_end"], r["strand"])
        if o is None:
            row.update(orf_id="", category="no_orf", outcome="no_orf")
            rows.append(row)
            continue
        oid, os_, oe, ost = o
        olen = oe - os_ + 1
        oseq = tseq[os_ - 1:oe] if ost == "+" else F_revcomp(tseq[os_ - 1:oe])
        sc = scores[oid]
        q = qvals.get(oid)
        p3, s1pos, cpos, q_ok = passes_rule3(sc, q)
        long_rule = olen >= lstar
        row.update(orf_id=oid, orf_type=orf_types.get(oid, "NA"), orf_start=os_, orf_end=oe,
                   orf_len=olen, ext5=(r["cds_start"] - os_) if ost == "+" else (oe - r["cds_end"]),
                   S1=sc[0], S2=sc[1], S3=sc[2], S4=sc[3], S5=sc[4], S6=sc[5],
                   contrast=sc[0] - max(sc[1], sc[2]), s1_per_base=sc[0] / olen,
                   contrast_rel=(sc[0] - max(sc[1], sc[2])) / (olen * mu), decoy_q=q if q is not None else float("nan"),
                   pass_scores=int(s1pos and cpos),
                   best_other_frame=("S2" if sc[1] >= sc[2] else "S3"),
                   rc_max_gt_S1=int(max(sc[3:6]) > sc[0]),
                   s1_pos=int(s1pos), contrast_pos=int(cpos), decoy_q_ok=int(q_ok),
                   pass_coding=int(p3), long_rule=int(long_rule), candidate=int(p3 or long_rule),
                   in_training=int(oid in train_ids))
        if p3:
            row["category"] = "pass_coding"
        elif long_rule:
            row["category"] = "fail_coding_long_rescued"
        else:
            row["category"] = "fail_coding_lost"
        if not s1pos:
            row["fail_reason"] = "S1<=0"
        elif not cpos:
            row["fail_reason"] = "shifted_frame>=S1"
        elif not q_ok:
            row["fail_reason"] = "decoy_q>0.05"
        else:
            row["fail_reason"] = ""
        # predictions
        tpreds = preds.get(tid, [])
        hit = match_stop([("p", a, b, s) for a, b, s in tpreds], r["cds_start"], r["cds_end"], r["strand"])
        row["predicted_stop"] = int(hit is not None)
        row["predicted_exact"] = int(any(a == r["cds_start"] and b == r["cds_end"] and s == r["strand"]
                                         for a, b, s in tpreds))
        if hit is not None:
            row["outcome"] = "predicted"
        elif row["candidate"]:
            row["outcome"] = "lost_outcompeted"
        else:
            row["outcome"] = "lost_not_candidate"
        if hit is None and tpreds:
            a, b, s = tpreds[0]
            w = match_stop(orfs.get(tid, []), a, b, s)
            row.update(winner_start=a, winner_end=b, winner_strand=s, winner_len=b - a + 1,
                       winner_S1=scores[w[0]][0] if w else float("nan"))
            ov = min(b, r["cds_end"]) - max(a, r["cds_start"]) + 1
            if s != r["strand"]:
                rel = "other_strand"
            elif ov <= 0:
                rel = "no_overlap_utr"
            elif (stop_end(a, b, s) - stop_end(r["cds_start"], r["cds_end"], r["strand"])) % 3 == 0:
                rel = "same_frame"
            else:
                rel = "other_frame_overlap"
            row["winner_relation"] = rel
        # strongest out-of-frame LongOrfs ORF overlapping the CDS
        best_alt = (0.0, float("nan"))
        for aid, a, b, s in orfs.get(tid, []):
            if aid == oid or s != r["strand"]:
                continue
            ov = min(b, r["cds_end"]) - max(a, r["cds_start"]) + 1
            if ov <= 0 or (b - r["cds_end"]) % 3 == 0:
                continue
            frac = ov / len(cds)
            if frac > best_alt[0]:
                best_alt = (frac, scores[aid][0])
        row["altframe_overlap_frac"], row["altframe_S1"] = best_alt
        # homology (LongOrfs >= 50 aa; see homology_source)
        ho = match_stop(h_orfs.get(tid, []), r["cds_start"], r["cds_end"], r["strand"])
        row["homology_orf"] = int(ho is not None)
        if ho is not None:
            row["has_pfam"] = int(ho[0] in h_pfam)
            row["pfam_domains"] = ",".join(h_pfam.get(ho[0], [])[:4])
            sp = h_sprot.get(ho[0])
            row["has_sprot"] = int(sp is not None)
            row["sprot_hit"] = sp[0] if sp else ""
            row["sprot_pident"] = sp[1] if sp else float("nan")
        # alternative Markov models (trained on annotated CDSs of the other fold)
        # (release rule: S1 > 0, C > 0 and q <= 0.05, q from decoys scored with the same model)
        for name in oracle[r["fold"]]:
            if not with_models:
                break
            if name.startswith("gc3q") and name != f"gc3q{row['gc3_quintile']}":
                continue
            key = "gc3match500" if name.startswith("gc3q") else name
            s1m, cm, qm = model_q[f"{name}.f{r['fold']}"][oid]
            row[f"pass_{key}"] = int(s1m > 0 and cm > 0 and qm is not None and qm <= DECOY_Q)
            row[f"S1_{key}"] = s1m
            row[f"C_{key}"] = cm
            row[f"q_{key}"] = qm if qm is not None else float("nan")
        if with_models:
            s1m, cm, qm = model_q["td"][oid]
            row["q_td_recomputed"] = qm if qm is not None else float("nan")
        rows.append(row)

    meta = dict(org=org, genetic_code=gcode, transcriptome_gc_pct=gc_log, L_star=lstar, mu=mu, mu_log=mu_log,
                run=RUN, decoy_q=DECOY_Q, homology_source=homology_source(org, inp), **check,
                n_truth=len(truth_rows), n_training=len(train), training_is_truth_cds=pure,
                training_gc3_mean=sum(train_gc3) / len(train_gc3),
                truth_gc3_mean=sum(r["gc3"] for r in truth_rows) / len(truth_rows),
                training_len_median=sorted(len(s) for s in train.values())[len(train) // 2],
                truth_cds_len_median=sorted(len(r["cds_seq"]) for r in truth_rows)[len(truth_rows) // 2],
                gc3_quintile_breaks=",".join(f"{q:.3f}" for q in gc3_qs))
    for f in (0, 1):
        for name, (m, mu_m) in oracle[f].items():
            meta[f"mu_{name}_fold{f}"] = mu_m
    tag = org if inp == "complete" else f"{org}.{inp}"
    meta["input"] = inp
    write_tsv(out_dir / "per_cds" / f"{tag}.per_cds.tsv.gz", rows)
    write_tsv(out_dir / "per_cds" / f"{tag}.meta.tsv", [meta])
    return org, len(rows)


def F_revcomp(s):
    from markov import revcomp
    return revcomp(s)


def write_tsv(path, rows, cols=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if cols is None:
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")


def fmt(v):
    if isinstance(v, float):
        return "NA" if v != v else ("inf" if v == float("inf") else f"{v:.6g}")
    return str(v)


# ------------------------------------------------------------------- exploratory: alternative rules

ALT_INPUTS = [(o, i) for o in ORGS for i in ("complete.wLnc", "complete.wRand")]


def alt_rules_one(org, inp, out_dir):
    """Exploratory: what would dropping the S1 > 0 requirement cost on noncoding negatives, under the
    release rule? For a *.wLnc / *.wRand input with a pytd_v7 run, count transcripts with at least one
    candidate ORF (any ORF passing the rule; homology off) under
      current: S1 > 0, C > 0 and q <= 0.05 (q from decoy_qvalues.tsv), or length >= L*
      A:       C > 0 and q_A <= 0.05, or length >= L*   (S1 > 0 dropped; q_A recomputed with the run's
               decoys and model, eligibility C > 0 only)
      B:       current, or C >= 0.5 * mu * length        (mu: training median per-base contrast)
    for negatives (transcripts not in the truth table) and, for positives, whether the truth ORF
    (stop within 2 nt) passes. Candidate status is an upper bound on calls: Predict keeps one ORF per
    transcript, ranked by S1."""
    run = P / f"work/bench_v3/runs/{org}/{inp}/{RUN}"
    if not (run / "preds.tsv").exists():
        return None
    from pytransdecoder.core.translator import Translator
    gcode = orgs_table()[org]["td_genetic_code"]
    stops = set(Translator(gcode).get_stop_codons())
    wd = run / "input.fa.transdecoder_dir"
    lstar, mu, _, nshuf_log = parse_predict_log(run / "predict.log")
    qvals = read_qvalues(wd / "decoy_qvalues.tsv")
    orfs = read_orf_gff3(wd / "longest_orfs.gff3")
    scores = read_scores(wd / "longest_orfs.cds.scores")
    seqs = dict(read_fasta(run / "input.fa"))
    tgt = target_orfs(wd, seqs)
    tids = sorted(tgt)
    n_shuf, dec = decoy_orfs(run, out_dir / "decoys" / f"{org}.{inp}.decoys.fa.gz", nshuf_log)
    model = MarkovModel.from_scores_file(wd / "hexamer.scores")
    sc = batch_frame_scores([tgt[k][3] for k in tids] + [d[3] for d in dec], {"td": model}, stops)["td"]
    nt = len(tids)
    t_meta, d_meta = [tgt[k][:3] for k in tids], [d[:3] for d in dec]
    q_cur = model_qvalues(t_meta, sc[:nt], d_meta, sc[nt:], n_shuf, require_s1=True)
    q_a = model_qvalues(t_meta, sc[:nt], d_meta, sc[nt:], n_shuf, require_s1=False)
    q_cur = {k: q_cur[i] for i, k in enumerate(tids)}
    q_a = {k: q_a[i] for i, k in enumerate(tids)}
    q_same = sum(1 for k in tids if (k in qvals) == (q_cur[k] == q_cur[k])
                 and (k not in qvals or abs(qvals[k] - q_cur[k]) <= 1e-9))
    truth = {}
    for line in open(P / f"data/bench_v3/inputs/{org}/{inp}.truth.tsv").read().splitlines()[1:]:
        p = line.split("\t")
        truth[p[0]] = (int(p[1]), int(p[2]), p[3])
    ids = list(seqs)

    def rules(o):
        oid, a, b, _ = o
        L = b - a + 1
        sc_ = scores[oid]
        c = sc_[0] - max(sc_[1], sc_[2])
        q = qvals.get(oid)
        qa = q_a.get(oid, float("nan"))
        cur = (sc_[0] > 0 and c > 0 and q is not None and q <= DECOY_Q) or L >= lstar
        alt_a = (c > 0 and qa == qa and qa <= DECOY_Q) or L >= lstar
        alt_b = cur or (mu is not None and c >= 0.5 * mu * L)
        return cur, alt_a, alt_b
    neg = [0, 0, 0]
    pos = [0, 0, 0]
    npos = nneg = 0
    for tid in ids:
        if tid in truth:
            s_, e_, st = truth[tid]
            o = match_stop(orfs.get(tid, []), s_, e_, st)
            npos += 1
            if o:
                for i, v in enumerate(rules(o)):
                    pos[i] += v
        else:
            nneg += 1
            flags = [False, False, False]
            for o in orfs.get(tid, []):
                for i, v in enumerate(rules(o)):
                    flags[i] |= v
            for i in range(3):
                neg[i] += flags[i]
    return dict(org=org, input=inp, L_star=lstar, mu=mu, n_pos=npos, n_neg=nneg,
                pos_pass_current=pos[0], pos_pass_A_no_S1=pos[1], pos_pass_B_contrast=pos[2],
                neg_with_candidate_current=neg[0], neg_with_candidate_A_no_S1=neg[1],
                neg_with_candidate_B_contrast=neg[2], n_decoy_orfs=len(dec), n_decoy_shuffles=n_shuf,
                check_q_identical=q_same, check_q_n=nt)


def alt_rules(out_dir, jobs=2):
    with Pool(jobs) as pool:
        rows = [r for r in pool.starmap(alt_rules_one, [(o, i, out_dir) for o, i in ALT_INPUTS]) if r]
    write_tsv(out_dir / "alt_rules_negatives.tsv", rows)
    return len(rows)


def length_retained(out_dir):
    """Final predictions whose ORF fails the coding test and is kept by the long-ORF rule alone:
    bench_v3 complete / wLnc / wRand inputs (release runs). GM12878 is not included in this copy."""
    def final_ids(gff3):
        out = []
        for line in open(gff3):
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "mRNA":
                out.append((re.search(r"ID=([^;]+)", f[8]).group(1), f[0]))
        return out

    def classify(wd, run, gff3):
        lstar = parse_predict_log(run / "predict.log")[0]
        qvals = read_qvalues(wd / "decoy_qvalues.tsv")
        orfs = {o[0]: o for v in read_orf_gff3(wd / "longest_orfs.gff3").values() for o in v}
        scores = read_scores(wd / "longest_orfs.cds.scores")
        res = []
        for oid, tid in final_ids(gff3):
            o = orfs[oid]
            L = o[2] - o[1] + 1
            p3 = passes_rule3(scores[oid], qvals.get(oid))[0]
            res.append((oid, tid, o, p3, L >= lstar))
        return res, lstar

    rows = []
    for org in ORGS:
        for inp in ("complete", "complete.wLnc", "complete.wRand"):
            run = P / f"work/bench_v3/runs/{org}/{inp}/{RUN}"
            if not (run / "input.fa.transdecoder.gff3").exists():
                continue
            truth = {}
            for line in open(P / f"data/bench_v3/inputs/{org}/{inp}.truth.tsv").read().splitlines()[1:]:
                p = line.split("\t")
                truth[p[0]] = (int(p[1]), int(p[2]), p[3])
            res, lstar = classify(run / "input.fa.transdecoder_dir", run, run / "input.fa.transdecoder.gff3")
            lo = [r for r in res if not r[3]]
            n_truth = n_neg = n_other = 0
            for oid, tid, o, p3, lr in lo:
                if tid not in truth:
                    n_neg += 1
                elif match_stop([o], *truth[tid]):
                    n_truth += 1
                else:
                    n_other += 1
            rows.append(dict(dataset=f"bench_v3:{org}:{inp}", L_star=lstar, predictions=len(res),
                             length_rule_only=len(lo), length_rule_only_pct=round(100 * len(lo) / len(res), 3),
                             lo_is_truth_cds=n_truth, lo_on_negative_transcript=n_neg,
                             lo_other_orf_on_coding_transcript=n_other))
    write_tsv(out_dir / "length_rule_only_predictions.tsv", rows)
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["bench", "altrules", "summarize", "all"])
    ap.add_argument("--orgs", nargs="+", default=ORGS)
    ap.add_argument("--input", default="complete", choices=["complete", "complete_fullutr"],
                    help="bench_v3 input set (complete: UTRs trimmed to <= 250 nt; complete_fullutr: full UTRs)")
    ap.add_argument("--jobs", type=int, default=2, help="processes (at most 2)")
    ap.add_argument("--out", type=Path, default=P / "work/long_orf_failures_v7")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    jobs = max(1, min(args.jobs, 2))
    if args.step in ("bench", "all"):
        inputs = ["complete", "complete_fullutr"] if args.step == "all" else [args.input]
        with Pool(min(jobs, len(args.orgs))) as pool:
            for inp in inputs:
                # truth-trained models only on the complete inputs (the full-UTR set is used for the
                # UTR comparison, which needs the categories only)
                for org, n in pool.starmap(process_org, [(o, args.out, inp, inp == "complete") for o in args.orgs]):
                    print(f"{org} {inp}: {n} truth CDSs", file=sys.stderr)
    if args.step in ("altrules", "all"):
        print(f"alternative rules: {alt_rules(args.out, jobs)} inputs", file=sys.stderr)
        print(f"length-rule-only predictions: {length_retained(args.out)} datasets", file=sys.stderr)
    if args.step in ("summarize", "all"):
        import importlib.util
        spec = importlib.util.spec_from_file_location("summarize_v7", HERE / "05_summarize_v7.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.main(args.out)


if __name__ == "__main__":
    main()
