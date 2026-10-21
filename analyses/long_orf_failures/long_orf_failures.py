#!/usr/bin/env python3
"""How often do true, annotated CDSs fail TransDecoder's coding criteria, and why?

For each bench_v3 organism (complete transcripts, strand-specific), every truth CDS is matched to the
LongOrfs ORF with the same strand and a stop within 2 nt. Using the final-default PyTD run
(pytd_dev_v2 = PY_dev 48c870c defaults) the ORF is classified by the Predict selection rules
(predict.py _select_best_orfs, no homology):

  rule 3 (coding)  S1 > 0 and C = S1 - max(S2, S3) > 0, and, if the ORF is < 300 nt, C >= 200 * mu
  rule 2 (length)  ORF length >= L*(GC)  (disabled for non-TAA/TAG/TGA codes, e.g. Tetrahymena)

and by outcome (a final prediction with the same stop, or not). Features of each CDS (GC, GC3, codon
usage relative to the training set, amino-acid composition and low-complexity measures, Pfam and
Swiss-Prot hits from the old pytd_homology runs, annotation names) are recorded, and each ORF is
rescored with alternative Markov models trained on annotated CDSs (two-fold, by gene) to test
whether failures come from the training set or from the sequences themselves.

The GM12878 part does the same for MANE Select CDSs fully contained in an assembled StringTie or
LRAA transcript (eval_genome_cds.py 'assembled'), using the final-default ab initio runs.

usage:
  long_orf_failures.py bench [--orgs human mouse ...] [--jobs 4]   # per-CDS tables
  long_orf_failures.py gm12878                                     # GM12878 MANE tables
  long_orf_failures.py altrules                                    # exploratory: rules without S1 > 0; length-only calls
  long_orf_failures.py summarize                                   # summary tables from the above
  long_orf_failures.py all [--jobs 4]       # bench for complete and complete_fullutr, then the rest
Outputs go to results/ next to this script (override with --out).
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
sys.path.insert(0, str(HERE))
P = Path(os.environ.get("TD_PAPER_ROOT", "/broad/hptmp/bhaas/TRANSDECODER_PAPER"))
T = P / "TransDecoderPaper"
PYTD = Path(os.environ.get("PYTD_TREE", P / "work/trees/PY_dev"))
sys.path.insert(0, str(PYTD))
sys.path.insert(0, str(T / "genome_annotation/scripts"))

from markov import MarkovModel, median_contrast, read_bg  # noqa: E402
import features as F  # noqa: E402

ORGS = ["athal", "creinhardtii", "dmel", "human", "mouse", "pfalciparum", "spombe", "tthermophila"]
SHORT_ORF_NT = 300
SHORT_ORF_MIN_CONTRAST = 200.0
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
    txt = open(path).read()
    m = re.search(r"min auto-accept length = (\S+) nt", txt)
    lstar = float(m.group(1)) if m else float("inf")
    if "no length-based auto-accept" in txt:
        lstar = float("inf")
    m = re.search(r"\(200 bases x ([0-9.]+)\)", txt)
    mu = float(m.group(1)) if m else None
    m = re.search(r"Transcriptome GC content: ([0-9.]+)%", txt)
    gc = float(m.group(1)) if m else float("nan")
    return lstar, mu, gc


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


def passes_rule3(sc, orf_len, mu):
    s1, c = sc[0], sc[0] - max(sc[1], sc[2])
    short_ok = orf_len >= SHORT_ORF_NT or (mu is not None and c >= SHORT_ORF_MIN_CONTRAST * mu)
    return s1 > 0 and c > 0 and short_ok, s1 > 0, c > 0, short_ok


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


def homology_lookup(org, tids, inp="complete"):
    """(transcript, stop) -> dict(pfam=..., sprot=..., pident=...) from the old pytd_homology run
    (LongOrfs ORFs >= 100 aa); matched by coordinates, not IDs."""
    d = P / f"work/bench_v3/runs/{org}/{inp}/pytd_homology"
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


def process_org(org, out_dir, inp="complete"):
    orgs = orgs_table()
    gcode = orgs[org]["td_genetic_code"]
    from pytransdecoder.core.translator import Translator
    tr = Translator(gcode)
    stops = set(tr.get_stop_codons())
    run = P / f"work/bench_v3/runs/{org}/{inp}/pytd_dev_v2"
    wd = run / "input.fa.transdecoder_dir"
    lstar, mu_log, gc_log = parse_predict_log(run / "predict.log")
    seqs = dict(read_fasta(run / "input.fa"))
    orfs = read_orf_gff3(wd / "longest_orfs.gff3")
    orf_types = read_fasta_headers(wd / "longest_orfs.cds")
    scores = read_scores(wd / "longest_orfs.cds.scores")
    train = dict(read_fasta(wd / "__checkpoints_TDpredict/top_training_orfs.cds"))
    preds = read_preds_tsv(run / "preds.tsv")
    model = MarkovModel.from_scores_file(wd / "hexamer.scores")
    bg = read_bg(wd / "base_freqs.dat")
    mu = median_contrast(model, list(train.values()), stops)
    assert mu_log is None or abs(mu - mu_log) < 1e-3, (mu, mu_log)
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
        row["pass_coding_cds"] = int(passes_rule3(sc_cds, len(cds), mu)[0])

        o = match_stop(orfs.get(tid, []), r["cds_start"], r["cds_end"], r["strand"])
        if o is None:
            row.update(orf_id="", category="no_orf", outcome="no_orf")
            rows.append(row)
            continue
        oid, os_, oe, ost = o
        olen = oe - os_ + 1
        oseq = tseq[os_ - 1:oe] if ost == "+" else F_revcomp(tseq[os_ - 1:oe])
        sc = scores[oid]
        p3, s1pos, cpos, short_ok = passes_rule3(sc, olen, mu)
        long_rule = olen >= lstar
        row.update(orf_id=oid, orf_type=orf_types.get(oid, "NA"), orf_start=os_, orf_end=oe,
                   orf_len=olen, ext5=(r["cds_start"] - os_) if ost == "+" else (oe - r["cds_end"]),
                   S1=sc[0], S2=sc[1], S3=sc[2], S4=sc[3], S5=sc[4], S6=sc[5],
                   contrast=sc[0] - max(sc[1], sc[2]), s1_per_base=sc[0] / olen,
                   contrast_rel=(sc[0] - max(sc[1], sc[2])) / (olen * mu),
                   best_other_frame=("S2" if sc[1] >= sc[2] else "S3"),
                   rc_max_gt_S1=int(max(sc[3:6]) > sc[0]),
                   s1_pos=int(s1pos), contrast_pos=int(cpos), short_guard_ok=int(short_ok),
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
        elif not short_ok:
            row["fail_reason"] = "short_orf_guard"
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
        # homology (old pytd_homology LongOrfs, >= 100 aa)
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
        for name, (m, mu_m) in oracle[r["fold"]].items():
            if name.startswith("gc3q") and name != f"gc3q{row['gc3_quintile']}":
                continue
            key = "gc3match500" if name.startswith("gc3q") else name
            s3 = [m.score(oseq, stops), m.score(oseq[1:], stops), m.score(oseq[2:], stops)]
            row[f"pass_{key}"] = int(passes_rule3(s3, olen, mu_m)[0])
            row[f"S1_{key}"] = s3[0]
            row[f"C_{key}"] = s3[0] - max(s3[1], s3[2])
        rows.append(row)

    meta = dict(org=org, genetic_code=gcode, transcriptome_gc_pct=gc_log, L_star=lstar, mu=mu,
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


# ------------------------------------------------------------------------------------ GM12878

GM_ASSEMBLIES = {
    "stringtie": "data/gm12878/star_stringtie/stringtie.gtf",
    "lraa": "data/gm12878/longread_lraa/lraa_denovo/LRAA.gtf",
}


def gm12878(out_dir):
    import eval_genome_cds as E
    from pytransdecoder.genome import read_structures
    from pytransdecoder.core.translator import Translator
    tr = Translator("Universal")
    stops = set(tr.get_stop_codons())
    mane_gtf = P / "data/annot/human/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
    mane_all, _ = E.read_ref_gtf(mane_gtf, None)
    mane = {t: r for t, r in mane_all.items() if "MANE_Select" in r["tags"]}
    # bench human per-CDS table (same MANE transcripts) for concordance
    bench = {}
    bp = out_dir / "per_cds/human.per_cds.tsv.gz"
    if bp.exists():
        with gzip.open(bp, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                bench[d["transcript_id"]] = d
    rows = []
    for asm, gtf in GM_ASSEMBLIES.items():
        run = P / f"work/gm12878_td/{asm}/final_abinitio"
        wd = run / "transcripts.fa.transdecoder_dir"
        lstar, mu, gc = parse_predict_log(run / "predict.log")
        status = {}
        for line in open(P / f"work/gm12878_td/eval/status/{asm}_final_abinitio.ref_status.tsv"):
            p = line.rstrip("\n").split("\t")
            if p[0] == "MANE":
                status[p[1]] = (p[2], int(p[3]), int(p[4]), int(p[5]))
        structs = read_structures(P / gtf)
        assembly = {tid: dict(chrom=s.seqid, strand=s.strand, segs=tuple(sorted(s.exons)))
                    for tid, s in structs.items() if s.strand in "+-" and s.exons}
        idx = E.Index(assembly)
        orfs = read_orf_gff3(wd / "longest_orfs.gff3")
        scores = read_scores(wd / "longest_orfs.cds.scores")
        train_ids = {line[1:].split()[0] for line in open(wd / "__checkpoints_TDpredict/top_training_orfs.cds")
                     if line.startswith(">")}
        preds = defaultdict(list)
        for line in open(run / "transcripts.fa.transdecoder.gff3"):
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "CDS":
                preds[f[0]].append((int(f[3]), int(f[4]), f[6]))
        need = set()
        cand = {}
        for tid, (gname, asmb, exact, comp) in status.items():
            if not asmb or tid not in mane:
                continue
            r = mane[tid]
            conts = [k for k in idx.overlapping(r) if E.contained(r["segs"], assembly[k]["segs"]) is not None]
            cand[tid] = conts
            need.update(conts)
        seqs = {k: v for k, v in read_fasta(run / "transcripts.fa") if k in need}
        order = ["predicted_exact", "predicted_stop", "lost_outcompeted", "fail_coding_lost", "no_orf"]
        for tid, conts in cand.items():
            gname, asmb, exact, comp = status[tid]
            r = mane[tid]
            g_stop = r["segs"][-1][1] if r["strand"] == "+" else r["segs"][0][0]
            g_start = r["segs"][0][0] if r["strand"] == "+" else r["segs"][-1][1]
            best = None
            for k in conts:
                s = structs[k]
                t_stop, t_start = to_tx(s, g_stop), to_tx(s, g_start)
                if t_stop is None or t_start is None:
                    continue
                ost = "+" if s.strand == r["strand"] else "-"
                lo, hi = min(t_start, t_stop), max(t_start, t_stop)
                o = match_stop(orfs.get(k, []), lo, hi, ost)
                info = dict(asm_transcript=k, cds_t_start=lo, cds_t_end=hi)
                tseq = seqs.get(k, "")
                cds = tseq[lo - 1:hi] if ost == "+" else F_revcomp(tseq[lo - 1:hi])
                info["cds_len"] = len(cds)
                if o is None:
                    st = "no_orf"
                else:
                    oid, a, b, _ = o
                    olen = b - a + 1
                    sc = scores[oid]
                    p3, s1pos, cpos, short_ok = passes_rule3(sc, olen, mu)
                    info.update(orf_id=oid, orf_len=olen, S1=sc[0], S2=sc[1], S3=sc[2],
                                contrast=sc[0] - max(sc[1], sc[2]),
                                contrast_rel=(sc[0] - max(sc[1], sc[2])) / (olen * mu),
                                pass_coding=int(p3), long_rule=int(olen >= lstar),
                                in_training=int(oid in train_ids),
                                fail_reason="" if p3 else ("S1<=0" if not s1pos else
                                                           "shifted_frame>=S1" if not cpos else "short_orf_guard"))
                    tp = preds.get(k, [])
                    if any(a2 == lo and b2 == hi and s2 == ost for a2, b2, s2 in tp):
                        st = "predicted_exact"
                    elif match_stop([("p", a2, b2, s2) for a2, b2, s2 in tp], lo, hi, ost):
                        st = "predicted_stop"
                    elif p3 or olen >= lstar:
                        st = "lost_outcompeted"
                        if tp:
                            a2, b2, s2 = tp[0]
                            info["winner"] = f"{a2}-{b2}({s2})"
                            w = match_stop(orfs.get(k, []), a2, b2, s2)
                            info["winner_S1"] = scores[w[0]][0] if w else float("nan")
                    else:
                        st = "fail_coding_lost"
                info["status"] = st
                if best is None or order.index(st) < order.index(best["status"]):
                    best = info
            row = dict(assembly=asm, transcript_id=tid, gene_name=gname, eval_exact=exact,
                       eval_compatible=comp, n_containing=len(conts), L_star=lstar, mu=mu)
            if best:
                row.update(best)
                tseq = seqs.get(best["asm_transcript"], "")
                cds = tseq[best["cds_t_start"] - 1:best["cds_t_end"]]
                if cds:
                    row.update(gc=F.gc_content(cds), gc3=F.gc3(cds))
                    pf = F.protein_features(tr.translate(cds, frame=0))
                    row.update({k: pf[k] for k in ("aa_len", "aa_entropy", "top1_aa", "top1_frac",
                                                   "top3_aa", "seg_lowc_frac", "disorder_frac",
                                                   "tandem3_frac", "wmer_unique_pct")})
            b = bench.get(tid)
            if b:
                row.update(bench_category=b.get("category", ""), bench_outcome=b.get("outcome", ""),
                           bench_has_pfam=b.get("has_pfam", ""), bench_has_sprot=b.get("has_sprot", ""))
            rows.append(row)
    write_tsv(out_dir / "gm12878_mane_assembled.tsv.gz", rows)
    return len(rows)


# ------------------------------------------------------------------- exploratory: alternative rules

def alt_rules(out_dir):
    """Exploratory: what would dropping the S1 > 0 requirement cost on noncoding negatives?
    For each *.wLnc / *.wRand input with a pytd_dev_v2 run, count transcripts with at least one
    candidate ORF (any ORF passing the rule; homology off) under
      current: S1 > 0 and C > 0 (+ short-ORF guard) or length >= L*
      A:       C > 0 (+ short-ORF guard) or length >= L*            (S1 > 0 dropped)
      B:       current, or C >= 0.5 * mu * length                     (strong frame contrast)
    for negatives (transcripts not in the truth table) and, for positives, whether the truth ORF
    (stop within 2 nt) passes. Candidate status is an upper bound on calls: Predict keeps one ORF per
    transcript, ranked by S1."""
    rows = []
    for org in ORGS:
        for inp in ("complete.wLnc", "complete.wRand"):
            run = P / f"work/bench_v3/runs/{org}/{inp}/pytd_dev_v2"
            if not (run / "preds.tsv").exists():
                continue
            wd = run / "input.fa.transdecoder_dir"
            lstar, mu, _ = parse_predict_log(run / "predict.log")
            orfs = read_orf_gff3(wd / "longest_orfs.gff3")
            scores = read_scores(wd / "longest_orfs.cds.scores")
            truth = {}
            for line in open(P / f"data/bench_v3/inputs/{org}/{inp}.truth.tsv").read().splitlines()[1:]:
                p = line.split("\t")
                truth[p[0]] = (int(p[1]), int(p[2]), p[3])
            ids = [line[1:].split()[0] for line in open(run / "input.fa") if line.startswith(">")]

            def rules(o):
                oid, a, b, _ = o
                L = b - a + 1
                sc = scores[oid]
                c = sc[0] - max(sc[1], sc[2])
                short_ok = L >= SHORT_ORF_NT or c >= SHORT_ORF_MIN_CONTRAST * mu
                cur = (sc[0] > 0 and c > 0 and short_ok) or L >= lstar
                alt_a = (c > 0 and short_ok) or L >= lstar
                alt_b = cur or c >= 0.5 * mu * L
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
            rows.append(dict(org=org, input=inp, L_star=lstar, mu=mu, n_pos=npos, n_neg=nneg,
                             pos_pass_current=pos[0], pos_pass_A_no_S1=pos[1], pos_pass_B_contrast=pos[2],
                             neg_with_candidate_current=neg[0], neg_with_candidate_A_no_S1=neg[1],
                             neg_with_candidate_B_contrast=neg[2]))
    write_tsv(out_dir / "alt_rules_negatives.tsv", rows)
    return len(rows)


def length_retained(out_dir):
    """Final predictions whose ORF fails the coding test and is kept by the long-ORF rule alone
    (the count E10 lists as missing): bench_v3 complete / wLnc / wRand inputs and GM12878."""
    def final_ids(gff3):
        out = []
        for line in open(gff3):
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "mRNA":
                out.append((re.search(r"ID=([^;]+)", f[8]).group(1), f[0]))
        return out

    def classify(wd, run, gff3):
        lstar, mu, _ = parse_predict_log(run / "predict.log")
        orfs = {o[0]: o for v in read_orf_gff3(wd / "longest_orfs.gff3").values() for o in v}
        scores = read_scores(wd / "longest_orfs.cds.scores")
        res = []
        for oid, tid in final_ids(gff3):
            o = orfs[oid]
            L = o[2] - o[1] + 1
            p3 = passes_rule3(scores[oid], L, mu)[0]
            res.append((oid, tid, o, p3, L >= lstar))
        return res, lstar

    rows = []
    for org in ORGS:
        for inp in ("complete", "complete.wLnc", "complete.wRand"):
            run = P / f"work/bench_v3/runs/{org}/{inp}/pytd_dev_v2"
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
    for asm in ("stringtie", "lraa", "trinity"):
        run = P / f"work/gm12878_td/{asm}/final_abinitio"
        if not run.exists():
            continue
        res, lstar = classify(run / "transcripts.fa.transdecoder_dir", run, run / "transcripts.fa.transdecoder.gff3")
        status = {}
        for line in open(P / f"work/gm12878_td/eval/status/{asm}_final_abinitio.pred_status.tsv"):
            p = line.rstrip("\n").split("\t")
            if p[0] != "model_id":
                status[p[0]] = p
        lo = [r for r in res if not r[3]]
        st = [status[r[0]] for r in lo if r[0] in status]
        rows.append(dict(dataset=f"gm12878:{asm}:final_abinitio", L_star=lstar, predictions=len(res),
                         length_rule_only=len(lo), length_rule_only_pct=round(100 * len(lo) / len(res), 3),
                         lo_distinct_in_eval=len(st), lo_mane_exact=sum(int(x[6]) for x in st),
                         lo_gencode_compatible=sum(int(x[4]) for x in st),
                         lo_no_gencode_cds_overlap=sum(int(x[5]) for x in st)))
    write_tsv(out_dir / "length_rule_only_predictions.tsv", rows)
    return len(rows)


def to_tx(struct, g):
    for g0, g1, t0, t1 in struct.blocks:
        if g0 <= g <= g1:
            return t0 + (g - g0) if t1 >= t0 else t0 - (g - g0)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["bench", "gm12878", "altrules", "summarize", "all"])
    ap.add_argument("--orgs", nargs="+", default=ORGS)
    ap.add_argument("--input", default="complete", choices=["complete", "complete_fullutr"],
                    help="bench_v3 input set (complete: UTRs trimmed to <= 250 nt; complete_fullutr: full UTRs)")
    ap.add_argument("--jobs", type=int, default=4, help="processes (the machine is shared: <= 4)")
    ap.add_argument("--out", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.step in ("bench", "all"):
        inputs = ["complete", "complete_fullutr"] if args.step == "all" else [args.input]
        with Pool(min(args.jobs, 4, len(args.orgs))) as pool:
            for inp in inputs:
                for org, n in pool.starmap(process_org, [(o, args.out, inp) for o in args.orgs]):
                    print(f"{org} {inp}: {n} truth CDSs", file=sys.stderr)
    if args.step in ("gm12878", "all"):
        print(f"GM12878: {gm12878(args.out)} MANE x assembly rows", file=sys.stderr)
    if args.step in ("altrules", "all"):
        print(f"alternative rules: {alt_rules(args.out)} inputs", file=sys.stderr)
        print(f"length-rule-only predictions: {length_retained(args.out)} datasets", file=sys.stderr)
    if args.step in ("summarize", "all"):
        import summarize
        summarize.main(args.out)


if __name__ == "__main__":
    main()
