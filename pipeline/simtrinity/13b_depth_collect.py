#!/usr/bin/env python3
"""Collect a read-depth ablation into small tables (E31 layout; generalised for E32).

usage: 13b_depth_collect.py [--data-root D --work-root W] <outdir> <level> [<level> ...]
  Default roots are the human E31 run ($P/data/simtrinity, $P/work/simtrinity), where d1 is the E30 run.
  With --data-root/--work-root (e.g. $P/data/simtrinity_athal, $P/work/simtrinity_athal, E32) every level,
  including d1, is read from <root>/depth/<level>/, and the E29 statistic is optional (written as NA if absent).
  Same outputs as 13_depth_collect.py.

original usage (E31): 13_depth_collect.py <outdir> <level> [<level> ...]
  levels: d1 (= the full-depth E30 run) and d2, d4, ... (1/2, 1/4, ... of the read pairs)
Inputs per level (P = /broad/hptmp/bhaas/TRANSDECODER_PAPER):
  d1:  $P/data/simtrinity/labels/sim_trinity.*, $P/data/simtrinity/trinity/, scores in data/simtrinity/,
       filtered predictions in $P/work/simtrinity/filtered_preds/
  dN:  $P/data/simtrinity/depth/dN/{labels,pair_counts.tsv,Trinity.fasta},
       $P/work/simtrinity/depth/dN/scores/ (score_simtrinity.py output, incl. _filtered_preds/)
Writes to <outdir>:
  levels.tsv            level, fraction of read pairs, pairs (observed), Trinity transcripts
  label_counts.tsv      level x label x reason: transcripts
  truth_counts.tsv      level: scored truth by type, negatives, ambiguous, short, unresolved, total, and each
                        type's share of the scored truth (share_*)
  scores_overall.tsv    level x tool: score_predictions.py metrics (+ n_pred_excluded)
  scores_by_type.tsv    level x tool x type: n, tp_stop, tp_exact, sensitivity, and precision / F1 over the
                        (scored) predictions on transcripts of that type; *_wneg columns also count, in the
                        precision denominator, the predictions on negatives whose own ORF is of that type
                        (ATG start / stop codon end, from the transcript sequence; same rule for every tool)
  neg_pred_types.tsv    level x tool: predictions on negatives by their own ORF type
  scores_by_expr.tsv    level x tool x TPM bin (score_simtrinity.py)
  negatives.tsv         level x tool: negatives with a prediction
  excluded_preds.tsv    level x tool x excluded label
  e29_five_prime.tsv    level: transcripts with a MANE stop + reference structure to their 5' end, share
                        beginning inside the CDS (03_utr_gm12878_truth.py)
"""
import collections
import gzip
import os
import sys

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
TOL = 2


DATA_ROOT = WORK_ROOT = None


def paths(level):
    if DATA_ROOT:
        d = f"{DATA_ROOT}/depth/{level}"
        w = f"{WORK_ROOT}/depth/{level}/scores"
        return dict(labels=f"{d}/labels/sim_trinity", fasta=f"{d}/Trinity.fasta", scores=w, fpreds=f"{w}/_filtered_preds",
                    e29=f"{d}/labels/sim_trinity_e29.pairs.tsv.gz",
                    pairs=f"{d}/pair_counts.tsv" if level != "d1" else None)
    if level == "d1":
        return dict(labels=f"{P}/data/simtrinity/labels/sim_trinity", fasta=f"{P}/data/simtrinity/trinity/Trinity.fasta",
                    scores=os.path.join(T, "data", "simtrinity"), fpreds=f"{P}/work/simtrinity/filtered_preds",
                    e29=f"{P}/data/simtrinity/labels/sim_trinity_e29.pairs.tsv.gz", pairs=None)
    d = f"{P}/data/simtrinity/depth/{level}"
    w = f"{P}/work/simtrinity/depth/{level}/scores"
    return dict(labels=f"{d}/labels/sim_trinity", fasta=f"{d}/Trinity.fasta", scores=w, fpreds=f"{w}/_filtered_preds",
                e29=f"{d}/labels/sim_trinity_e29.pairs.tsv.gz", pairs=f"{d}/pair_counts.tsv")


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return hdr, [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]


STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path, keep):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in keep:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in keep:
                buf.append(line.strip())
    if name in keep:
        seqs[name] = "".join(buf).upper()
    return seqs


def orf_type(seq, s, e, strand):
    """type of a predicted ORF from its own ends: starts with ATG / ends with a stop codon (preds.tsv
    coordinates include the stop when the tool reports one)"""
    c = seq[s - 1:e]
    if strand == "-":
        c = c.translate(COMP)[::-1]
    st, sp = c[:3] == "ATG", c[-3:] in STOPS
    return {(True, True): "complete", (False, True): "5prime_partial",
            (True, False): "3prime_partial", (False, False): "internal"}[(st, sp)]


def ends(s, e, strand):
    return (s, e) if strand == "+" else (e, s)


def main():
    global DATA_ROOT, WORK_ROOT
    argv = sys.argv[1:]
    while argv and argv[0].startswith("--"):
        if argv[0] == "--data-root":
            DATA_ROOT = argv[1]
        elif argv[0] == "--work-root":
            WORK_ROOT = argv[1]
        argv = argv[2:]
    outdir, levels = argv[0], argv[1:]
    os.makedirs(outdir, exist_ok=True)
    out = collections.defaultdict(list)
    hdrs = {}
    full_pairs = None
    for level in levels:
        p = paths(level)
        ntx = sum(1 for line in open(p["fasta"]) if line.startswith(">"))
        if p["pairs"]:
            _, r = read_tsv(p["pairs"])
            frac, pairs, full_pairs = r[0]["fraction"], r[0]["pairs_sub"], r[0]["pairs_full"]
            obs = r[0]["observed_fraction"]
        else:
            frac, pairs, obs = "1", "", "1"
        out["levels"].append([level, frac, pairs, obs, ntx])
        # labels
        _, labs = read_tsv(p["labels"] + ".labels.tsv")
        c = collections.Counter((r["label"], r["reason"] if r["label"] in ("unresolved", "ambiguous") else ".") for r in labs)
        for (l, r), n in sorted(c.items()):
            out["label_counts"].append([level, l, r, n])
        cl = collections.Counter(r["label"] for r in labs)
        tot = sum(cl[t] for t in TYPES)
        out["truth_counts"].append([level] + [cl[t] for t in TYPES] + [tot, cl["noncoding"],
                                   cl["ambiguous"], cl["short_complete"] + cl["short_partial"], cl["unresolved"], len(labs)]
                                   + [f"{cl[t] / tot:.6f}" for t in TYPES])
        # scores from score_simtrinity.py
        for name in ("scores_overall", "scores_by_expr", "negatives", "excluded_preds"):
            h, rows = read_tsv(os.path.join(p["scores"], name + ".tsv"))
            hdrs[name] = h
            for r in rows:
                out[name].append([level] + [r[k] for k in h])
        # by type with per-type precision, from the filtered (scored) predictions
        _, truth = read_tsv(p["labels"] + ".truth.tsv")
        tt = {r["transcript_id"]: r for r in truth}
        negs = {l.strip() for l in open(p["labels"] + ".negatives.txt") if l.strip()}
        seqs = read_fasta(p["fasta"], negs)
        for tool in [r["tool"] for r in read_tsv(os.path.join(p["scores"], "scores_overall.tsv"))[1]]:
            preds = collections.defaultdict(list)
            for line in open(os.path.join(p["fpreds"], tool + ".preds.tsv")):
                x = line.rstrip("\n").split("\t")
                preds[x[0]].append((int(x[1]), int(x[2]), x[3]))
            acc = {t: [0, 0, 0, 0] for t in TYPES}  # n_truth, n_pred, tp_stop, tp_exact
            negt = collections.Counter()  # predictions on negatives, typed by their own ends
            for tid in negs:
                for q in preds.get(tid, []):
                    negt[orf_type(seqs[tid], *q)] += 1
            out["neg_pred_types"].append([level, tool, sum(negt.values())] + [negt[t] for t in TYPES])
            for tid, r in tt.items():
                a = acc[r["type"]]
                a[0] += 1
                pl = preds.get(tid, [])
                a[1] += len(pl)
                t5, t3 = ends(int(r["cds_start"]), int(r["cds_end"]), r["strand"])
                best = None
                for q in pl:
                    if q[2] != r["strand"]:
                        continue
                    q5, q3 = ends(*q)
                    if abs(q3 - t3) <= TOL:
                        lvl = "exact" if abs(q5 - t5) <= TOL else "stop"
                        if best != "exact":
                            best = lvl
                a[2] += best is not None
                a[3] += best == "exact"
            for t in TYPES:
                n, npred, ts, te = acc[t]
                row = [level, tool, t, n, npred, ts, te]
                for tp in (ts, te):
                    sens = tp / n if n else 0.0
                    prec = tp / npred if npred else 0.0
                    f1 = 2 * sens * prec / (sens + prec) if sens + prec else 0.0
                    row += [f"{sens:.4f}", f"{prec:.4f}", f"{f1:.4f}"]
                # with negatives: predictions on negatives typed as t also count in the denominator
                row += [negt[t]]
                for tp in (ts, te):
                    sens = tp / n if n else 0.0
                    den = npred + negt[t]
                    prec = tp / den if den else 0.0
                    f1 = 2 * sens * prec / (sens + prec) if sens + prec else 0.0
                    row += [f"{prec:.4f}", f"{f1:.4f}"]
                out["scores_by_type"].append(row)
        # E29
        if os.path.exists(p["e29"]):
            _, e = read_tsv(p["e29"])
            n5 = sum(1 for r in e if r["status"] == "5prime_partial")
            out["e29_five_prime"].append([level, len(e), n5, f"{n5 / len(e):.4f}" if e else "NA"])

    heads = {
        "levels": ["level", "fraction", "pairs", "observed_fraction", "trinity_transcripts"],
        "label_counts": ["level", "label", "reason", "n"],
        "truth_counts": ["level"] + TYPES + ["truth_total", "negatives", "ambiguous", "short", "unresolved", "transcripts"]
                        + [f"share_{t}" for t in TYPES],
        "scores_by_type": ["level", "tool", "type", "n_truth", "n_pred_on_type", "tp_stop", "tp_exact",
                           "sensitivity_stop", "precision_stop", "f1_stop", "sensitivity_exact", "precision_exact", "f1_exact",
                           "n_pred_neg_typed", "precision_stop_wneg", "f1_stop_wneg", "precision_exact_wneg", "f1_exact_wneg"],
        "neg_pred_types": ["level", "tool", "n_pred_on_negatives"] + TYPES,
        "e29_five_prime": ["level", "n_with_mane_stop", "n_5prime_partial", "frac_5prime_partial"],
    }
    for name in ("scores_overall", "scores_by_expr", "negatives", "excluded_preds"):
        heads[name] = ["level"] + hdrs[name]
    for name, rows in out.items():
        if name == "levels" and full_pairs:
            rows[0][2] = full_pairs
        with open(os.path.join(outdir, name + ".tsv"), "w") as fh:
            fh.write("\t".join(heads[name]) + "\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    sys.exit(main())
