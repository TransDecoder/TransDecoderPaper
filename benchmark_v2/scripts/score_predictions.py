#!/usr/bin/env python3
"""Score standardized CDS predictions against a truth set.

Predictions: TSV with columns transcript_id, start, end, strand (1-based, inclusive,
transcript coordinates, CDS including the stop codon when present); extra columns
are ignored; a header line starting with 'transcript_id' is optional.

Matching follows the 2018 benchmark (score_TP_FP_FN.py): a prediction matches a truth
CDS on the same transcript and strand when its 3' (stop) end lies within 2 nt of the
truth 3' end ('stop' level); it is an exact match when the 5' (start) end is also
within 2 nt ('exact' level). Unlike 2018, the ends are taken strand-aware and a truth
CDS is credited at most once (extra matching predictions count as false positives).

Writes a one-row TSV of metrics to stdout and, with --detail, per-prediction calls.
Precision counts all predictions, including those on random negative sequences
('rand_' prefix); fp_random is reported separately. 'rand_'-prefixed transcripts never
have truth rows.
"""
import argparse
import collections
import sys

TOL = 2


def ends(start, end, strand):
    """Return (five_prime, three_prime) coordinates."""
    return (start, end) if strand == "+" else (end, start)


def read_truth(path):
    truth = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            truth[r["transcript_id"]] = (int(r["cds_start"]), int(r["cds_end"]), r["strand"], r.get("type", ""))
    return truth


def read_preds(path):
    preds = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "transcript_id\t")):
                continue
            x = line.rstrip("\n").split("\t")
            preds[x[0]].append((int(x[1]), int(x[2]), x[3]))
    return preds


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--preds", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--detail")
    args = ap.parse_args()

    truth = read_truth(args.truth)
    preds = read_preds(args.preds)
    c = collections.Counter()
    detail = open(args.detail, "w") if args.detail else None
    if detail:
        detail.write("call\ttranscript_id\tpred_start\tpred_end\tpred_strand\ttruth_start\ttruth_end\ttruth_strand\n")

    for tid, plist in preds.items():
        t = truth.get(tid)
        is_rand = tid.startswith("rand_")
        # best match first: exact, then stop-level
        calls = []
        for p in plist:
            call = "FP"
            if t and p[2] == t[2]:
                p5, p3 = ends(*p)
                t5, t3 = ends(*t[:3])
                if abs(p3 - t3) <= TOL:
                    call = "TP_exact" if abs(p5 - t5) <= TOL else "TP_stop"
            calls.append([call, p])
        credited = False
        for rank in ("TP_exact", "TP_stop"):
            for cp in calls:
                if cp[0] == rank and not credited:
                    credited = True
                elif cp[0] == rank:
                    cp[0] = "FP_dup"
        for call, p in calls:
            c["n_pred"] += 1
            c["n_pred_random" if is_rand else "n_pred_real"] += 1
            c[call] += 1
            if call.startswith("FP"):
                c["FP"] += 1
                if is_rand:
                    c["FP_random"] += 1
                elif t and p[2] != t[2]:
                    c["FP_wrong_strand"] += 1
            if detail:
                tt = t[:3] if t else (".", ".", ".")
                detail.write("\t".join(map(str, (call, tid) + p + tt)) + "\n")
        if t and not credited:
            c["FN"] += 1
    for tid, t in truth.items():
        if tid not in preds:
            c["FN"] += 1
            if detail:
                detail.write("\t".join(map(str, ("FN", tid, ".", ".", ".") + t[:3])) + "\n")

    n_truth = len(truth)
    tp_stop = c["TP_exact"] + c["TP_stop"]
    tp_exact = c["TP_exact"]
    out = collections.OrderedDict(label=args.label, n_truth=n_truth, n_pred=c["n_pred"],
                                  n_pred_real=c["n_pred_real"], n_pred_random=c["n_pred_random"],
                                  tp_stop=tp_stop, tp_exact=tp_exact, fp=c["FP"],
                                  fp_random=c["FP_random"], fp_wrong_strand=c["FP_wrong_strand"],
                                  fp_dup=c["FP_dup"], fn=c["FN"])
    for lvl, tp in (("stop", tp_stop), ("exact", tp_exact)):
        sens = tp / n_truth if n_truth else 0.0
        prec = tp / c["n_pred"] if c["n_pred"] else 0.0
        f1 = 2 * sens * prec / (sens + prec) if sens + prec else 0.0
        out[f"sensitivity_{lvl}"] = f"{sens:.4f}"
        out[f"precision_{lvl}"] = f"{prec:.4f}"
        out[f"f1_{lvl}"] = f"{f1:.4f}"
    print("\t".join(out))
    print("\t".join(str(v) for v in out.values()))


if __name__ == "__main__":
    sys.exit(main())
