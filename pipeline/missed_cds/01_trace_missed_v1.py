#!/usr/bin/env python3
"""E34: trace every truth CDS of a simulated Trinity assembly through a finished PyTD run.

usage: PYTHONPATH=$P/work/trees/PY_hstart python 01_trace_missed_v1.py <manifest.tsv> <dataset> <big_out_dir>

Reads only the run's own files (input.fa, longest_orfs.{gff3,cds}, longest_orfs.cds.scores, decoy_qvalues.tsv,
hexamer.scores, top_training_orfs.cds, start_refinement.kozak.tsv, best_candidates.gff3, preds.tsv, predict.log)
and, for homology runs, the BLASTP / Pfam files that run was given. Nothing in the run directory is written.
PyTD code (a5b24c3, unchanged) is imported for the frame contrast, Kozak score and start refinement.

Outputs in <big_out_dir>/<dataset>/:
  per_cds.tsv.gz        every scored truth CDS: found / missed, miss class (a-e), candidate scores, winner
  start5open.tsv.gz     truth 3' partials and complete CDSs whose candidate ORF is 5'-open: start decision
  counterfactuals.tsv   selection variants re-run offline (q cutoffs, top candidate on empty transcripts)
  checks.txt            offline selection / refinement reproduce the run (must be 0 differences)

Miss classes (a CDS is missed when no prediction has its stop within 2 nt on its strand):
  d  a prediction overlaps the CDS in its frame but ends at another stop (checked first)
  a  no LongOrfs ORF has the true stop (reason: below_min_length / internal_below_min / other)
  b  the candidate (the LongOrfs ORF with the true stop) fails acceptance: not homology, shorter than the
     long-ORF length, and not (in-frame score > 0, frame contrast > 0, q <= 0.05); sub-reason s0<=0,
     contrast<=0 or q>0.05
  c  the candidate is accepted but another ORF on the transcript ranks above it (single best ORF per
     transcript: homology count, then in-frame score, then length)
  e  anything else
"""
import collections
import gzip
import math
import os
import random
import re
import sys

import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
sys.path.insert(0, f"{P}/TransDecoderPaper/benchmark_v2/scripts")
from score_predictions import TOL, ends  # noqa: E402

from pytransdecoder.predict import TransDecoderPredict, SENSITIVE_DECOY_Q, DEFAULT_DECOY_Q  # noqa: E402
from pytransdecoder.core import kozak as kozak_mod  # noqa: E402

STOPS = {"TAA", "TAG", "TGA"}
EXPR = [(0, 1, "<1"), (1, 5, "1-5"), (5, 20, "5-20"), (20, 100, "20-100"), (100, float("inf"), ">100")]
TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]


def expr_bin(tpm):
    if tpm != tpm:
        return "NA"
    return next(lab for lo, hi, lab in EXPR if lo <= tpm < hi)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def gc3(seq):
    c = seq[2::3]
    return (c.count("G") + c.count("C")) / len(c) if c else float("nan")


def read_mrna(path):
    """ORF id -> dict(tid, start, end, strand, line) from mRNA lines of a TransDecoder gff3."""
    d = {}
    with open(path) as fh:
        for line in fh:
            x = line.rstrip("\n").split("\t")
            if len(x) == 9 and x[2] == "mRNA":
                oid = re.search(r"ID=([^;]+)", x[8]).group(1)
                d[oid] = dict(tid=x[0], start=int(x[3]), end=int(x[4]), strand=x[6], line=line)
    return d


def read_preds(path):
    preds = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(("transcript_id\t", "#")) or not line.strip():
                continue
            x = line.rstrip("\n").split("\t")
            preds[x[0].split()[0]].append((int(x[1]), int(x[2]), x[3]))
    return preds


def main():
    man_path, ds, big = sys.argv[1], sys.argv[2], sys.argv[3]
    man = pd.read_csv(man_path, sep="\t").set_index("dataset")
    m = man.loc[ds]
    run = m.run_dir
    W = f"{run}/input.fa.transdecoder_dir"
    out = f"{big}/{ds}"
    os.makedirs(out, exist_ok=True)
    checks = open(f"{out}/checks.txt", "w")

    def log(*a):
        print(*a, file=checks, flush=True)
        print(ds, *a, file=sys.stderr, flush=True)

    # ---------- inputs ----------
    tx = read_fasta(f"{run}/input.fa")
    lab = m.labels_prefix
    excl = {l.split("\t")[0].strip() for l in open(f"{lab}.excluded.txt") if l.strip()}
    negs = {l.split("\t")[0].strip() for l in open(f"{lab}.negatives.txt") if l.strip()}
    truth = pd.read_csv(f"{lab}.truth.tsv", sep="\t")
    truth = truth[~truth.transcript_id.isin(excl)].set_index("transcript_id")
    log(f"truth {len(truth)} negatives {len(negs)} excluded {len(excl)}")

    lo = read_mrna(f"{W}/longest_orfs.gff3")
    cds = read_fasta(f"{W}/longest_orfs.cds")
    lo_type = {}
    with open(f"{W}/longest_orfs.cds") as fh:
        for line in fh:
            if line.startswith(">"):
                x = line[1:].split()
                lo_type[x[0]] = x[1].split(":", 1)[1]
    sc = {}
    with open(f"{W}/longest_orfs.cds.scores") as fh:
        for line in fh:
            if not line.startswith("#"):
                x = line.rstrip("\n").split("\t")
                sc[x[0]] = [float(v) for v in x[3:9]]
    qv, strat = {}, {}
    with open(f"{W}/decoy_qvalues.tsv") as fh:
        for line in fh:
            if not line.startswith("#"):
                x = line.rstrip("\n").split("\t")
                qv[x[0]] = float(x[3])
                strat[x[0]] = x[1]
    auto = None
    for line in open(f"{run}/predict.log"):
        mm = re.search(r"min auto-accept length = (\S+) nt", line)
        if mm:
            auto = float(mm.group(1))
    params = eval(open(f"{W}/longorfs.params.json").read().replace("true", "True").replace("false", "False"))
    min_aa = params["min_protein_length"]
    log(f"auto-accept {auto} nt; LongOrfs min_protein_length {min_aa}; {len(lo)} LongOrfs ORFs")

    blast_ids, pfam_ids = set(), set()
    homology = m.blast != "-"
    td = TransDecoderPredict(transcripts_file=f"{run}/input.fa", output_dir=f"{out}/td",
                             genetic_code="Universal",
                             retain_blastp_hits=(m.blast if homology else None),
                             retain_pfam_hits=(m.pfam if homology else None))
    os.makedirs(td.workdir, exist_ok=True)
    if homology:
        blast_ids = set(td._parse_blastp_hits())
        pfam_ids = set(td._parse_pfam_hits())
        log(f"homology: {len(blast_ids)} ORFs with BLASTP hit, {len(pfam_ids)} with Pfam")

    by_tid = collections.defaultdict(list)
    for oid, o in lo.items():
        s = sc.get(oid, [0.0] * 6)
        hom = (oid in blast_ids) + (oid in pfam_ids)
        L = o["end"] - o["start"] + 1
        o.update(id=oid, s0=s[0], contrast=s[0] - max(s[1], s[2]), length=L, hom=hom, q=qv.get(oid, 1.0),
                 key=(-hom, -s[0], -L, min(o["start"], o["end"]), max(o["start"], o["end"]), o["strand"]))
        by_tid[o["tid"]].append(o)

    def accepted(o, qthr):
        return o["hom"] > 0 or o["length"] >= auto or (o["s0"] > 0 and o["contrast"] > 0 and o["q"] <= qthr)

    def select(qthr, fill=None):
        """single best ORF per transcript; fill='eligible'/'any': on transcripts with no accepted ORF take the
        top-ranked ORF among eligible (s0>0, contrast>0) / all ORFs."""
        sel = {}
        for tid, orfs in by_tid.items():
            c = [o for o in orfs if accepted(o, qthr)]
            if not c and fill == "eligible":
                c = [o for o in orfs if o["s0"] > 0 and o["contrast"] > 0]
            elif not c and fill == "any":
                c = list(orfs)
            if c:
                sel[tid] = min(c, key=lambda o: o["key"])
        return sel

    # ---------- check offline selection against the run ----------
    sel0 = select(DEFAULT_DECOY_Q)
    run_best = read_mrna(f"{W}/best_candidates.gff3")
    a, b = {o["id"] for o in sel0.values()}, set(run_best)
    log(f"CHECK selection q0.05 offline vs best_candidates.gff3: offline {len(a)} run {len(b)} "
        f"only_offline {len(a - b)} only_run {len(b - a)}")

    # ---------- start refinement with PyTD's own code ----------
    hexf = f"{W}/hexamer.scores"
    trainf = f"{W}/__checkpoints_TDpredict/top_training_orfs.cds"
    hexs = td._load_hexamer_scores(hexf)
    cc = td._coding_contrast(trainf, hexf)
    kz = kozak_mod.read_model(f"{W}/start_refinement.kozak.tsv")
    hom_ev = td._homology_start_evidence() if homology else None
    log(f"median training contrast per base {cc:.4f}")

    def refine(sel, tag):
        src, dst = f"{out}/td/{tag}.gff3", f"{out}/td/{tag}.revised.gff3"
        with open(src, "w") as fh:
            fh.write("##gff-version 3\n")
            for tid in sorted(sel):
                fh.write(sel[tid]["line"])
        td._revise_start_sites(src, dst, tx, hexs, cc, kz, hom_ev)
        r = read_mrna(dst)
        os.replace(src, f"{out}/td/{tag}.selected.gff3")
        return {o["tid"]: (o["start"], o["end"], o["strand"], oid) for oid, o in r.items()}

    tdict = {tid: (int(r.cds_start), int(r.cds_end), r.strand) for tid, r in truth.iterrows()}

    def score(preds):
        """preds: tid -> (start, end, strand, oid). Scored set = truth + negatives (excluded removed)."""
        tp_s = tp_e = n_pred = n_negp = 0
        hit = {}
        for tid, p in preds.items():
            if tid in excl:
                continue
            n_pred += 1
            if tid in negs:
                n_negp += 1
                continue
            if tid in tdict:
                tcs, tce, tst = tdict[tid]
                t5, t3 = ends(tcs, tce, tst)
                if p[2] == tst:
                    p5, p3 = ends(p[0], p[1], p[2])
                    if abs(p3 - t3) <= TOL:
                        tp_s += 1
                        ex = abs(p5 - t5) <= TOL
                        tp_e += ex
                        hit[tid] = "exact" if ex else "stop"
        n = len(truth)
        return dict(n_truth=n, n_pred=n_pred, n_pred_neg=n_negp, n_neg=len(negs), tp_stop=tp_s, tp_exact=tp_e,
                    f1_stop=2 * tp_s / (n + n_pred), f1_exact=2 * tp_e / (n + n_pred)), hit

    base_preds = refine(sel0, "q0.05")
    runp = read_preds(f"{run}/preds.tsv")
    runp1 = {t: v[0] for t, v in runp.items()}
    diff = sum(1 for t in set(runp1) | set(base_preds)
               if runp1.get(t) != (base_preds[t][:3] if t in base_preds else None))
    log(f"CHECK refined offline q0.05 vs run preds.tsv: run {sum(map(len, runp.values()))} offline {len(base_preds)}"
        f" transcripts differing {diff}")
    base_score, base_hit = score(base_preds)
    log("base score", base_score)

    # ---------- per-CDS trace ----------
    train = read_fasta(trainf)
    tr_gc3 = sorted(gc3(s) for s in train.values())
    q05, q95 = tr_gc3[int(0.05 * len(tr_gc3))], tr_gc3[int(0.95 * len(tr_gc3)) - 1]
    log(f"training GC3 median {tr_gc3[len(tr_gc3) // 2]:.3f} 5-95% {q05:.3f}-{q95:.3f} (n={len(tr_gc3)})")
    rng = random.Random(f"E34:{ds}")

    rows, s5rows = [], []
    for tid, t in truth.iterrows():
        cs, ce, st, typ = int(t.cds_start), int(t.cds_end), t.strand, t.type
        t5, t3 = ends(cs, ce, st)
        seq = tx[tid]
        tseq = seq[cs - 1:ce]
        r = dict(dataset=ds, transcript_id=tid, type=typ, n_codons=int(t.n_codons), cds_nt=ce - cs + 1,
                 tx_len=len(seq), tpm=float(t.tpm), expr_bin=expr_bin(float(t.tpm)),
                 gc3=round(gc3(tseq), 4), status=base_hit.get(tid, "missed"),
                 true_internal_stop=any(tseq[i:i + 3] in STOPS for i in range(0, len(tseq) - 3, 3)))
        # candidate: LongOrfs ORF with the true stop
        cand = None
        for o in by_tid.get(tid, []):
            if o["strand"] == st and abs(ends(o["start"], o["end"], o["strand"])[1] - t3) <= TOL:
                cand = o
        win = sel0.get(tid)
        if cand:
            r.update(cand_id=cand["id"], cand_type=lo_type[cand["id"]], cand_len=cand["length"],
                     cand_s0=cand["s0"], cand_contrast=round(cand["contrast"], 2),
                     cand_contrast_pb=round(cand["contrast"] / cand["length"], 5),
                     cand_s0_pb=round(cand["s0"] / cand["length"], 5),
                     cand_q=qv.get(cand["id"], float("nan")), cand_stratum=strat.get(cand["id"], "."),
                     cand_eligible=cand["s0"] > 0 and cand["contrast"] > 0, cand_hom=cand["hom"],
                     cand_blast=cand["id"] in blast_ids, cand_pfam=cand["id"] in pfam_ids,
                     cand_auto=cand["length"] >= auto, cand_accepted=accepted(cand, DEFAULT_DECOY_Q),
                     cand_accept_q02=accepted(cand, SENSITIVE_DECOY_Q), cand_accept_q05=accepted(cand, 0.5),
                     cand_rank=sorted(by_tid[tid], key=lambda o: o["key"]).index(cand) + 1,
                     n_orfs_tx=len(by_tid[tid]))
        if win:
            r.update(win_id=win["id"], win_len=win["length"], win_s0=win["s0"],
                     win_contrast=round(win["contrast"], 2), win_hom=win["hom"], win_q=win["q"])
        cls = sub = ""
        if r["status"] == "missed":
            p = base_preds.get(tid)
            same_frame = False
            if p and p[2] == st:
                p5, p3 = ends(p[0], p[1], p[2])
                lo_, hi_ = min(p[0], p[1]), max(p[0], p[1])
                if (p5 - t5) % 3 == 0 and lo_ <= max(cs, ce) and hi_ >= min(cs, ce):
                    same_frame = True
            if same_frame:
                cls, sub = "d", "same_frame_other_stop"
            elif cand is None:
                cls = "a"
                # the ORF LongOrfs would report for this stop: back to the codon after the previous in-frame
                # stop; there the first ATG, or the open frame start when no stop lies upstream
                j = t5 - 1
                while j - 3 >= 0 and seq[j - 3:j] not in STOPS:
                    j -= 3
                if j < 3:
                    start = j
                else:
                    start = next((x for x in range(j, t5, 3) if seq[x:x + 3] == "ATG"), t5 - 1)
                orf_nt = t3 - start
                r["a_orf_nt"] = orf_nt
                has_stop = typ in ("complete", "5prime_partial")
                minnt = 3 * min_aa if has_stop else 3 * min_aa - 3
                if orf_nt < minnt:
                    sub = "below_min_length" if has_stop else "internal_or_3open_below_min"
                else:
                    sub = "other"
            elif not accepted(cand, DEFAULT_DECOY_Q):
                cls = "b"
                sub = ("s0<=0" if cand["s0"] <= 0 else "contrast<=0" if cand["contrast"] <= 0 else "q>0.05")
            elif win is not None and win["id"] != cand["id"]:
                cls = "c"
                w = win
                same_strand = w["strand"] == st
                wl, wh = min(w["start"], w["end"]), max(w["start"], w["end"])
                cl, ch = min(cand["start"], cand["end"]), max(cand["start"], cand["end"])
                ov = max(0, min(wh, ch) - max(wl, cl) + 1)
                if not same_strand:
                    sub = "antisense"
                elif w["hom"] > cand["hom"]:
                    sub = "homology_ranked"
                else:
                    sub = "other_frame_overlap" if ov > 0 else "nonoverlapping_orf"
                r.update(win_overlap_nt=ov, win_longer=w["length"] > cand["length"],
                         win_higher_contrast=w["contrast"] > cand["contrast"],
                         win_same_frame=same_strand and (ends(w["start"], w["end"], w["strand"])[1] - t3) % 3 == 0)
            else:
                cls, sub = "e", "selected_but_no_pred" if win and win["id"] == cand["id"] else "other"
        r.update(miss_class=cls, miss_sub=sub)
        # true-region contrast (missed + a sample of found)
        if r["status"] == "missed" or rng.random() < 0.15:
            r["true_contrast_pb"] = round(td._frame_contrast(tseq, hexs) / len(tseq), 5)
            r["true_s0_pb"] = round(td._score_cds_via_markov(tseq, 5, hexs) / len(tseq), 5)
        rows.append(r)

        # ---------- part 2: 5'-open candidates of 3' partials and complete CDSs ----------
        if cand and typ in ("3prime_partial", "complete") and lo_type[cand["id"]] in ("5prime_partial", "internal"):
            o5 = ends(cand["start"], cand["end"], cand["strand"])[0]
            orf = seq[cand["start"] - 1:cand["end"]]
            atg = next((i for i in range(3, len(orf) - 3, 3) if orf[i:i + 3] == "ATG"), None)
            true_off = t5 - o5
            s = dict(dataset=ds, transcript_id=tid, type=typ, cand_id=cand["id"], cand_type=lo_type[cand["id"]],
                     u_true=true_off, atg_first=atg, tx_len=len(seq), cand_len=cand["length"])
            if atg is not None:
                coding = td._frame_contrast(orf[:atg], hexs) / cc
                K = kozak_mod.score(kz, seq, cand["start"] - 1 + atg)
                utr = td.start_contrast_alpha * atg + td.start_kozak_weight * K
                s.update(coding_u=round(coding, 2), kozak=round(K, 3), utr_u=round(utr, 2),
                         ab_initio=("trim" if coding < utr else "keep"),
                         true_is_first_atg=abs(true_off - atg) <= TOL,
                         true_atg_rank=sum(1 for i in range(3, true_off + 1, 3) if orf[i:i + 3] == "ATG"),
                         kozak_true=round(kozak_mod.score(kz, seq, t5 - 1), 3) if seq[t5 - 1:t5 + 2] == "ATG" else None)
            p = base_preds.get(tid)
            outcome = "missed"
            if r["status"] != "missed":
                p5 = ends(p[0], p[1], p[2])[0]
                if abs(p5 - t5) <= TOL:
                    outcome = "exact"
                elif abs(p5 - o5) <= TOL:
                    outcome = "kept_open"
                else:
                    outcome = "trimmed_other_atg"
            s["outcome"] = outcome
            s["miss_class"] = cls
            s5rows.append(s)

    D = pd.DataFrame(rows)
    D.to_csv(f"{out}/per_cds.tsv.gz", sep="\t", index=False)
    S = pd.DataFrame(s5rows)
    S.to_csv(f"{out}/start5open.tsv.gz", sep="\t", index=False)
    log("classes", D[D.status == "missed"].groupby(["type", "miss_class"]).size().to_dict())

    # homology start-assist decisions (B): attach source from the run's table
    if homology:
        H = pd.read_csv(f"{W}/start_refinement.homology.tsv", sep="\t", dtype=str)
        H[["orf_id", "rules_fired", "decision", "source", "ab_initio_decision"]].to_csv(
            f"{out}/homology_decisions.tsv.gz", sep="\t", index=False)

    # ---------- counterfactuals ----------
    missed0 = set(D.transcript_id[D.status == "missed"])
    cf = []

    def add(name, preds):
        s, hit = score(preds)
        rec = set(hit) & missed0
        lost = {t for t in base_hit if t not in hit}
        lost_exact = {t for t in base_hit if base_hit[t] == "exact" and hit.get(t) != "exact"}
        s.update(variant=name, recovered_stop=len(rec), recovered_exact=sum(hit[t] == "exact" for t in rec),
                 lost_stop=len(lost), lost_exact=len(lost_exact),
                 neg_tx_with_pred=sum(1 for t in preds if t in negs))
        for typ in TYPES:
            ids = set(D.transcript_id[(D.type == typ) & (D.status == "missed")])
            s[f"recovered_stop_{typ}"] = len(set(hit) & ids)
        cf.append(s)
        log(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()})

    add("q0.05 (release)", base_preds)
    for qthr, name in ((0.1, "q0.1"), (SENSITIVE_DECOY_Q, "q0.2 (--sensitive)"), (0.5, "q0.5"), (1.0, "q1.0")):
        add(name, refine(select(qthr), name.split()[0]))
    add("top eligible ORF on empty transcripts", refine(select(DEFAULT_DECOY_Q, "eligible"), "fill_eligible"))
    add("top ORF on empty transcripts", refine(select(DEFAULT_DECOY_Q, "any"), "fill_any"))
    pd.DataFrame(cf).to_csv(f"{out}/counterfactuals.tsv", sep="\t", index=False)
    log("DONE")


if __name__ == "__main__":
    main()
