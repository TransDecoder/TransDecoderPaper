#!/usr/bin/env python3
"""E33 step 4: anatomy of the CDSs that the homology start assist (PyTD a5b24c3) still gets wrong.

usage: python step4_02_error_anatomy_v1.py <out_dir>

Runs (assist on, read only):
  sim.human.{d1,d8,d32}, sim.athal.d1: $P/work/start_classifier/hstart_runs/<ds>/on (step-3 candidate
    searches; homology does not change selection there);
  sim.athal.d1.B: the athal depth series condition B (full DIAMOND/Pfam on all LongOrfs, assist on),
    $P/work/start_classifier/athal_depth/B/depth/d1/runs/pytd_v7 (another agent's output; read only).
Every truth CDS (excluded contigs removed) of type 5' partial, 3' partial, internal or complete that is
not exact (stop and start within 2 nt, as score_predictions.py) gets exactly one class:
  selection        no prediction with the true stop; reason no_prediction (none on the transcript),
                   wrong_strand, wrong_stop_or_orf (a prediction on the strand, other stop)
  not_5open        stop right, but the predicted ORF is not a 5'-open candidate with an in-frame ATG
                   (5'-complete ORF with a wrong start, or 5'-open with no ATG to trim to)
  start_elsewhere  stop right, 5'-open candidate, but the true start is neither the open end nor the
                   first ATG (a later ATG, or upstream of the transcript): neither decision is exact
  no_hit           wrong keep/trim decision; the ORF has no DIAMOND hit and no Pfam domain
  hit_no_rule      wrong decision; a hit exists but no rule fired (ab initio rule decided); reason below
  rule_wrong       wrong decision made by a homology rule (R1, R2, R3, R1+R3)
  conflict         wrong decision after an R2 + R3 conflict fell back to the ab initio rule
  other            decision equals the needed one but still not exact (should not occur)
hit_no_rule reasons (a = peptide position of the first ATG; BLAST = best hit, first HSP):
  blast_up_sstart1     qstart <= a-5 but sstart = 1 (R1 requires sstart > 1)
  blast_up_sstart2_3_near?  (cannot happen: R1 would fire)
  blast_ss4_10         qstart > a-5 and 4 <= sstart <= 10 (subject N-terminus just missed, R2 needs <= 3)
  blast_ss11_50        qstart > a-5 and 11 <= sstart <= 50
  blast_ss_gt50        qstart > a-5 and sstart > 50 (alignment starts deep inside the subject)
  blast_nterm_off      sstart <= 3 but the subject N-terminus maps > 5 aa from the ATG (and R1 not met)
  pfam_only_down       no BLAST hit; Pfam envelope starts at or after the ATG
The needed decision: keep open if the true 5' end is within 2 nt of the open ORF's 5' end; trim if it
is within 2 nt of the first in-frame ATG.
Writes anatomy_counts.tsv (dataset x truth type x class x reason x needed decision) and
anatomy_per_cds.tsv.gz.
"""
import collections
import gzip
import os
import sys

import pandas as pd

P = "/broad/hptmp/bhaas/TRANSDECODER_PAPER"
sys.path.insert(0, f"{P}/TransDecoderPaper/benchmark_v2/scripts")
from score_predictions import TOL, ends, read_truth  # noqa: E402

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
man = pd.read_csv(f"{P}/TransDecoderPaper/data/start_classifier/step3_manifest.tsv", sep="\t").set_index("dataset")
RUNS = [(ds, ds, f"{P}/work/start_classifier/hstart_runs/{ds}/on")
        for ds in ("sim.human.d1", "sim.human.d8", "sim.human.d32", "sim.athal.d1")]
RUNS.append(("sim.athal.d1.B", "sim.athal.d1", f"{P}/work/start_classifier/athal_depth/B/depth/d1/runs/pytd_v7"))


def read_cds(path, final):
    """orf_id -> (tid, start, end, strand) from CDS lines."""
    d = {}
    for line in open(path):
        x = line.rstrip("\n").split("\t")
        if len(x) < 9 or x[2] != "CDS":
            continue
        oid = x[8].split("ID=")[1].split(";")[0]
        if oid.startswith("cds."):
            oid = oid[4:]
        d[oid] = (x[0], int(x[3]), int(x[4]), x[6])
    return d


def reason(a, qs, ss, pe):
    if qs is not None:
        if qs <= a - 5 and ss == 1:
            return "blast_up_sstart1"
        if ss <= 3:
            return "blast_nterm_off"
        if ss <= 10:
            return "blast_ss4_10"
        if ss <= 50:
            return "blast_ss11_50"
        return "blast_ss_gt50"
    return "pfam_only_down"


rows = []
for name, ds, run in RUNS:
    m = man.loc[ds]
    excl = {l.split("\t")[0] for l in open(m.excluded)}
    truth = {k: v for k, v in read_truth(m.truth).items() if k not in excl}
    fin = read_cds(f"{run}/input.fa.transdecoder.gff3", True)
    lo = read_cds(f"{run}/input.fa.transdecoder_dir/longest_orfs.gff3", False)
    H = pd.read_csv(f"{run}/input.fa.transdecoder_dir/start_refinement.homology.tsv", sep="\t",
                    dtype=str).set_index("orf_id")
    by_tid = collections.defaultdict(list)
    for oid, (tid, s, e, st) in fin.items():
        if tid not in excl:
            by_tid[tid].append((oid, s, e, st))
    for tid, (ts, te, tst, typ) in truth.items():
        t5, t3 = ends(ts, te, tst)
        best, exact = None, False
        for oid, s, e, st in by_tid.get(tid, []):
            if st != tst:
                continue
            p5, p3 = ends(s, e, st)
            if abs(p3 - t3) <= TOL:
                if abs(p5 - t5) <= TOL:
                    exact = True
                best = best or (oid, p5)
        if exact:
            continue
        r = dict(dataset=name, transcript_id=tid, type=typ, cls=None, reason="", needed="", decision="",
                 rules="", orf_id="")
        if best is None:
            r["cls"] = "selection"
            same = [x for x in by_tid.get(tid, []) if x[3] == tst]
            r["reason"] = "no_prediction" if not by_tid.get(tid) else ("wrong_strand" if not same else "wrong_stop_or_orf")
        else:
            oid, p5 = best
            r["orf_id"] = oid
            if oid not in H.index:
                r["cls"] = "not_5open"
            else:
                h = H.loc[oid]
                o5 = ends(*lo[oid][1:])[0]  # open 5' end (longest ORF)
                sgn = 1 if tst == "+" else -1
                atg5 = o5 + sgn * int(h.atg_offset_nt)
                need = "keep" if abs(t5 - o5) <= TOL else ("trim" if abs(t5 - atg5) <= TOL else None)
                r.update(decision=h.decision, rules=h.rules_fired)
                if need is None:
                    r["cls"] = "start_elsewhere"
                    r["needed"] = "downstream_atg" if sgn * (t5 - atg5) > 0 else "upstream_of_atg"
                else:
                    r["needed"] = need
                    if h.decision == need:
                        r["cls"] = "other"
                    elif h.source == "homology":
                        r["cls"] = "rule_wrong"
                    elif h.source == "ab_initio(conflict)":
                        r["cls"] = "conflict"
                    elif h.blast_qstart == "." and h.pfam_env_from == ".":
                        r["cls"] = "no_hit"
                    else:
                        r["cls"] = "hit_no_rule"
                        a = int(h.atg_aa)
                        qs = None if h.blast_qstart == "." else int(h.blast_qstart)
                        ss = None if h.blast_sstart == "." else int(h.blast_sstart)
                        pe = None if h.pfam_env_from == "." else int(h.pfam_env_from)
                        r["reason"] = reason(a, qs, ss, pe)
        rows.append(r)
    print(name, "done", file=sys.stderr)

D = pd.DataFrame(rows)
D.to_csv(f"{out}/anatomy_per_cds.tsv.gz", sep="\t", index=False)
C = D.groupby(["dataset", "type", "cls", "reason", "needed", "rules"]).size().rename("n").reset_index()
C.to_csv(f"{out}/anatomy_counts.tsv", sep="\t", index=False)
print(D.groupby(["dataset", "type", "cls"]).size().unstack("cls", fill_value=0).to_string())
