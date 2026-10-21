#!/usr/bin/env python
"""Tabulate out/*.signal.json and out/*.sim.json into markdown tables (stdout)."""
import json, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pwmlib import ORGS
OUT = "/broad/hptmp/bhaas/TRANSDECODER_PAPER/work/analysis/pwm_audit/out"
what = sys.argv[1] if len(sys.argv) > 1 else "signal"
variants = ["complete", "complete_fullutr"]


def load(kind):
    d = {}
    for o in ORGS:
        for v in variants:
            p = f"{OUT}/{o}.{v}.{kind}.json"
            if os.path.exists(p):
                d[(o, v)] = json.load(open(p))
    return d


if what == "signal":
    S = load("signal")
    print("\n### Q1. Training positives (what the self-trained PWM learns from)\n")
    print("| org | variant | train ORFs | no ATG (skipped) | positives | = true start | true start, UTR added | upstream in-frame ATG (N-ext) | other ORF | enhanced kept |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for (o, v), r in S.items():
        c = r["q1_counts"]
        ts = c.get("pos_true_start", 0)
        tsa = c.get("pos_true_start|utr5_added", 0)
        upn = sum(val for k, val in c.items() if k.startswith("pos_upstream_inframe"))
        oth = sum(val for k, val in c.items() if k.startswith("pos_") and "true_start" not in k and "upstream_inframe" not in k)
        print(f"| {o} | {v} | {c['train_orfs']} | {c.get('no_ATG_skipped',0)} | {c['positives']} | {ts} | {tsa} | {upn} | {oth} | {c['enhanced_n']} |")
    print("\n### Q1. Training negatives composition (ATGs downstream of the training-ORF start)\n")
    print("| org | variant | n neg | in CDS, in-frame | in CDS, out-of-frame | 3' UTR |")
    print("|---|---|---|---|---|---|")
    for (o, v), r in S.items():
        n = r["q1_negatives_composition"]
        tot = sum(n.values())
        print(f"| {o} | {v} | {tot} | {n.get('in_CDS_inframe',0)} ({n.get('in_CDS_inframe',0)/tot:.0%}) | {n.get('in_CDS_oof',0)} ({n.get('in_CDS_oof',0)/tot:.0%}) | {n.get('3UTR',0)} ({n.get('3UTR',0)/tot:.0%}) |")
    print("\n### Q2. Selected range/threshold and what passes\n")
    print("| org | variant | enh range (up,down) | enh thr/base | reported AUC (enh) | best F1 | grid step | raw positives passing thr | train negatives passing thr | raw range | raw AUC |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for (o, v), r in S.items():
        q = r["q2"]
        print(f"| {o} | {v} | {q['enh_range']} | {q['enh_threshold']:.3f} | {q['enh_auc_reported']:.3f} | {q['enh_bestF1']:.2f} | {q['threshold_grid_step']:.3f} | {q['frac_raw_positives_pass_threshold']:.0%} | {q['frac_training_negatives_pass_threshold']:.1%} | {q['raw_range']} | {q['raw_auc_reported']:.3f} |")
    print("\n### Q3. AUC on truth-defined contrasts (positives = annotated true starts; negatives as named)\n")
    print("Self-trained = the PWM files from the pytd_dev_pwm run; truth-trained = 2-fold CV on true starts, full 20+3+10 window unless noted.\n")
    print("| org | variant | n pos | contrast | self enh, chosen range | self enh, full | self raw, full | self raw, Kozak(6,3) | truth (neg=TD down_any) full | truth (neg=in-frame down) full | truth (neg=mixed up/oof/down) full | truth mixed, Kozak(6,3) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for (o, v), r in S.items():
        for k in ("in_frame_down", "in_frame_up", "utr_oof", "down_any"):
            g = lambda d: d.get(k, float('nan'))
            print(f"| {o} | {v} | {r['q3_set_sizes']['pos']} | {k} (n={r['q3_set_sizes'].get(k,0)}) | {g(r['q3_auc_self_enhanced_chosen_range']):.3f} | {g(r['q3_auc_self_enhanced_full']):.3f} | {g(r['q3_auc_self_raw_full']):.3f} | {g(r['q3_auc_self_raw_kozak_6_3']):.3f} | {g(r['q3_auc_truth_TDneg_down_any']['20,10']):.3f} | {g(r['q3_auc_truth_neg_inframe_down']['20,10']):.3f} | {g(r['q3_auc_truth_neg_mixed_up_oof']['20,10']):.3f} | {g(r['q3_auc_truth_neg_mixed_up_oof']['6,3']):.3f} |")
    print("\n### Q3. Kozak positions in truth-trained PWM (annotated UTRs): base probabilities at -3 and +4, and per-position KL (pos vs in-frame-down ATGs)\n")
    print("| org | variant | -3 A | -3 G | -3 A+G (in-frame-down ATG: A+G) | +4 G (in-frame-down: G) | KL -3 | KL +4 | max KL at other pos (pos idx) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for (o, v), r in S.items():
        p = r["q3_truth_pos_pwm_key_positions"]; d = r["q3_inframe_down_pwm_key_positions"]
        kl = r["q3_kl_by_pos"]["pos_vs_inframe_down"]
        others = [(x, i) for i, x in enumerate(kl) if i not in (17, 23, 20, 21, 22)]
        mx = max(others)
        print(f"| {o} | {v} | {p['-3']['A']:.2f} | {p['-3']['G']:.2f} | {p['-3']['A']+p['-3']['G']:.2f} ({d['-3']['A']+d['-3']['G']:.2f}) | {p['+4']['G']:.2f} ({d['+4']['G']:.2f}) | {kl[17]:.3f} | {kl[23]:.3f} | {mx[0]:.3f} ({mx[1]-20:+d}) |")
    print("\n### Q4. Truth: where is the true start among in-frame ATG candidates? (annotated 5' UTRs; closed = in-frame stop upstream in the transcript)\n")
    print("| org | variant | n | closed, true start = first ATG after stop | closed, later ATG | open, true start = first in-frame ATG | open, later ATG | later total |")
    print("|---|---|---|---|---|---|---|---|")
    for (o, v), r in S.items():
        q = r["q4_truth_candidate_structure_annotated_utr"]
        n = q["n"]
        later = q.get("closed_later", 0) + q.get("open_later", 0)
        print(f"| {o} | {v} | {n} | {q.get('closed_first',0)} ({q.get('closed_first',0)/n:.1%}) | {q.get('closed_later',0)} ({q.get('closed_later',0)/n:.1%}) | {q.get('open_first',0)} ({q.get('open_first',0)/n:.1%}) | {q.get('open_later',0)} ({q.get('open_later',0)/n:.1%}) | {later} ({later/n:.1%}) |")
    print("\n### Q3 extra. Truth-trained AUC by range (neg = mixed in-frame-up / UTR out-of-frame / in-frame-down), contrast = in_frame_down / in_frame_up / utr_oof\n")
    print("| org | variant | 20,10 | 10,5 | 6,3 | 3,3 | 3,1 | 6,1 | 10,1 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for (o, v), r in S.items():
        t = r["q3_auc_truth_neg_mixed_up_oof"]
        cells = []
        for rg in ("20,10", "10,5", "6,3", "3,3", "3,1", "6,1", "10,1"):
            x = t[rg]
            cells.append("/".join(f"{x.get(k, float('nan')):.2f}" for k in ("in_frame_down", "in_frame_up", "utr_oof")))
        print(f"| {o} | {v} | " + " | ".join(cells) + " |")
else:
    S = load("sim")
    rules = ["first", "current_pwm", "self_raw_argmax", "self_raw_margin2", "self_raw_margin3", "self_enh_chosen_argmax",
             "truth_argmax", "truth_margin2", "truth_margin3", "truth_td_argmax", "oracle"]
    for mode in ("open", "closed_ext"):
        print(f"\n### Q5. Simulation, mode = {mode}: exact starts among stop-correct ORFs (count; n = ORFs in this mode)\n")
        print("| org | variant | n | multi-cand | truth in cands | truth is first | " + " | ".join(rules) + " |")
        print("|---|---|---|---|---|---|" + "---|" * len(rules))
        for (o, v), r in S.items():
            m = r[mode]
            print(f"| {o} | {v} | {m['n']} | {m.get('multi_cand',0)} | {m.get('truth_in_cands',0)} | {m.get('truth_is_first',0)} | " + " | ".join(str(m.get(k, 0)) for k in rules) + " |")
