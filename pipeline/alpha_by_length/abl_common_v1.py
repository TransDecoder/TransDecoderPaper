"""E36: shared code. Per-candidate labels and a fast exact rescoring of start-rule variants.

The release rule trims a 5'-open ORF to its first in-frame ATG when alpha*|u| + 17 K > coding(u).
With every other quantity fixed (E35 candidate tables), candidate i is trimmed iff alpha_i > a*_i, where
    a*_i = (coding_i - 17 K_i) / |u_i|          (|u| >= 3 always).
A per-transcript alpha_i = g(L, O, ...) therefore only changes which candidates fall on each side of a*_i.
With the homology start assist, the run's homology call (trim / keep) overrides where it fires (as PyTD a5b24c3).

Scoring (identical to E35 02_evaluate_v1.py / score_predictions.py): stop within 2 nt, exact = start too,
one credit per truth CDS; overall F1 = 2 tp / (n_truth + n_pred); per type f1 = 2 tp / (n + npred of that
truth type); f1w ("with negatives") also counts predictions on negative contigs typed by their own ends.
Only starts change, so n_pred, tp_stop and the truth-typed prediction counts are constant; tp_exact and the
negative-contig type counts are recomputed from per-candidate flags (labels table, 02_labels_v1.py).
"""
import numpy as np
import pandas as pd

TYPES = ["complete", "5prime_partial", "3prime_partial", "internal"]
W = 17.0


def astar(lab):
    return (lab.coding.values - W * lab.kozak.values) / lab.seg.values.astype(float)


def decide(lab, alpha, homology=True):
    """trim decision per candidate for per-candidate alpha (array or scalar)."""
    ab = lab.coding.values < alpha * lab.seg.values.astype(float) + W * lab.kozak.values
    if not homology:
        return ab
    h = lab.hcall.values
    return np.where(h == "trim", True, np.where(h == "keep", False, ab))


class FastScorer:
    """Exact rescoring of one dataset for any decision vector (checked against full substitution)."""

    def __init__(self, lab, const):
        self.c = const
        m = lab.on_excl.values == 0
        self.lab = lab
        self.m = m
        # transcript groups of scored candidates with a truth CDS
        tt = lab.truth_type.values
        has = m & (tt != "-")
        self.has = has
        codes, uniq = pd.factorize(lab.transcript_id.values[has])
        self.grp = codes
        self.ng = len(uniq)
        self.gtype = np.empty(self.ng, dtype=object)
        self.gtype[codes] = tt[has]
        self.other = lab.other_exact.values[has].astype(bool)
        self.ek = lab.e_keep.values[has].astype(bool)
        self.et = lab.e_trim.values[has].astype(bool)
        self.neg = m & (lab.is_neg.values == 1)
        self.negstop = lab.orf_has_stop.values[self.neg].astype(bool)
        d0 = lab.trim_run.values.astype(bool)
        self.tp0_g = self._tp_groups(d0[has])
        self.neg0 = self._negtypes(d0[self.neg])

    def _tp_groups(self, dh):
        ex = self.other | np.where(dh, self.et, self.ek)
        return np.bincount(self.grp, weights=ex, minlength=self.ng) > 0

    def _negtypes(self, dn):
        # keep -> 5prime_partial (stop) / internal; trim -> complete (stop) / 3prime_partial
        s = self.negstop
        return {"complete": int((dn & s).sum()), "3prime_partial": int((dn & ~s).sum()),
                "5prime_partial": int((~dn & s).sum()), "internal": int((~dn & ~s).sum())}

    def tp_by_type(self, d):
        g = self._tp_groups(d[self.has])
        out = {}
        for t in TYPES:
            sel = self.gtype == t
            out[t] = int(g[sel].sum()) - int(self.tp0_g[sel].sum())
        return out  # change vs run

    def score(self, d, full=True):
        c = self.c
        dtp = self.tp_by_type(d)
        tp = c["tp_exact"] + sum(dtp.values())
        row = dict(tp_exact=tp, f1_exact=2 * tp / (c["n_truth"] + c["n_pred"]),
                   n_changed=int((d != self.lab.trim_run.values.astype(bool))[self.m].sum()))
        if not full:
            return row
        dn = self._negtypes(d[self.neg])
        for t in TYPES:
            n, npr = c[f"n_{t}"], c[f"npred_{t}"]
            tpt = c[f"tp_{t}"] + dtp[t]
            negt = c[f"negpred_{t}"] - self.neg0[t] + dn[t]
            row[f"tp_{t}"] = tpt
            row[f"f1_{t}"] = 2 * tpt / (n + npr) if n + npr else np.nan
            sens = tpt / n if n else 0.0
            den = npr + negt
            prec = tpt / den if den else 0.0
            row[f"f1w_{t}"] = 2 * sens * prec / (sens + prec) if sens + prec else (np.nan if not n else 0.0)
            row[f"tpstop_{t}"] = c.get(f"tpstop_{t}", np.nan)
        return row


def load(labels_dir, ds):
    lab = pd.read_csv(f"{labels_dir}/{ds}.labels.tsv.gz", sep="\t", keep_default_na=False,
                      dtype={"coding": float, "kozak": float, "truth_type": str, "hcall": str})
    const = pd.read_csv(f"{labels_dir}/{ds}.const.tsv", sep="\t").iloc[0].to_dict()
    return lab, const


# ---------------------------------------------------------------- alpha forms
L0 = 1000.0   # nt, reference transcript length for the log form
O0 = 1000.0   # nt, reference ORF length
R0 = 0.8      # reference O/L


def alpha_form(form, p, lab):
    """per-candidate alpha for form `form` with parameter dict p; clipped to [0, 5]."""
    L = lab.tx_len.values.astype(float)
    O = lab.orf_len.values.astype(float)
    U = lab.seg.values.astype(float)
    if form == "const":
        a = np.full(len(L), p["a0"])
    elif form == "logL":
        a = p["a0"] + p["b"] * np.log2(L / L0)
    elif form == "logO":
        a = p["a0"] + p["b"] * np.log2(O / O0)
    elif form == "OL":
        a = p["a0"] + p["b"] * (O / L - R0)
    elif form == "uL":
        a = p["a0"] + p["b"] * (U / L)
    elif form == "stepOL":
        a = np.where(O / L >= p["r"], p["ahi"], p["alo"])
    elif form == "stepL":
        a = np.where(L >= p["l"], p["ahi"], p["alo"])
    else:
        raise ValueError(form)
    return np.clip(a, 0.0, 5.0)
