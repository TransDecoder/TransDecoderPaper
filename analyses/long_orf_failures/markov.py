"""Vectorised re-implementation of PyTransDecoder's coding model (predict.py at PY_dev 48c870c).

Reproduces, without modifying PyTD:
  - _train_hexamer_model: three-periodic 5th-order Markov chain, pseudocounts (n+1)/(n_prefix+4),
    log-odds against the base_freqs.dat background (3-decimal ratios, as PyTD reads them);
    k-mers never seen in training are absent from hexamer.scores and score 0;
  - _score_cds_via_markov: sum of framed k-mer log-odds, k = min(i,5)+1, stopping before a
    terminal stop codon in frame 0;
  - the six frame scores written to longest_orfs.cds.scores (S1..S3 forward offsets 0/1/2,
    S4..S6 reverse complement offsets 0/1/2).
"""
import numpy as np

_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _LUT[ord(_b)] = _i
_RC = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(_RC)[::-1]


def encode(seq):
    return _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def _n_used(seq, stops):
    """Positions scored: PyTD stops at i == L-3 when that is frame 0 and a stop codon."""
    L = len(seq)
    if L >= 3 and (L - 3) % 3 == 0 and seq[L - 3:] in stops:
        return L - 3
    return L


def _kmer_index(codes, k, n):
    """For positions i in [k-1, n): index of codes[i-k+1..i] (base 4) and validity."""
    m = n - (k - 1)
    if m <= 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=bool)
    idx = np.zeros(m, dtype=np.int64)
    valid = np.ones(m, dtype=bool)
    for j in range(k):
        c = codes[j:j + m]
        valid &= c < 4
        idx = idx * 4 + np.where(c < 4, c, 0)
    return idx, valid


class MarkovModel:
    def __init__(self, tables):
        # tables[k] : float array (3, 4**k), k = 1..6; 0 where the framed k-mer is unknown
        self.tables = tables

    @classmethod
    def from_scores_file(cls, path):
        tables = {k: np.zeros((3, 4 ** k)) for k in range(1, 7)}
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                kmer, frame = p[0].rsplit("-", 1)
                idx = 0
                for b in kmer:
                    idx = idx * 4 + "ACGT".index(b)
                tables[len(kmer)][int(frame), idx] = float(p[3])
        return cls(tables)

    @classmethod
    def train(cls, seqs, bg, stops):
        """seqs: iterable of uppercase CDS strings; bg: {base: freq}; stops: set of stop codons."""
        counts = {k: np.zeros((3, 4 ** k), dtype=np.int64) for k in range(1, 7)}
        frame_counts = np.zeros(3, dtype=np.int64)
        for seq in seqs:
            n = _n_used(seq, stops)
            if n <= 0:
                continue
            codes = encode(seq)
            frame_counts += np.bincount(np.arange(n) % 3, minlength=3)
            for k in range(1, 7):
                idx, valid = _kmer_index(codes, k, n)
                if not len(idx):
                    continue
                frames = (np.arange(k - 1, n) % 3)[valid]
                counts[k] += np.bincount(frames * 4 ** k + idx[valid],
                                         minlength=3 * 4 ** k).reshape(3, 4 ** k)
        logbg = np.log(np.array([bg[b] for b in "ACGT"]))
        tables = {}
        for k in range(1, 7):
            c = counts[k].astype(float)
            if k == 1:
                prefix = np.repeat(frame_counts[[2, 0, 1]][:, None], 4, axis=1).astype(float)
            else:
                pc = counts[k - 1][[2, 0, 1]]           # prefix in frame (f-1) mod 3
                prefix = np.repeat(pc, 4, axis=1).astype(float)   # idx // 4 -> prefix index
            last = np.tile(np.arange(4), 4 ** (k - 1))
            ll = np.log((c + 1) / (prefix + 4)) - logbg[last][None, :]
            tables[k] = np.where(c > 0, ll, 0.0)
        return cls(tables)

    def score(self, seq, stops):
        n = _n_used(seq, stops)
        if len(seq) < 6 or n <= 0:
            return 0.0
        codes = encode(seq)
        total = 0.0
        for i in range(min(5, n)):          # short k-mers at the start
            k = i + 1
            sub = codes[0:i + 1]
            if (sub < 4).all():
                idx = 0
                for c in sub:
                    idx = idx * 4 + int(c)
                total += self.tables[k][i % 3, idx]
        if n > 5:
            idx, valid = _kmer_index(codes, 6, n)
            frames = np.arange(5, n) % 3
            total += float(self.tables[6][frames[valid], idx[valid]].sum())
        return total

    def frame_scores(self, seq, stops):
        rc = revcomp(seq)
        return [self.score(seq, stops), self.score(seq[1:], stops), self.score(seq[2:], stops),
                self.score(rc, stops), self.score(rc[1:], stops), self.score(rc[2:], stops)]

    def contrast(self, seq, stops):
        return self.score(seq, stops) - max(self.score(seq[1:], stops), self.score(seq[2:], stops))


def read_bg(path):
    bg = {}
    for line in open(path):
        p = line.split()
        if len(p) == 3:
            bg[p[0]] = float(p[2])
        elif len(p) == 2:
            bg[p[0]] = float(p[1])
    return bg


def median_contrast(model, seqs, stops):
    """PyTD _coding_contrast: median per-base frame contrast over sequences >= 30 nt."""
    vals = sorted(model.contrast(s, stops) / len(s) for s in seqs if len(s) >= 30)
    if not vals:
        return None
    m = len(vals) // 2
    return vals[m] if len(vals) % 2 else (vals[m - 1] + vals[m]) / 2


def per_kmer_dict(model):
    """framed_kmer -> score, for comparison with hexamer.scores."""
    out = {}
    for k, t in model.tables.items():
        for f in range(3):
            nz = np.nonzero(t[f])[0]
            for idx in nz:
                s = ""
                x = int(idx)
                for _ in range(k):
                    s = "ACGT"[x % 4] + s
                    x //= 4
                out[f"{s}-{f}"] = t[f, idx]
    return out

