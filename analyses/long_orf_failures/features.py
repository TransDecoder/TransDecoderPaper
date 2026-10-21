"""Sequence features for annotated CDSs: GC, GC3, codon usage, amino-acid composition and simple
low-complexity measures, plus annotation (gene name / description / chromosome) parsing."""
import gzip
import math
import re
from urllib.parse import unquote

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
_AA_IDX = {a: i for i, a in enumerate(AA)}
# disorder-promoting residues (Dunker et al. 2001)
DISORDER = set("ARGQSPEK")
LOG2 = math.log(2)


def gc_content(seq):
    acgt = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / acgt if acgt else float("nan")


def gc3(cds):
    third = cds[2::3]
    return gc_content(third)


def codon_counts(cds, stops):
    """Sense codons (ACGT only) in frame 0, excluding stop codons."""
    out = {}
    for i in range(0, len(cds) - 2, 3):
        c = cds[i:i + 3]
        if c in stops or any(b not in "ACGT" for b in c):
            continue
        out[c] = out.get(c, 0) + 1
    return out


def aa_entropy(prot):
    """Shannon entropy (bits) of the amino-acid composition."""
    n = len(prot)
    if not n:
        return float("nan")
    counts = {}
    for a in prot:
        counts[a] = counts.get(a, 0) + 1
    return -sum(c / n * math.log2(c / n) for c in counts.values())


def top_residues(prot):
    n = len(prot)
    counts = {}
    for a in prot:
        counts[a] = counts.get(a, 0) + 1
    ranked = sorted(counts.items(), key=lambda x: -x[1])
    top1 = ranked[0] if ranked else ("", 0)
    top2 = sum(c for _, c in ranked[:2])
    top3 = sum(c for _, c in ranked[:3])
    return top1[0], top1[1] / n if n else 0.0, top2 / n if n else 0.0, top3 / n if n else 0.0, \
        "".join(a for a, _ in ranked[:3])


def seg_like_fraction(prot, window=12, max_bits=2.2):
    """Fraction of residues covered by a window of `window` residues whose composition entropy
    (bits) is below `max_bits` (the SEG trigger-window idea, K1 = 2.2 bits, W = 12)."""
    n = len(prot)
    if n < window:
        return 0.0
    idx = np.array([_AA_IDX.get(a, 20) for a in prot])
    onehot = np.zeros((n + 1, 21), dtype=np.int32)
    onehot[np.arange(1, n + 1), idx] = 1
    cs = np.cumsum(onehot, axis=0)
    wc = (cs[window:] - cs[:-window]).astype(float) / window       # (n-window+1, 21)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(wc > 0, wc * np.log2(wc), 0.0), axis=1)
    low = ent < max_bits
    covered = np.zeros(n + 1, dtype=np.int32)
    starts = np.nonzero(low)[0]
    np.add.at(covered, starts, 1)
    np.add.at(covered, starts + window, -1)
    return float((np.cumsum(covered)[:n] > 0).mean())


def wmer_unique_pct(prot, w=5):
    """TransDecoder's training-set low-complexity measure (exclude_similar_proteins): percent of
    distinct non-overlapping 5-residue words among all words; < 30 is skipped as low complexity."""
    words = {prot[i:i + w] for i in range(0, len(prot) - w + 1, w)}
    if not words:
        return float("nan")
    return len(words) / (len(prot) / w) * 100


def tandem_repeat_fraction(prot, k=3):
    """Fraction of k-mer positions whose k-mer recurs within the next 30 residues: a crude repeat
    measure that is high for collagen-like, zinc-finger-array and other repetitive proteins."""
    n = len(prot) - k + 1
    if n <= 0:
        return 0.0
    last = {}
    hits = 0
    for i in range(n):
        km = prot[i:i + k]
        j = last.get(km)
        if j is not None and i - j <= 30:
            hits += 1
        last[km] = i
    return hits / n


def aa_composition(prot):
    n = len(prot)
    comp = np.zeros(20)
    for a in prot:
        j = _AA_IDX.get(a)
        if j is not None:
            comp[j] += 1
    return comp / n if n else comp


def protein_features(prot):
    prot = prot.rstrip("*")
    t1, t1f, t2f, t3f, top3 = top_residues(prot)
    comp = aa_composition(prot)
    f = dict(aa_len=len(prot), aa_entropy=aa_entropy(prot), top1_aa=t1, top1_frac=t1f,
             top2_frac=t2f, top3_frac=t3f, top3_aa=top3, seg_lowc_frac=seg_like_fraction(prot),
             wmer_unique_pct=wmer_unique_pct(prot), tandem3_frac=tandem_repeat_fraction(prot),
             disorder_frac=sum(prot.count(a) for a in DISORDER) / len(prot) if prot else float("nan"),
             internal_stops=prot.count("*"))
    for a, v in zip(AA, comp):
        f[f"aa_{a}"] = v
    return f


# ---------------------------------------------------------------- annotation

ORGANELLE_RE = re.compile(r"^(chrM|MT|M|ChrC|ChrM|mitochondrion|chloroplast|mito|Pt|NC_037304|"
                          r"NC_036769|NC_001799|NC_005353)", re.I)


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _gff3_attrs(text):
    out = {}
    for kv in text.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = unquote(v)
    return out


def _gtf_attrs(text):
    return {m.group(1): m.group(2) for m in re.finditer(r'(\S+) "([^"]*)"', text)}


def read_annotation(path):
    """transcript_id -> (gene_name, description, seqid). Works for the GTF/GFF3 files of the
    bench_v3 truth sets (MANE GTF, GENCODE GFF3, FlyBase GTF, TAIR/PomBase/RefSeq/TGD GFF3)."""
    is_gtf = ".gtf" in str(path)
    genes, tx = {}, {}
    for line in _open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = _gtf_attrs(f[8]) if is_gtf else _gff3_attrs(f[8])
        if f[2] == "gene" and not is_gtf and "ID" in a:
            genes[a["ID"]] = a
        if f[2] in ("mRNA", "transcript"):
            tid = a.get("transcript_id") if is_gtf else a.get("ID")
            if tid:
                tx[tid] = (a, f[0])
            if not is_gtf and a.get("transcript_id"):
                tx.setdefault(a["transcript_id"], (a, f[0]))
    out = {}
    for tid, (a, seqid) in tx.items():
        g = genes.get(a.get("Parent", ""), {})
        name = (a.get("gene_name") or a.get("gene_symbol") or g.get("Name") or a.get("locus_tag")
                or g.get("locus_tag") or a.get("Parent", ""))
        desc = a.get("product") or g.get("product") or a.get("Note") or g.get("Note") or \
            g.get("description") or ""
        out[tid] = (name, desc, seqid)
    return out


def is_organelle(seqid):
    return bool(ORGANELLE_RE.match(seqid or ""))
