#!/usr/bin/env python3
"""Kozak (start-codon context) logos for the benchmark organisms.

usage: make_kozak_figure.py <inputs_dir> <runs_dir> <out_dir>
       (PYTHONPATH must include a PyTransDecoder tree, for pytransdecoder.core.kozak)

For each organism, on the stranded 'complete' benchmark input:
  self:  the model TransDecoder trains: starts of its 500 training ORFs (the run's
         longest_orfs.cds.top_500_longest) vs in-frame ATGs inside them
  truth: annotated CDS starts of genes with an annotated (not added) 5' UTR vs in-frame ATGs
         inside those CDSs

Outputs (PDF and PNG):
  kozak_logos                 log2-odds logos of the self-trained models, one y scale
  kozak_logos_self_vs_truth   self-trained and annotation-trained log2-odds logos side by side
  kozak_logos_bits            information content of the self-trained start contexts
  kozak_stats.tsv             training starts, Pearson r of self vs truth log-odds (8 positions x
                              4 bases), and the -3 A/G and +4 G log2-odds
"""
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

from pytransdecoder.core import kozak

ORGS = [("human", "Homo sapiens"), ("mouse", "Mus musculus"), ("dmel", "Drosophila melanogaster"),
        ("athal", "Arabidopsis thaliana"), ("spombe", "Schizosaccharomyces pombe"),
        ("tthermophila", "Tetrahymena thermophila"), ("pfalciparum", "Plasmodium falciparum"),
        ("creinhardtii", "Chlamydomonas reinhardtii")]
COLORS = {"A": "#109648", "C": "#255C99", "G": "#F7B32B", "T": "#D62839"}
COMP = str.maketrans("ACGT", "TGCA")
FONT = FontProperties(family="DejaVu Sans", weight="bold")
# x positions: -6..-1 at 0..5, ATG at 6..8, +4/+5 at 9..10
XPOS = [0, 1, 2, 3, 4, 5, 9, 10]
XLABELS = ["-6", "-5", "-4", "-3", "-2", "-1", "A", "T", "G", "+4", "+5"]


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def train(sites, seqs):
    """sites: (transcript_id, start0, end_excl, strand) of ORFs; those starting with ATG are used."""
    pos, neg, n = kozak.new_counts(), kozak.new_counts(), 0
    for tid, start, end, strand in sites:
        seq = seqs[tid]
        if strand == "-":
            seq = seq.translate(COMP)[::-1]
            start, end = len(seq) - end, len(seq) - start
        if seq[start:start + 3] != "ATG":
            continue
        n += 1
        kozak.add_site(pos, seq, start)
        for i in range(start + 3, end - 3, 3):
            if seq[i:i + 3] == "ATG":
                kozak.add_site(neg, seq, i)
    return pos, neg, n


def self_sites(training_file):
    pat = re.compile(r"(\S+):(\d+)-(\d+)\(([+-])\)")
    for line in open(training_file):
        if line.startswith(">"):
            m = pat.search(line)
            if m:
                lo, hi = sorted((int(m.group(2)), int(m.group(3))))
                yield m.group(1), lo - 1, hi, m.group(4)


def truth_sites(truth_file):
    for line in open(truth_file):
        x = line.rstrip("\n").split("\t")
        if x[0] == "transcript_id" or x[8] != "annotated":
            continue
        yield x[0], int(x[1]) - 1, int(x[2]), x[3]


def draw_letter(ax, base, x, y, height, width=0.9):
    if height <= 0:
        return
    tp = TextPath((0, 0), base, size=1, prop=FONT)
    bb = tp.get_extents()
    t = Affine2D().translate(-bb.x0, -bb.y0).scale(width / bb.width, height / bb.height) \
        .translate(x - width / 2, y)
    ax.add_patch(PathPatch(t.transform_path(tp), facecolor=COLORS[base], edgecolor="none"))


def stacks_log2odds(model):
    out = []
    for k in range(len(model)):
        lo = {b: model[k][b] / math.log(2) for b in "ACGT"}
        up = sorted([(v, b) for b, v in lo.items() if v > 0])
        dn = sorted([(-v, b) for b, v in lo.items() if v < 0])
        out.append((up, dn))
    return out


def stacks_bits(pos_counts):
    out = []
    for counts in pos_counts:
        n = sum(counts.values())
        f = {b: (counts[b] + 0.5) / (n + 2) for b in "ACGT"}
        ic = 2 + sum(f[b] * math.log2(f[b]) for b in "ACGT")
        out.append(sorted([(f[b] * ic, b) for b in "ACGT"]))
    return out


def style_axis(ax, ylim, ylabel):
    ax.set_xlim(-0.6, 10.6)
    ax.set_ylim(*ylim)
    ax.set_xticks(range(11))
    ax.set_xticklabels(XLABELS, fontsize=8)
    for i, lab in enumerate(ax.get_xticklabels()):
        if i in (6, 7, 8):
            lab.set_color("#888888")
            lab.set_fontweight("bold")
    ax.tick_params(axis="y", labelsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def plot_log2odds(ax, model, ylim):
    for x, (up, dn) in zip(XPOS, stacks_log2odds(model)):
        y = 0
        for h, b in up:
            draw_letter(ax, b, x, y, h)
            y += h
        y = 0
        for h, b in dn:
            y -= h
            draw_letter(ax, b, x, y, h)
    ax.axhline(0, color="black", lw=0.6)
    ax.axvspan(5.5, 8.5, color="#f2f2f2", zorder=0)
    style_axis(ax, ylim, "log2-odds")


def plot_bits(ax, pos_counts, ymax):
    for x, stack in zip(XPOS, stacks_bits(pos_counts)):
        y = 0
        for h, b in stack:
            draw_letter(ax, b, x, y, h)
            y += h
    ax.axvspan(5.5, 8.5, color="#f2f2f2", zorder=0)
    style_axis(ax, (0, ymax), "bits")


def extent_log2odds(model):
    return max(max(sum(h for h, _ in up), sum(h for h, _ in dn)) for up, dn in stacks_log2odds(model))


def italic(name):
    return "$\\it{" + name.replace(" ", "~") + "}$"


def save(fig, out, stem):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    inputs, runs, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)

    trained, stats = {}, []
    for org, _ in ORGS:
        seqs = fasta(inputs / org / "complete.fa")
        training = runs / org / "complete" / "perl_td" / "input.fa.transdecoder_dir" / "longest_orfs.cds.top_500_longest"
        m = {}
        for kind, sites in [("self", self_sites(training)),
                            ("truth", truth_sites(inputs / org / "complete.truth.tsv"))]:
            pos, neg, n = train(sites, seqs)
            m[kind] = {"model": kozak.log_odds_model(pos, neg), "pos": pos, "n": n}
        a = [m["self"]["model"][k][b] for k in range(8) for b in "ACGT"]
        t = [m["truth"]["model"][k][b] for k in range(8) for b in "ACGT"]
        ma, mt = sum(a) / len(a), sum(t) / len(t)
        m["r"] = sum((x - ma) * (y - mt) for x, y in zip(a, t)) / math.sqrt(
            sum((x - ma) ** 2 for x in a) * sum((y - mt) ** 2 for y in t))
        trained[org] = m

        def purine_minus3(kind):
            pos, model = m[kind]["pos"][3], m[kind]["model"][3]
            return sum(pos[b] * model[b] for b in "AG") / (pos["A"] + pos["G"]) / math.log(2)
        stats.append([org, m["self"]["n"], m["truth"]["n"], f"{m['r']:.3f}",
                      f"{purine_minus3('self'):.2f}", f"{purine_minus3('truth'):.2f}",
                      f"{m['self']['model'][6]['G'] / math.log(2):.2f}",
                      f"{m['truth']['model'][6]['G'] / math.log(2):.2f}"])

    with open(out / "kozak_stats.tsv", "w") as fh:
        fh.write("org\tn_self_starts\tn_truth_starts\tpearson_self_truth\tself_-3_AG_log2odds\t"
                 "truth_-3_AG_log2odds\tself_+4_G_log2odds\ttruth_+4_G_log2odds\n")
        for row in stats:
            fh.write("\t".join(map(str, row)) + "\n")

    # main: self-trained log2-odds, one y scale
    lim = math.ceil(max(extent_log2odds(trained[o]["self"]["model"]) for o, _ in ORGS) * 2) / 2
    fig, axes = plt.subplots(2, 4, figsize=(13, 5.2))
    for ax, (org, name) in zip(axes.flat, ORGS):
        plot_log2odds(ax, trained[org]["self"]["model"], (-lim, lim))
        ax.set_title(f"{italic(name)}\n{trained[org]['self']['n']} training starts; "
                     f"r = {trained[org]['r']:.2f} vs annotated", fontsize=9)
    fig.tight_layout()
    save(fig, out, "kozak_logos")

    # supplement: self-trained vs annotation-trained log2-odds
    lim2 = math.ceil(max(extent_log2odds(trained[o][k]["model"])
                         for o, _ in ORGS for k in ("self", "truth")) * 2) / 2
    fig, axes = plt.subplots(4, 4, figsize=(13, 11))
    for i, (org, name) in enumerate(ORGS):
        for j, kind in enumerate(("self", "truth")):
            ax = axes[i // 2][(i % 2) * 2 + j]
            plot_log2odds(ax, trained[org][kind]["model"], (-lim2, lim2))
            what = "self-trained" if kind == "self" else "annotated starts"
            ax.set_title(f"{italic(name)}\n{what} (n = {trained[org][kind]['n']:,})", fontsize=9)
    fig.tight_layout()
    save(fig, out, "kozak_logos_self_vs_truth")

    # supplement: information content of the self-trained start contexts
    ymax = math.ceil(max(sum(h for h, _ in s) for o, _ in ORGS
                         for s in stacks_bits(trained[o]["self"]["pos"])) * 4) / 4
    fig, axes = plt.subplots(2, 4, figsize=(13, 4.6))
    for ax, (org, name) in zip(axes.flat, ORGS):
        plot_bits(ax, trained[org]["self"]["pos"], ymax)
        ax.set_title(f"{italic(name)} ({trained[org]['self']['n']} training starts)", fontsize=9)
    fig.tight_layout()
    save(fig, out, "kozak_logos_bits")

    for row in stats:
        print("\t".join(map(str, row)))


if __name__ == "__main__":
    main()
