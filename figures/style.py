"""Shared figure style for the manuscript (Genome Biology: max width 17 cm, 7-8 pt text).

usage (from figures/<dir>/make_*.py):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import style
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1 / 25.4
FULL_W = 170 * MM          # double column
HALF_W = 82 * MM           # single column

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 1.2,
    "lines.markersize": 3.5, "pdf.fonttype": 42, "ps.fonttype": 42,
})

# GM12878 assemblies (as in the existing drafts)
ASM_COLOR = {"lraa": "#1b9e77", "stringtie": "#d95f02", "trinity": "#7570b3"}
ASM_LABEL = {"lraa": "LRAA (long reads)", "stringtie": "StringTie", "trinity": "Trinity (genome-guided)"}

# Benchmark organisms, ordered by transcriptome G+C (E07 Table 1)
ORGS = ["pfalciparum", "tthermophila", "spombe", "athal", "dmel", "mouse", "human", "creinhardtii"]
ORG_LABEL = {"pfalciparum": "P. falciparum", "tthermophila": "T. thermophila", "spombe": "S. pombe",
             "athal": "A. thaliana", "dmel": "D. melanogaster", "mouse": "M. musculus",
             "human": "H. sapiens", "creinhardtii": "C. reinhardtii"}
ORG_GC = {"pfalciparum": 22.1, "tthermophila": 25.9, "spombe": 37.9, "athal": 42.2, "dmel": 50.9,
          "mouse": 51.7, "human": 52.3, "creinhardtii": 68.3}
ORG_COLOR = dict(zip(ORGS, ["#a6761d", "#e7298a", "#666666", "#66a61e", "#e6ab02", "#7570b3",
                            "#1b9e77", "#d95f02"]))

# Input types
COMPLETE_C, FRAG_C = "#2c7fb8", "#f03b20"
TD_C = "#08519c"           # TransDecoder (final) highlight
GREY = "#999999"


def panel_label(ax, letter, x=-0.18, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")


def save(fig, out_dir, stem):
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_dir / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
