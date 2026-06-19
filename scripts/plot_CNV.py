import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import re

# ── Config ────────────────────────────────────────────────────────────────────

INPUT_FILE = "out_mergeCNVR"
TARGET_GENE = "HzeaTryp076"
CNV_THRESHOLD = 2

COLORS = {
    "low":  "#5B9BD5",
    "high": "#E06C5B",
}

# Original 9-population grouping (figs 1 & 2)
POPULATIONS_SPLIT = [
    ("Hz_LA_2002",      r"^Hz_LA_2002_"),
    ("Taylor2021_2002", r"^Taylor2021_2002_"),
    ("Hz_LA_2008",      r"^Hz_LA_2008_"),
    ("Hz_LA_2010",      r"^Hz_LA_2010_"),
    ("Hz_LA_2012",      r"^Hz_LA_2012_"),
    ("Taylor2021_2012", r"^Taylor2021_2012_"),
    ("Hz_LA_2017",      r"^Hz_LA_2017_"),
    ("Taylor2021_2017", r"^Taylor2021_2017_"),
    ("Hz_LA_2019",      r"^Hz_LA_2019_"),
]

POPULATIONS_MERGED = [
    ("2002",       [r"^Hz_LA_2002_",  r"^Taylor2021_2002_"]),
    ("2008", [r"^Hz_LA_2008_"]),
    ("2010", [r"^Hz_LA_2010_"]),
    ("2012",       [r"^Hz_LA_2012_",  r"^Taylor2021_2012_"]),
    ("2017",       [r"^Hz_LA_2017_",  r"^Taylor2021_2017_"]),
    ("2019", [r"^Hz_LA_2019_"]),
]

# ── Load & subset ─────────────────────────────────────────────────────────────

df = pd.read_csv(INPUT_FILE, sep="\t")

META_TAIL = ["average", "sd", "gene_name"]
META_HEAD = ["chr", "start", "end", "number", "gap", "repeat", "gc", "kmer"]

row = df[df["gene_name"] == TARGET_GENE]
if row.empty:
    raise ValueError(f"Gene '{TARGET_GENE}' not found in file.")
row = row.iloc[0]

sample_cols = [c for c in df.columns if c not in META_HEAD + META_TAIL]

# ── Build population → values mapping ────────────────────────────────────────

def build_pop_data(sample_cols, row, population_defs):
    """
    population_defs: list of (label, pattern_or_list_of_patterns)
    Returns an ordered dict of label → np.array of CNV values,
    and warns about any unassigned samples.
    """
    pop_data = {}
    assigned = set()

    for entry in population_defs:
        label, patterns = entry[0], entry[1]
        # Accept either a single string pattern or a list of patterns
        if isinstance(patterns, str):
            patterns = [patterns]
        members = [c for c in sample_cols if any(re.match(p, c) for p in patterns)]
        if members:
            pop_data[label] = row[members].astype(float).values
            assigned.update(members)

    unassigned = [c for c in sample_cols if c not in assigned]
    if unassigned:
        print(f"[Warning] {len(unassigned)} sample(s) not assigned to any population:")
        print("  ", unassigned)

    return pop_data


# Normalise POPULATIONS_SPLIT into the same (label, [patterns]) format
POPULATIONS_SPLIT_NORM = [(label, [pat]) for label, pat in POPULATIONS_SPLIT]

pop_data_split  = build_pop_data(sample_cols, row, POPULATIONS_SPLIT_NORM)
pop_data_merged = build_pop_data(sample_cols, row, POPULATIONS_MERGED)

pop_order_split  = [label for label, _ in POPULATIONS_SPLIT_NORM  if label in pop_data_split]
pop_order_merged = [label for label, _ in POPULATIONS_MERGED       if label in pop_data_merged]

# ── Plotting functions ────────────────────────────────────────────────────────

def plot_stacked_bar(pop_order, pop_data, title_suffix, out_path):
    fig, ax = plt.subplots(figsize=(11, 5))
    bar_w = 0.6

    for i, pop in enumerate(pop_order):
        vals    = pop_data[pop]
        n_total = len(vals)
        n_high  = (vals > CNV_THRESHOLD).sum()
        n_low   = n_total - n_high
        prop    = n_high / n_total if n_total > 0 else 0

        ax.bar(i, n_low,  bar_w, color=COLORS["low"])
        ax.bar(i, n_high, bar_w, color=COLORS["high"], bottom=n_low)

        ax.text(i, n_total + 0.3, f"n={n_total}",
                ha="center", va="bottom", fontsize=8, color="0.35")
        ax.text(i, n_total - 1.1, f"{prop:.0%}",
                ha="center", va="top", fontsize=7.5, color="0.25")

    ax.set_xticks(np.arange(len(pop_order)))
    ax.set_xticklabels(pop_order, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Number of individuals")
    ax.set_ylim(0, 30)
    ax.set_title(f"{TARGET_GENE} — CNV call distribution by population {title_suffix}")
    ax.yaxis.grid(True, linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)

    legend_patches = [
        mpatches.Patch(color=COLORS["low"],  label=f"CNV ≤ {CNV_THRESHOLD}"),
        mpatches.Patch(color=COLORS["high"], label=f"CNV > {CNV_THRESHOLD}"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", framealpha=0.9)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")
    return fig


def plot_boxplot(pop_order, pop_data, title_suffix, out_path):
    fig, ax = plt.subplots(figsize=(11, 5))
    plot_vals = [pop_data[p] for p in pop_order]

    ax.boxplot(
        plot_vals,
        positions=np.arange(len(pop_order)),
        widths=0.5,
        patch_artist=True,
        showmeans=True,
        showfliers=False,
        meanprops=dict(marker="D", markerfacecolor="white",
                       markeredgecolor="black", markersize=6, zorder=5),
        boxprops=dict(facecolor="#A1A4A7", color="#2C5F8A"),
        medianprops=dict(color="#BFD0DE", linewidth=2),
        whiskerprops=dict(color="#BBC9D4"),
        capprops=dict(color="#B4C6D5"),
    )

    rng = np.random.default_rng(42)
    for i, vals in enumerate(plot_vals):
        jitter = rng.uniform(-0.18, 0.18, size=len(vals))
        ax.scatter(i + jitter, vals, color="#E06C5B", s=18,
                   alpha=0.8, zorder=3, linewidths=0)

    ax.axhline(CNV_THRESHOLD, color="grey", linestyle="--", linewidth=1, alpha=0.7)

    ax.set_xticks(np.arange(len(pop_order)))
    ax.set_xticklabels(pop_order, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("CNV call")
    ax.set_title(f"{TARGET_GENE} — CNV call distribution per population {title_suffix}")
    ax.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")
    return fig

def plot_combined(pop_order, pop_data, out_path):
    fig, (ax_bar, ax_box) = plt.subplots(
        2, 1,
        figsize=(11, 9),
        sharex=True,
        gridspec_kw={"hspace": 0.08},   # tight gap between panels
    )

    # ── Top panel: stacked bar ────────────────────────────────────────────────
    bar_w = 0.6
    for i, pop in enumerate(pop_order):
        vals    = pop_data[pop]
        n_total = len(vals)
        n_high  = (vals > CNV_THRESHOLD).sum()
        n_low   = n_total - n_high
        prop    = n_high / n_total if n_total > 0 else 0

        ax_bar.bar(i, n_low,  bar_w, color=COLORS["low"])
        ax_bar.bar(i, n_high, bar_w, color=COLORS["high"], bottom=n_low)

        ax_bar.text(i, n_total + 0.3, f"n={n_total}",
                    ha="center", va="bottom", fontsize=8, color="0.35")
        ax_bar.text(i, n_total - 1.1, f"{prop:.0%}",
                    ha="center", va="top", fontsize=7.5, color="0.25")

    ax_bar.set_ylabel("Number of individuals")
    ax_bar.set_ylim(0, 25)
    ax_bar.set_title(f"{TARGET_GENE} — CNV call distribution by population (merged)")
    ax_bar.yaxis.grid(True, linestyle="--", alpha=0.5)
    ax_bar.set_axisbelow(True)
    ax_bar.legend(
        handles=[
            mpatches.Patch(color=COLORS["low"],  label=f"CNV ≤ {CNV_THRESHOLD}"),
            mpatches.Patch(color=COLORS["high"], label=f"CNV > {CNV_THRESHOLD}"),
        ],
        loc="upper right", framealpha=0.9,
    )

    # ── Bottom panel: boxplot ─────────────────────────────────────────────────
    plot_vals = [pop_data[p] for p in pop_order]

    ax_box.boxplot(
        plot_vals,
        positions=np.arange(len(pop_order)),
        widths=0.5,
        patch_artist=True,
        showmeans=True,
        showfliers=False,
        meanprops=dict(marker="D", markerfacecolor="white",
                       markeredgecolor="black", markersize=6, zorder=5),
        boxprops=dict(facecolor="#A8C8E8", color="#2C5F8A"),
        medianprops=dict(color="#2C5F8A", linewidth=2),
        whiskerprops=dict(color="#2C5F8A"),
        capprops=dict(color="#2C5F8A"),
    )

    rng = np.random.default_rng(42)
    for i, vals in enumerate(plot_vals):
        jitter = rng.uniform(-0.18, 0.18, size=len(vals))
        ax_box.scatter(i + jitter, vals, color="#E06C5B", s=18,
                       alpha=0.55, zorder=3, linewidths=0)

    ax_box.axhline(CNV_THRESHOLD, color="grey", linestyle="--", linewidth=1, alpha=0.7)
    ax_box.set_ylabel("CNV call")
    ax_box.yaxis.grid(True, linestyle="--", alpha=0.4)
    ax_box.set_axisbelow(True)

    # Shared x-axis labels only on the bottom panel
    ax_box.set_xticks(np.arange(len(pop_order)))
    ax_box.set_xticklabels(pop_order, rotation=35, ha="right", fontsize=9)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")
    return fig


# ── Generate all four figures ─────────────────────────────────────────────────

plot_stacked_bar(pop_order_split,  pop_data_split,  "(split)",  f"fig1_{TARGET_GENE}_color.png")
plot_boxplot    (pop_order_split,  pop_data_split,  "(split)",  f"fig2_{TARGET_GENE}_color.png")
plot_stacked_bar(pop_order_merged, pop_data_merged, "(merged)", f"fig3_{TARGET_GENE}_color.png")
plot_boxplot    (pop_order_merged, pop_data_merged, "(merged)", f"fig4_{TARGET_GENE}_color.png")
plot_combined   (pop_order_merged, pop_data_merged,             "fig5.png")

plt.show()