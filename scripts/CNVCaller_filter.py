import pandas as pd
import numpy as np
import re
from collections import defaultdict
from itertools import combinations
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

# =============================================================================
# GLOBAL PARAMETERS
# =============================================================================

INPUT_FILE = "mergeCNVR"
OUTPUT_FILE = f"filtered_{INPUT_FILE}_14"
GFF_FILE = "/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/Liftoff_Analysis/test1.gff"

# Fixation threshold — used for start/end filter and filter 3
FIXATION_POINT = 1.25

# --- Filter 1: Start/End Fixation Check ---
USE_START_END_FILTER = True
# Exclude if >50% of 2002 samples are above FIXATION_POINT (started too high)
# Exclude if >50% of 2019 samples are below FIXATION_POINT (never reached fixation)

# --- Filter 2: Minimum Delta ---
USE_DELTA_FILTER = False
MIN_DELTA = 0.5  # Require median(2019) - median(2002) > MIN_DELTA

# --- Filter 3: Spearman Rank Correlation ---
USE_SPEARMAN_FILTER = False
SPEARMAN_RHO_MIN = 0.7  # Require Spearman rho across time point medians > this value
SPEARMAN_P_MAX = 0.05   # Optionally also require significance of the correlation

# --- Filter 4: Monotonically Non-Decreasing Fraction Above Fixation ---
USE_MONOTONIC_FILTER = True  # No threshold — binary monotonicity check across all 6 time points

# --- Filter 5: Mann-Whitney U Test (2002 vs 2019) ---
USE_MANNWHITNEY_FILTER = False
MANNWHITNEY_P_MAX = 0.05  # Require p-value < this value

# =============================================================================
# TIME POINT SAMPLE GROUPS
# =============================================================================

def get_timepoint_columns(df_columns):
    """
    Returns a dict mapping time point label -> list of column names present in df.
    """
    timepoints = {
        '2002': ['Hz_LA_2002_', 'Taylor2021_2002_'],
        '2008': ['Hz_LA_2008_'],
        '2010': ['Hz_LA_2010_'],
        '2012': ['Hz_LA_2012_', 'Taylor2021_2012_'],
        '2017': ['Hz_LA_2017_', 'Taylor2021_2017_'],
        '2019': ['Hz_LA_2019_'],
    }

    result = {}
    for tp, prefixes in timepoints.items():
        cols = [c for c in df_columns if any(c.startswith(p) for p in prefixes)]
        result[tp] = cols

    return result

# =============================================================================
# GFF PARSING
# =============================================================================

def parse_gff_genes(gff_file):
    """
    Parse GFF file and extract gene features with their coordinates and names.
    Returns dict: chrom -> list of (start, end, gene_name)
    """
    genes = defaultdict(list)

    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]

            note_match = re.search(r'Note=([^;]+)', attributes)
            if not note_match:
                continue

            gene_name = note_match.group(1)
            genes[chrom].append((start, end, gene_name))

    return genes


def find_overlapping_genes(chrom, window_start, window_end, genes_by_chrom):
    """
    Find all genes that overlap with the given window.
    Returns list of gene names.
    """
    if chrom not in genes_by_chrom:
        return []

    overlapping = []
    for gene_start, gene_end, gene_name in genes_by_chrom[chrom]:
        if not (window_end < gene_start or window_start > gene_end):
            overlapping.append(gene_name)

    return overlapping

# =============================================================================
# FILTERING FUNCTIONS
# =============================================================================

def passes_start_end_filter(row, tp_cols):
    if not USE_START_END_FILTER:
        return True, ""

    vals_2002 = row[tp_cols['2002']].values.astype(float)
    vals_2019 = row[tp_cols['2019']].values.astype(float)

    frac_2002_above = np.mean(vals_2002 > FIXATION_POINT)
    frac_2019_below = np.mean(vals_2019 < FIXATION_POINT)

    if frac_2002_above > 0.5:
        return False, f"start filter: {frac_2002_above:.2f} of 2002 samples above fixation"
    if frac_2019_below > 0.5:
        return False, f"end filter: {frac_2019_below:.2f} of 2019 samples below fixation"

    return True, ""


def passes_delta_filter(row, tp_cols):
    """
    Filter 2: Require median(2019) - median(2002) > MIN_DELTA.
    """
    if not USE_DELTA_FILTER:
        return True, ""

    median_2002 = np.median(row[tp_cols['2002']].values.astype(float))
    median_2019 = np.median(row[tp_cols['2019']].values.astype(float))
    delta = median_2019 - median_2002

    if delta <= MIN_DELTA:
        return False, f"delta filter: median delta = {delta:.3f} <= {MIN_DELTA}"

    return True, ""


def passes_spearman_filter(row, tp_cols):
    """
    Filter 3: Compute Spearman rho between time order [1..6] and per-timepoint
    medians. Require rho > SPEARMAN_RHO_MIN (and optionally p < SPEARMAN_P_MAX).
    """
    if not USE_SPEARMAN_FILTER:
        return True, ""

    ordered_tps = ['2002', '2008', '2010', '2012', '2017', '2019']
    time_order = list(range(len(ordered_tps)))
    medians = [np.median(row[tp_cols[tp]].values.astype(float)) for tp in ordered_tps]

    rho, p_val = stats.spearmanr(time_order, medians)

    if rho <= SPEARMAN_RHO_MIN:
        return False, f"Spearman filter: rho = {rho:.3f} <= {SPEARMAN_RHO_MIN}"
    if p_val >= SPEARMAN_P_MAX:
        return False, f"Spearman filter: p = {p_val:.4f} >= {SPEARMAN_P_MAX}"

    return True, ""


def passes_monotonic_filter(row, tp_cols):
    """
    Filter 4: Fraction of samples above FIXATION_POINT must be monotonically
    non-decreasing across all 6 time points.
    """
    if not USE_MONOTONIC_FILTER:
        return True, ""

    ordered_tps = ['2002', '2008', '2010', '2012', '2017', '2019']
    fractions = [
        np.mean(row[tp_cols[tp]].values.astype(float) > FIXATION_POINT)
        for tp in ordered_tps
    ]

    violations = []
    for i in range(1, len(fractions)):
        if fractions[i] < fractions[i - 1]:
            violations.append(
                f"{ordered_tps[i-1]}->{ordered_tps[i]}: {fractions[i-1]:.2f}->{fractions[i]:.2f}"
            )

    if violations:
        return False, f"monotonic filter violations: {'; '.join(violations)}"

    return True, ""


def passes_mannwhitney_filter(row, tp_cols):
    """
    Filter 5: Mann-Whitney U test between 2002 and 2019 samples.
    Require p < MANNWHITNEY_P_MAX.
    """
    if not USE_MANNWHITNEY_FILTER:
        return True, ""

    vals_2002 = row[tp_cols['2002']].values.astype(float)
    vals_2019 = row[tp_cols['2019']].values.astype(float)

    # alternative='less': tests if 2002 < 2019 (one-sided, directional)
    _, p_val = stats.mannwhitneyu(vals_2002, vals_2019, alternative='less')

    if p_val >= MANNWHITNEY_P_MAX:
        return False, f"Mann-Whitney filter: p = {p_val:.4f} >= {MANNWHITNEY_P_MAX}"

    return True, ""


def apply_all_filters(row, tp_cols):
    """
    Run all filters in order. Returns (bool, str) — passes all, reason if failed.
    """
    for filter_fn in [
        passes_start_end_filter,
        passes_delta_filter,
        passes_spearman_filter,
        passes_monotonic_filter,
        passes_mannwhitney_filter,
    ]:
        passed, reason = filter_fn(row, tp_cols)
        if not passed:
            return False, reason

    return True, ""

# =============================================================================
# INDEPENDENT FILTER EVALUATION
# =============================================================================

# Short display names for each filter function (used in UpSet plot labels)
FILTER_LABELS = {
    passes_start_end_filter: 'Endpoint Frequency Check',
    passes_delta_filter:     'delta',
    passes_spearman_filter:  'spearman',
    passes_monotonic_filter: 'Monotonic Increase of Frequency',
    passes_mannwhitney_filter: 'mann-whitney',
}

ALL_FILTER_FNS = list(FILTER_LABELS.keys())


def evaluate_filters_independently(df, tp_cols):
    """
    For every row, evaluate each active filter independently.
    Only filters with their USE_* flag set to True are included.

    Returns:
        fail_sets     — dict {filter_label: set of row indices that fail it}
        active_labels — list of labels for filters that are currently enabled
    """
    active_fns = [fn for fn in ALL_FILTER_FNS if _filter_is_active(fn)]
    active_labels = [FILTER_LABELS[fn] for fn in active_fns]

    fail_sets = {lbl: set() for lbl in active_labels}

    for idx, row in df.iterrows():
        for fn, lbl in zip(active_fns, active_labels):
            passed, _ = fn(row, tp_cols)
            if not passed:
                fail_sets[lbl].add(idx)

    return fail_sets, active_labels


def _filter_is_active(fn):
    """Return True if the USE_* flag for this filter function is enabled."""
    flags = {
        passes_start_end_filter:   USE_START_END_FILTER,
        passes_delta_filter:       USE_DELTA_FILTER,
        passes_spearman_filter:    USE_SPEARMAN_FILTER,
        passes_monotonic_filter:   USE_MONOTONIC_FILTER,
        passes_mannwhitney_filter: USE_MANNWHITNEY_FILTER,
    }
    return flags[fn]


def compute_upset_intersections(fail_sets, active_labels):
    """
    Compute exclusive intersection counts for all non-empty subsets of filters.

    A row belongs to intersection S (exclusive) if it fails every filter in S
    and passes every filter not in S.

    Returns list of (frozenset_of_labels, count) sorted descending by count.
    """
    all_indices = set().union(*fail_sets.values()) if fail_sets else set()

    intersections = []
    for r in range(1, len(active_labels) + 1):
        for combo in combinations(active_labels, r):
            combo_set = frozenset(combo)
            # rows failing all filters in combo
            exclusive = fail_sets[combo[0]].copy()
            for lbl in combo[1:]:
                exclusive &= fail_sets[lbl]
            # minus rows that also fail filters outside this combo
            for lbl in active_labels:
                if lbl not in combo_set:
                    exclusive -= fail_sets[lbl]
            if exclusive:
                intersections.append((combo_set, len(exclusive)))

    intersections.sort(key=lambda x: -x[1])
    return intersections

# =============================================================================
# UPSET PLOT
# =============================================================================

def plot_upset(intersections, fail_sets, active_labels, output_path):
    """
    Draw a publication-quality UpSet plot and save as PNG.

    Layout (top → bottom):
      [top bar chart]   — exclusive intersection counts
      [dot matrix]      — which filters participate in each intersection
      [bottom bar chart] — per-filter independent totals (horizontal, on the left)
    """
    if not intersections:
        print("No filter rejections to plot — skipping UpSet plot.")
        return

    # Show at most 15 intersections (sorted by count, already done)
    top = intersections[:15]
    n_cols = len(top)
    n_rows = len(active_labels)

    filter_totals = {lbl: len(fail_sets[lbl]) for lbl in active_labels}

    # --- Colour scheme ---
    def bar_color(n):
        if n == 1:  return '#3B8BD4'   # blue  — single filter
        if n == 2:  return '#7F77DD'   # purple — 2-way
        return             '#1D9E75'   # teal  — 3+

    # --- Figure geometry ---
    col_w   = 0.7          # inches per intersection column
    row_h   = 0.50         # inches per filter row in the dot matrix
    top_h   = 2.8          # inches for top bar chart
    side_w  = 2.2          # inches for left-side horizontal bars + labels
    pad     = 0.15         # inches between panels
    right_m = 0.3          # right margin

    fig_w = side_w + pad + n_cols * col_w + right_m
    fig_h = top_h + pad + n_rows * row_h + 0.8

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)
    fig.patch.set_facecolor('white')

    # --- Axes ---
    # left edge of the main plot area (as fraction of fig width)
    l = (side_w + 0.05) / fig_w
    b_dot = 0.55 / fig_h
    h_dot = (n_rows * row_h) / fig_h
    h_top = top_h / fig_h
    b_top = b_dot + h_dot + pad / fig_h
    w_main = 1.0 - l - right_m / fig_w

    ax_top  = fig.add_axes([l, b_top, w_main, h_top])          # top bars
    ax_dot  = fig.add_axes([l, b_dot, w_main, h_dot])          # dot matrix
    # side bars sit to the left; leave ~0.25in at the very left for the x-axis label
    ax_side = fig.add_axes([0.25 / fig_w, b_dot,
                             (side_w - 0.3) / fig_w, h_dot])   # left side bars

    # ── Top bar chart ──────────────────────────────────────────────────────────
    counts = [c for _, c in top]
    colors = [bar_color(len(s)) for s, _ in top]

    bars = ax_top.bar(range(n_cols), counts, color=colors,
                      width=0.55, zorder=3, linewidth=0)
    ax_top.set_xlim(-0.5, n_cols - 0.5)
    ax_top.set_ylim(0, max(counts) * 1.25)
    ax_top.set_xticks([])
    ax_top.spines[['right', 'top', 'bottom']].set_visible(False)
    ax_top.spines['left'].set_color('#cccccc')
    ax_top.tick_params(axis='y', labelsize=8, color='#cccccc')
    ax_top.yaxis.set_tick_params(length=3)
    ax_top.set_ylabel('CNVs removed\n(exclusive)', fontsize=8, color='#555555')
    ax_top.yaxis.grid(True, color='#eeeeee', linewidth=0.5, zorder=0)
    ax_top.set_axisbelow(True)

    for bar, count in zip(bars, counts):
        ax_top.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(counts) * 0.01,
                    str(count), ha='center', va='bottom', fontsize=7.5, color='#333333', fontweight='500')

    # ── Dot matrix ─────────────────────────────────────────────────────────────
    ax_dot.set_xlim(-0.5, n_cols - 0.5)
    ax_dot.set_ylim(-0.5, n_rows - 0.5)
    ax_dot.set_xticks([])
    ax_dot.set_yticks(range(n_rows))
    ax_dot.set_yticklabels(active_labels[::-1], fontsize=9.5, fontweight='500')
    ax_dot.tick_params(axis='y', length=0, pad=8)
    ax_dot.spines[:].set_visible(False)

    # Horizontal guide lines
    for fi in range(n_rows):
        ax_dot.axhline(fi, color='#f0f0f0', linewidth=0.8, zorder=0)

    for col_i, (combo, _) in enumerate(top):
        active_rows = sorted(
            [n_rows - 1 - active_labels.index(lbl) for lbl in combo]
        )

        # Vertical connector line through active dots
        if len(active_rows) > 1:
            color = bar_color(len(combo))
            ax_dot.plot(
                [col_i, col_i],
                [active_rows[0], active_rows[-1]],
                color=color, linewidth=2.5, zorder=2, solid_capstyle='round'
            )

        # Draw all dots (active and inactive)
        for fi in range(n_rows):
            lbl = active_labels[::-1][fi]
            active = lbl in combo
            color = bar_color(len(combo)) if active else '#e0e0e0'
            zorder = 3 if active else 1
            ax_dot.scatter(col_i, fi, s=55, color=color, zorder=zorder,
                           linewidths=0)

    # ── Left side bars (per-filter independent totals) ─────────────────────────
    max_total = max(filter_totals.values()) if filter_totals else 1
    reversed_labels = active_labels[::-1]

    ax_side.barh(range(n_rows),
                 [filter_totals[lbl] for lbl in reversed_labels],
                 color='#3B8BD4', height=0.45, linewidth=0)
    ax_side.set_ylim(-0.5, n_rows - 0.5)
    ax_side.set_yticks([])
    ax_side.set_xlim(max_total * 1.35, 0)   # reversed so bars grow rightward
    ax_side.spines[['left', 'top', 'bottom']].set_visible(False)
    ax_side.spines['right'].set_color('#cccccc')
    ax_side.tick_params(axis='x', labelsize=7, color='#cccccc')
    ax_side.xaxis.grid(True, color='#eeeeee', linewidth=0.5, zorder=0)
    ax_side.set_axisbelow(True)
    ax_side.set_xlabel('Independent\ntotal removed', fontsize=7.5, color='#555555')

    for fi, lbl in enumerate(reversed_labels):
        val = filter_totals[lbl]
        ax_side.text(val + max_total * 0.02, fi, str(val),
                     ha='left', va='center', fontsize=7, color='#333333')

    # ── Legend ─────────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color='#3B8BD4', label='Single filter'),
        mpatches.Patch(color='#7F77DD', label='2-filter overlap'),
        mpatches.Patch(color='#1D9E75', label='3+ filter overlap'),
    ]
    fig.legend(handles=legend_patches, loc='lower right',
               fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.98, 0.01))

    fig.suptitle('CNV filter efficacy — UpSet plot', fontsize=10,
                 fontweight='500', x=0.55, y=0.99, ha='center', color='#222222')

    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"UpSet plot saved to: {output_path}")


# =============================================================================
# VENN DIAGRAM  (used when ≤ 3 filters are active)
# =============================================================================

def plot_venn(intersections, fail_sets, active_labels, output_path, total_rows=None):
    """
    Publication-quality Venn diagram for 1–3 active filters.
    No title block. Labels wrap to two lines when they would otherwise overlap.
    """
    import matplotlib.colors as mcolors

    n = len(active_labels)
    if n == 0:
        print("No active filters — skipping Venn diagram.")
        return

    inter_dict    = {combo: count for combo, count in intersections}
    total_removed = len(set().union(*fail_sets.values())) if fail_sets else 0

    # ── Design tokens ─────────────────────────────────────────────────────────
    PALETTE  = ['#4878cf', '#d65f5f', '#59a14f']
    FILL_A   = 0.22
    EDGE_A   = 0.90
    EDGE_LW  = 1.5
    R        = 1.0

    COUNT_FS = 12
    NAME_FS  = 10
    TOTAL_FS = 8.5

    # Character threshold above which a label is wrapped to two lines.
    # Tune this down if your labels are still colliding.
    WRAP_THRESHOLD = 9

    plt.rcParams.update({
        'font.family':     'sans-serif',
        'font.sans-serif': ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
    })

    # ── Label wrapping ────────────────────────────────────────────────────────

    def wrap_label(lbl):
        """
        If lbl exceeds WRAP_THRESHOLD characters, split it at the separator
        character (space / - _) nearest to the midpoint.
        Separators '/' and '-' are kept at the end of the first line.
        Returns the original string or a two-line string with '\n'.
        """
        if len(lbl) <= WRAP_THRESHOLD:
            return lbl
        mid = len(lbl) // 2
        best_pos, best_sep, best_dist = -1, None, len(lbl)
        for i, ch in enumerate(lbl):
            if ch in (' ', '/', '-', '_') and abs(i - mid) < best_dist:
                best_pos, best_sep, best_dist = i, ch, abs(i - mid)
        if best_pos < 0:                          # no separator — hard break
            return lbl[:mid] + '\n' + lbl[mid:]
        if best_sep in ('/', '-'):                # keep separator on line 1
            return lbl[:best_pos + 1] + '\n' + lbl[best_pos + 1:]
        return lbl[:best_pos] + '\n' + lbl[best_pos + 1:]

    # ── Drawing helpers ───────────────────────────────────────────────────────

    def draw_circle(cx, cy, color):
        r, g, b, _ = mcolors.to_rgba(color)
        ax.add_patch(plt.Circle((cx, cy), R,
                                facecolor=(r, g, b, FILL_A),
                                linewidth=0, zorder=2))
        ax.add_patch(plt.Circle((cx, cy), R,
                                facecolor='none',
                                edgecolor=(r, g, b, EDGE_A),
                                linewidth=EDGE_LW, zorder=3))

    def count_label(x, y, value):
        if value > 0:
            ax.text(x, y, str(value), ha='center', va='center',
                    fontsize=COUNT_FS, fontweight='bold',
                    color='#111111', zorder=5)
        else:
            ax.text(x, y, '0', ha='center', va='center',
                    fontsize=COUNT_FS - 1, fontweight='normal',
                    color='#cccccc', zorder=5)

    def filter_label(x, y, lbl, color, position='above'):
        """
        Two-line label: bold colored filter name + italic gray n= total.
        When the name itself wraps to two lines, the n= is shifted further
        out so it doesn't crowd the circle edge.
        """
        n_ind   = len(fail_sets[lbl])
        wrapped = wrap_label(lbl)
        extra   = 0.24 if '\n' in wrapped else 0.0   # extra clearance for 2-line names

        if position == 'above':
            ax.text(x, y + 0.26 + extra, wrapped,
                    ha='center', va='bottom',
                    fontsize=NAME_FS, fontweight='bold',
                    color=color, linespacing=1.25, zorder=5)
            ax.text(x, y + 0.02, f'n = {n_ind}',
                    ha='center', va='bottom',
                    fontsize=TOTAL_FS, style='italic',
                    color='#555555', zorder=5)
        else:                                          # position == 'below'
            ax.text(x, y - 0.02, wrapped,
                    ha='center', va='top',
                    fontsize=NAME_FS, fontweight='bold',
                    color=color, linespacing=1.25, zorder=5)
            ax.text(x, y - 0.26 - extra, f'n = {n_ind}',
                    ha='center', va='top',
                    fontsize=TOTAL_FS, style='italic',
                    color='#555555', zorder=5)

    # ── 1 filter ──────────────────────────────────────────────────────────────
    if n == 1:
        lbl = active_labels[0]
        fig, ax = plt.subplots(figsize=(4.5, 4.0))
        fig.patch.set_facecolor('white')
        ax.set_aspect('equal')
        ax.axis('off')

        draw_circle(0, 0, PALETTE[0])
        count_label(0, 0, inter_dict.get(frozenset([lbl]), 0))
        filter_label(0, R + 0.08, lbl, PALETTE[0], position='above')

        ax.set_xlim(-1.8, 1.8)
        ax.set_ylim(-1.4, 1.9)

    # ── 2 filters ─────────────────────────────────────────────────────────────
    elif n == 2:
        A, B = active_labels
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        fig.patch.set_facecolor('white')
        ax.set_aspect('equal')
        ax.axis('off')

        off = 0.65
        draw_circle(-off, 0, PALETTE[0])
        draw_circle( off, 0, PALETTE[1])

        count_label(-off - 0.62, 0, inter_dict.get(frozenset([A]),    0))
        count_label(0,           0, inter_dict.get(frozenset([A, B]), 0))
        count_label( off + 0.62, 0, inter_dict.get(frozenset([B]),    0))

        filter_label(-off - 0.28, R + 0.08, A, PALETTE[0], position='above')
        filter_label( off + 0.28, R + 0.08, B, PALETTE[1], position='above')

        ax.set_xlim(-2.3, 2.3)
        ax.set_ylim(-1.4, 1.9)

    # ── 3 filters ─────────────────────────────────────────────────────────────
    elif n == 3:
        A, B, C = active_labels
        fig, ax = plt.subplots(figsize=(6.0, 5.5))
        fig.patch.set_facecolor('white')
        ax.set_aspect('equal')
        ax.axis('off')

        cAx, cAy = -0.55,  0.38
        cBx, cBy =  0.55,  0.38
        cCx, cCy =  0.00, -0.52

        draw_circle(cAx, cAy, PALETTE[0])
        draw_circle(cBx, cBy, PALETTE[1])
        draw_circle(cCx, cCy, PALETTE[2])

        regions = {
            frozenset([A]):       (-1.28,  0.68),
            frozenset([B]):       ( 1.28,  0.68),
            frozenset([C]):       ( 0.00, -1.62),
            frozenset([A, B]):    ( 0.00,  0.82),
            frozenset([A, C]):    (-0.70, -0.42),
            frozenset([B, C]):    ( 0.70, -0.42),
            frozenset([A, B, C]): ( 0.00, -0.06),
        }
        for combo_set, (rx, ry) in regions.items():
            count_label(rx, ry, inter_dict.get(combo_set, 0))

        filter_label(cAx - 0.48, cAy + R + 0.08, A, PALETTE[0], position='above')
        filter_label(cBx + 0.48, cBy + R + 0.08, B, PALETTE[1], position='above')
        filter_label(cCx,        cCy - R - 0.08, C, PALETTE[2], position='below')

        ax.set_xlim(-2.4, 2.4)
        ax.set_ylim(-2.4, 2.1)

    plt.savefig(output_path, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print(f"Venn diagram saved to: {output_path}")

# =============================================================================
# MAIN
# =============================================================================

# Read the TSV file
df = pd.read_csv(INPUT_FILE, sep='\t')

# Identify time point columns
tp_cols = get_timepoint_columns(df.columns)
print("Sample counts per time point:")
for tp, cols in tp_cols.items():
    print(f"  {tp}: {len(cols)} samples")

# Calculate mean and standard deviation of the 'average' column
mean_avg = df['average'].mean()
std_avg = df['average'].std()
threshold = 0

print(f"\nMean of 'average' column: {mean_avg:.4f}")
print(f"Standard deviation of 'average' column: {std_avg:.4f}")
print(f"Threshold: {threshold:.4f}")
print(f"Number of rows above threshold: {(df['average'] > threshold).sum()}")

# Pre-filter on average
filtered_df = df[df['average'] > threshold].copy()

# Parse GFF file
print("\nParsing GFF file...")
genes_by_chrom = parse_gff_genes(GFF_FILE)

# ── Independent filter evaluation (for UpSet plot) ────────────────────────────
print("\nEvaluating filters independently across all pre-filtered rows...")
fail_sets, active_labels = evaluate_filters_independently(filtered_df, tp_cols)

print("Independent rejection counts per filter:")
for lbl in active_labels:
    pct = 100 * len(fail_sets[lbl]) / len(filtered_df) if len(filtered_df) else 0
    print(f"  {lbl}: {len(fail_sets[lbl])} rows ({pct:.1f}%)")

intersections = compute_upset_intersections(fail_sets, active_labels)

print("\nExclusive intersection counts (top 15):")
for combo, count in intersections[:15]:
    print(f"  {' ∩ '.join(sorted(combo))}: {count}")

plot_path = f"{OUTPUT_FILE}.png"
if len(active_labels) <= 3:
    plot_venn(intersections, fail_sets, active_labels, plot_path)
else:
    plot_upset(intersections, fail_sets, active_labels, plot_path)

# ── Sequential filtering (for TSV output) ─────────────────────────────────────
filter_counts = defaultdict(int)
rows_with_genes = []

for idx, row in filtered_df.iterrows():
    chrom = row['chr']
    window_start = row['start']
    window_end = row['end']

    overlapping_genes = find_overlapping_genes(chrom, window_start, window_end, genes_by_chrom)

    for gene_name in overlapping_genes:
        passed, reason = apply_all_filters(row, tp_cols)

        if passed:
            row_dict = row.to_dict()
            row_dict['gene_name'] = gene_name
            rows_with_genes.append(row_dict)
        else:
            filter_counts[reason.split(':')[0].strip()] += 1

# Summary
print(f"\nSequential filter rejection summary:")
for reason, count in sorted(filter_counts.items(), key=lambda x: -x[1]):
    print(f"  {reason}: {count} rows rejected")

# Create output dataframe
if rows_with_genes:
    output_df = pd.DataFrame(rows_with_genes)
    print(f"\nRows passing all filters: {len(output_df)}")
else:
    output_df = filtered_df.copy()
    output_df['gene_name'] = []
    output_df = output_df[0:0]
    print("\nRows passing all filters: 0")

# Write TSV output
output_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
print(f"Filtered results written to: {OUTPUT_FILE}")