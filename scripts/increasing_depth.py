#!/usr/bin/env python3
"""
find_increasing_depth.py
------------------------
Finds genomic positions whose mean read depth is non-decreasing across all
time points.  Uses DuckDB to JOIN and filter the TSV files entirely on disk —
no full DataFrames are ever loaded into RAM.

Outputs:
  - nondecreasing_positions.tsv   (qualifying positions + all depth columns)
  - nondecreasing_depth.png / .pdf (spaghetti plot)
  - stdout summary

Usage:
    python find_increasing_depth.py                        # whole genome
    python find_increasing_depth.py --chrom CM040203.1     # single chrom
    python find_increasing_depth.py --datadir /data --outdir results/

Dependencies:
    pip install duckdb pandas matplotlib numpy
"""

import argparse
import os
import sys
import time

import duckdb
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


TIMEPOINTS    = [2002, 2008, 2010, 2012, 2017, 2019]
FILE_TEMPLATE = "merged_depth_{year}.tsv"


# ── DuckDB query ──────────────────────────────────────────────────────────────

def build_query(file_paths, chrom, min_pct_increase=0.0):
    """
    Build a SQL query that:
      1. Reads all TSV files directly from disk (DuckDB scans them as views)
      2. INNER JOINs on (chromosome, position) — only positions in ALL files
      3. Filters for non-decreasing depth across time points
      4. Optionally restricts to a single chromosome

    Nothing is pulled into Python RAM until the final small result set.
    """
    years = TIMEPOINTS

    # Reference each file as a read_csv scan
    froms = {
        year: f"read_csv_auto('{path}', delim='\\t', header=true)"
        for year, path in zip(years, file_paths)
    }

    # Build the SELECT clause
    select_cols = ["t0.chromosome", "t0.position"] + [
        f"t{i}.mean_depth AS depth_{year}"
        for i, year in enumerate(years)
    ]

    # Build the JOIN chain (t0 is base; join t1…t5 on chrom+pos)
    from_clause  = f"{froms[years[0]]} AS t0"
    join_clauses = "\n    ".join(
        f"INNER JOIN {froms[year]} AS t{i} "
        f"ON t0.chromosome = t{i}.chromosome AND t0.position = t{i}.position"
        for i, year in enumerate(years[1:], start=1)
    )

    # Non-decreasing: each consecutive pair >=
    depth_aliases    = [f"t{i}.mean_depth" for i in range(len(years))]
    # Each step must increase by at least min_pct_increase %
    multiplier = 1.0 + (min_pct_increase / 100.0)
    nondec_conditions = " AND ".join(
        f"{depth_aliases[i+1]} >= {depth_aliases[i]} * {multiplier}"
        for i in range(len(years) - 1)
    )

    chrom_condition = f"AND t0.chromosome = '{chrom}'" if chrom else ""

    query = f"""
    SELECT
        {', '.join(select_cols)}
    FROM
        {from_clause}
    {join_clauses}
    WHERE
        {nondec_conditions}
        {chrom_condition}
    ORDER BY
        t0.chromosome, t0.position
    """
    return query


def run_query(file_paths, chrom, n_threads, min_pct_increase=0.0):
    """Execute the DuckDB query and return a DataFrame of qualifying rows only."""
    con = duckdb.connect()                       # in-memory session; no temp DB on disk
    con.execute(f"PRAGMA threads={n_threads};")  # parallelise DuckDB's own execution

    query = build_query(file_paths, chrom, min_pct_increase)

    print("[duckdb] Running on-disk JOIN + non-decreasing filter …")
    print("[duckdb] (Full files are never loaded into RAM)")

    result_df = con.execute(query).df()          # only qualifying rows come back
    con.close()
    return result_df


def count_total_positions(file_paths, chrom, n_threads):
    """Count positions shared across all time points (for the summary %)."""
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={n_threads};")

    froms = {
        year: f"read_csv_auto('{path}', delim='\\t', header=true)"
        for year, path in zip(TIMEPOINTS, file_paths)
    }
    from_clause  = f"{froms[TIMEPOINTS[0]]} AS t0"
    join_clauses = "\n    ".join(
        f"INNER JOIN {froms[year]} AS t{i} "
        f"ON t0.chromosome = t{i}.chromosome AND t0.position = t{i}.position"
        for i, year in enumerate(TIMEPOINTS[1:], start=1)
    )
    chrom_condition = f"WHERE t0.chromosome = '{chrom}'" if chrom else ""

    query = f"""
    SELECT COUNT(*) AS n
    FROM   {from_clause}
    {join_clauses}
    {chrom_condition}
    """
    total = con.execute(query).fetchone()[0]
    con.close()
    return total


# ── Plotting ──────────────────────────────────────────────────────────────────

def spaghetti_plot(qualifying_df, depth_cols, years, outdir, chrom_label,min_pct_increase):
    n = len(qualifying_df)
    print(f"\n[plot] Drawing spaghetti plot for {n:,} qualifying positions …")

    fig, ax = plt.subplots(figsize=(10, 6))

    cmap         = cm.get_cmap("tab20" if n <= 20 else "hsv")
    alpha        = max(0.01, min(0.3, 50 / max(n, 1)))

    depth_matrix = qualifying_df[depth_cols].values

    for i in range(n):
        ax.plot(years, depth_matrix[i],
                color=cmap(i / max(n - 1, 1)), alpha=alpha, linewidth=0.6)

    # Median overlay
    median_depths = np.median(depth_matrix, axis=0)
    ax.plot(years, median_depths, color="white",   linewidth=3.0, zorder=5)
    ax.plot(years, median_depths, color="crimson", linewidth=2.0, zorder=6, label="Median")

    ax.set_xticks(years)
    ax.set_xticklabels([str(y) for y in years], fontsize=11)
    ax.set_xlabel("Year", fontsize=13)
    ax.set_ylabel("Mean Read Depth", fontsize=13)

    title_suffix = chrom_label if chrom_label else "Whole Genome"
    ax.set_title(
        f"Non-decreasing Depth Positions — {title_suffix}\n"
        f"({n:,} positions  |  crimson = median)",
        fontsize=13, pad=12,
    )

    # Dark theme
    ax.set_facecolor("#1a1a2e")
    fig.patch.set_facecolor("#12122a")
    for obj in [ax.xaxis.label, ax.yaxis.label, ax.title]:
        obj.set_color("white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#555577")

    fig.tight_layout()

    png_path = os.path.join(outdir, f"nondecreasing_depth_{min_pct_increase}.png")
    # pdf_path = os.path.join(outdir, "nondecreasing_depth.pdf")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    # fig.savefig(pdf_path,           bbox_inches="tight")
    plt.close(fig)

    print(f"[plot] Saved → {png_path}")
    # print(f"[plot] Saved → {pdf_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Find genomic positions with non-decreasing depth (DuckDB on-disk)."
    )
    parser.add_argument("--chrom",   default=None,
                        help="Restrict to one chromosome, e.g. CM040203.1")
    parser.add_argument("--datadir", default=".",
                        help="Directory containing merged_depth_YYYY.tsv files.")
    parser.add_argument("--outdir",  default=".",
                        help="Directory for output files.")
    parser.add_argument("--min-pct-increase", type=float, default=0.0,
                    help="Minimum %% increase in depth between each consecutive time point (e.g. 5 = 5%%).")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    n_cpu = os.cpu_count() or 1
    print(f"[threads] Detected {n_cpu} logical CPU(s) — using all of them.")

    # Resolve + validate file paths
    file_paths = []
    for year in TIMEPOINTS:
        path = os.path.join(args.datadir, FILE_TEMPLATE.format(year=year))
        if not os.path.isfile(path):
            sys.exit(f"[error] File not found: {path}")
        file_paths.append(path)

    t0 = time.perf_counter()

    # ── On-disk JOIN + filter via DuckDB ──────────────────────────────────────
    qualifying   = run_query(file_paths, args.chrom, n_cpu, args.min_pct_increase)
    n_qualifying = len(qualifying)

    # ── Count total shared positions (second DuckDB pass, still on-disk) ──────
    print("[duckdb] Counting total shared positions for summary …")
    total_positions = count_total_positions(file_paths, args.chrom, n_cpu)

    elapsed     = time.perf_counter() - t0
    depth_cols  = [f"depth_{y}" for y in TIMEPOINTS]
    chrom_label = args.chrom or "whole genome"

    # ── Stdout summary ────────────────────────────────────────────────────────
    print("\n" + "═" * 58)
    print("  SUMMARY")
    print("═" * 58)
    print(f"  Region analysed     : {chrom_label}")
    print(f"  Total positions     : {total_positions:,}")
    print(f"  Non-decreasing      : {n_qualifying:,}  "
          f"({100 * n_qualifying / max(total_positions, 1):.2f} %)")
    print(f"  Min % increase/step : {args.min_pct_increase:.1f} %")
    if n_qualifying > 0:
        chromosomes = qualifying["chromosome"].unique()
        print(f"  Chromosomes hit     : {len(chromosomes)}  "
              f"({', '.join(sorted(chromosomes)[:5])}"
              f"{'…' if len(chromosomes) > 5 else ''})")
        dm = qualifying[depth_cols].values
        print(f"  Median depth range  : "
              f"{np.median(dm[:, 0]):.3f}  →  {np.median(dm[:, -1]):.3f}")
        print(f"  Max depth at {TIMEPOINTS[-1]}   : {dm[:, -1].max():.3f}")
    print(f"  Wall time           : {elapsed:.1f} s")
    print("═" * 58)

    if n_qualifying == 0:
        print("\n[info] No non-decreasing positions found — no output files written.")
        sys.exit(0)

    # ── Save TSV ──────────────────────────────────────────────────────────────
    tsv_path = os.path.join(args.outdir, f"nondecreasing_positions_{args.min_pct_increase}.tsv")
    qualifying.to_csv(tsv_path, sep="\t", index=False)
    print(f"\n[tsv]  Saved → {tsv_path}")

    # ── Plot ──────────────────────────────────────────────────────────────────
    spaghetti_plot(qualifying, depth_cols, TIMEPOINTS, args.outdir, args.chrom, args.min_pct_increase)

    print(f"\n[done] Total wall time: {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()