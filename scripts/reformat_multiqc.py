#!/usr/bin/env python3
"""
reformat_qc_tsv.py
------------------
Reformats a MultiQC-style QC TSV so there is one row per sample with:
  - Pre_trim_reads   : (Seqs from pre-trim R1 row) * 2 * 1,000,000  [integer]
  - Post_trim_reads  : (Seqs from post-trim R1_paired row) * 2 * 1,000,000  [integer]
  - Alignment_pct    : % Aligned from the trimmed summary row
  - Mean_depth       : mean per-base depth across entire genome (mosdepth)
  - Breadth_pct      : % of genome covered at >=1x depth (mosdepth)

Usage:
    python reformat_qc_tsv.py <input.tsv> [options] > output.tsv

Options:
    --bam-dir DIR         Directory containing trimmed_<sample>.bam files [default: .]
    --threads N           Threads per mosdepth job [default: 4]
    --jobs N              Number of mosdepth jobs to run in parallel [default: 4]
                          Set to (total CPUs allocated) / (--threads).
    --resume PARTIAL_TSV  Existing (incomplete) output TSV. Rows with no NAs are
                          kept as-is; only rows containing NA are recomputed.

Example (resuming a previous run with 16 CPUs):
    srun --cpus-per-task=16 python reformat_qc_tsv.py qc.tsv \\
        --jobs 4 --threads 4 --resume partial_output.tsv > out.tsv
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


# ---------------------------------------------------------------------------
# mosdepth
# ---------------------------------------------------------------------------

def _mosdepth_worker(args: tuple) -> tuple:
    """
    Run mosdepth on one BAM. Returns (bam_stem, mean_depth, breadth_pct).
    """
    bam_path, threads, tmpdir = args
    bam = Path(bam_path)
    prefix = os.path.join(tmpdir, bam.stem)

    cmd = [
        "mosdepth",
        "--fast-mode",         # skip per-base depth output file (much faster)
        "--quantize", "0:1:",  # label regions: uncovered=0:1, covered=anything else
        "-t", str(threads),
        prefix,
        bam_path,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[mosdepth ERROR] {bam.name}:\n{result.stderr.strip()}", file=sys.stderr)
        return bam.stem, "NA", "NA"

    # --- Mean depth from summary file ------------------------------------
    mean_depth = "NA"
    total_length = 0
    summary_path = f"{prefix}.mosdepth.summary.txt"
    try:
        with open(summary_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row["chrom"] == "total":
                    mean_depth = f"{float(row['mean']):.4f}"
                    total_length = int(row["length"])
                    break
    except (FileNotFoundError, KeyError, ValueError) as e:
        print(f"[mosdepth summary parse error] {bam.name}: {e}", file=sys.stderr)

    # --- Breadth of coverage from quantized BED --------------------------
    # Any label other than "0:1" means covered at >=1x
    breadth_pct = "NA"
    quantized_path = f"{prefix}.quantized.bed.gz"
    if os.path.exists(quantized_path) and total_length > 0:
        covered_bases = 0
        try:
            with gzip.open(quantized_path, "rt") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 4 and parts[3] != "0:1":
                        covered_bases += int(parts[2]) - int(parts[1])
            breadth_pct = f"{covered_bases / total_length * 100:.2f}"
        except Exception as e:
            print(f"[mosdepth breadth parse error] {bam.name}: {e}", file=sys.stderr)

    return bam.stem, mean_depth, breadth_pct


def run_mosdepth_parallel(bam_paths, threads, jobs, tmpdir):
    """
    Run mosdepth on all BAMs in parallel.
    Returns dict: bam_stem -> (mean_depth, breadth_pct)
    """
    results = {}
    work = [(bam, threads, tmpdir) for bam in bam_paths]
    total = len(work)

    with ProcessPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(_mosdepth_worker, w): w[0] for w in work}
        done = 0
        for future in as_completed(futures):
            bam_path = futures[future]
            done += 1
            try:
                stem, mean, breadth = future.result()
                results[stem] = (mean, breadth)
                print(
                    f"  [{done}/{total}] {Path(bam_path).name}  "
                    f"mean={mean}x  breadth={breadth}%",
                    file=sys.stderr,
                )
            except Exception as e:
                stem = Path(bam_path).stem
                results[stem] = ("NA", "NA")
                print(f"  [{done}/{total}] {Path(bam_path).name} FAILED: {e}", file=sys.stderr)

    return results


# ---------------------------------------------------------------------------
# TSV helpers
# ---------------------------------------------------------------------------

def seqs_to_reads(seqs_str):
    """Convert a Seqs value (millions, one strand) to total paired reads."""
    if seqs_str in (".", "", "NA"):
        return "NA"
    try:
        return str(round(float(seqs_str) * 2 * 1_000_000))
    except ValueError:
        return "NA"


def load_tsv(tsv_file):
    rows = {}
    with open(tsv_file, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        # Auto-detect the sample column: MultiQC may name it "Sample",
        # "SampleName", "# Sample", or similar.
        sample_col = None
        for candidate in reader.fieldnames or []:
            if "sample" in candidate.strip().lstrip("#").strip().lower():
                sample_col = candidate
                break

        if sample_col is None:
            print(
                f"ERROR: Could not find a sample name column.\n"
                f"Columns found: {reader.fieldnames}\n"
                f"Tip: check 'head -1 {tsv_file} | cat -A' for hidden characters.",
                file=sys.stderr,
            )
            sys.exit(1)

        print(f"Using column '{sample_col}' as sample names.", file=sys.stderr)

        for row in reader:
            rows[row[sample_col]] = row
    return rows


def discover_base_names(rows):
    """
    Return base sample names in input order, handling samples that have
    no summary row (only _R1/_R2 rows).
    """
    base_names = []
    seen = set()

    for sample in rows:
        if sample.startswith("trimmed_"):
            continue
        if sample.endswith("_R1") or sample.endswith("_R2"):
            base = sample[:-3]  # strip _R1 or _R2
        else:
            base = sample

        if base not in seen:
            base_names.append(base)
            seen.add(base)

    return base_names


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Reformat QC TSV to one row per sample with mosdepth coverage."
    )
    parser.add_argument("tsv", help="Input MultiQC general stats TSV file")
    parser.add_argument("--bam-dir", default=".",
                        help="Directory with trimmed_*.bam files [.]")
    parser.add_argument("--threads", type=int, default=4,
                        help="mosdepth threads per job [4] (saturates ~4; raise --jobs instead)")
    parser.add_argument("--jobs", type=int, default=4,
                        help="Parallel mosdepth jobs [4]. Tip: set to (total CPUs) / (--threads).")
    parser.add_argument("--resume", metavar="PARTIAL_TSV",
                        help="Existing (incomplete) output TSV. Rows with no NAs are kept "
                             "as-is; only rows containing NA in any column are recomputed.")
    args = parser.parse_args()

    if not os.path.isfile(args.tsv):
        print(f"ERROR: TSV not found: {args.tsv}", file=sys.stderr)
        sys.exit(1)

    # --- Load partial results if resuming --------------------------------
    completed = {}  # sample -> completed output row (no NAs)
    if args.resume:
        if not os.path.isfile(args.resume):
            print(f"ERROR: Resume file not found: {args.resume}", file=sys.stderr)
            sys.exit(1)
        with open(args.resume, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if "NA" not in row.values():
                    completed[row["Sample"]] = row
        print(
            f"Resuming: {len(completed)} sample(s) already complete, will skip them.",
            file=sys.stderr,
        )

    # --- Load TSV and find sample names ----------------------------------
    rows = load_tsv(args.tsv)
    base_names = discover_base_names(rows)
    if not base_names:
        print("ERROR: No base sample names found. Check the 'Sample' column.", file=sys.stderr)
        sys.exit(1)

    # --- Collect BAMs — skip samples already completed -------------------
    bam_paths = []
    missing_bams = []
    for base in base_names:
        if base in completed:
            continue  # already have good data for this sample
        bam = os.path.join(args.bam_dir, f"trimmed_{base}.bam")
        if os.path.isfile(bam):
            bam_paths.append(bam)
        else:
            missing_bams.append(base)

    if missing_bams:
        print(
            f"[warning] No BAM found for {len(missing_bams)} sample(s); "
            f"depth will be NA: {', '.join(missing_bams)}",
            file=sys.stderr,
        )

    # --- Run mosdepth in parallel ----------------------------------------
    depth_results = {}
    if bam_paths:
        print(
            f"Running mosdepth on {len(bam_paths)} BAMs "
            f"({args.jobs} jobs x {args.threads} threads) ...",
            file=sys.stderr,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            depth_results = run_mosdepth_parallel(bam_paths, args.threads, args.jobs, tmpdir)
    elif not completed:
        print("[warning] No BAMs found; depth columns will all be NA.", file=sys.stderr)

    # --- Write output ----------------------------------------------------
    out_fields = [
        "Sample",
        "Pre_trim_reads",
        "Post_trim_reads",
        "Alignment_pct",
        "Mean_depth",
        "Breadth_pct",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=out_fields, delimiter="\t")
    writer.writeheader()

    for base in base_names:
        # Use cached row if already complete
        if base in completed:
            writer.writerow(completed[base])
            continue

        pre_reads  = seqs_to_reads(rows.get(f"{base}_R1", {}).get("Seqs", "."))
        post_reads = seqs_to_reads(rows.get(f"trimmed_{base}_R1_paired", {}).get("Seqs", "."))

        align_pct = rows.get(f"trimmed_{base}", {}).get("% Aligned", ".")
        if align_pct in (".", ""):
            align_pct = "NA"

        bam_stem = f"trimmed_{base}"
        mean_depth, breadth_pct = depth_results.get(bam_stem, ("NA", "NA"))

        writer.writerow({
            "Sample":          base,
            "Pre_trim_reads":  pre_reads,
            "Post_trim_reads": post_reads,
            "Alignment_pct":   align_pct,
            "Mean_depth":      mean_depth,
            "Breadth_pct":     breadth_pct,
        })

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()