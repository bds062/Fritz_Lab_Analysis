#!/usr/bin/env python3
import sys
import math
import os
import csv

def load_chrom_map(chromref_file):
    """Load mapping from GenBank accession -> friendly name (chromN, MT, UnidentifiedContig)."""
    mapping = {}
    if not os.path.exists(chromref_file):
        return mapping

    with open(chromref_file, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            chrom_num = row.get('Chromosome name') or row.get('Chromosome') or row.get('chromosome') or row.get('Chr') or row.get('chrom')
            accession = (row.get('GenBank seq accession') or row.get('GenBank accession')
                         or row.get('GenBank') or row.get('GenBank seq acc') or row.get('GenBank Accession'))
            if not chrom_num or not accession:
                continue
            chrom_num = chrom_num.strip()
            accession = accession.strip()

            # Special cases
            if chrom_num == "Un":
                friendly = "UnidentifiedContig"
            elif chrom_num == "MT":
                friendly = "MT"
            else:
                friendly = f"chrom{chrom_num}"

            mapping[accession] = friendly
    return mapping

def compute_mean_stdev(filename):
    """One-pass streaming mean and stdev using Welford's algorithm"""
    n = 0
    mean = 0.0
    M2 = 0.0

    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            try:
                depth = int(parts[2])
            except ValueError:
                continue
            n += 1
            delta = depth - mean
            mean += delta / n
            delta2 = depth - mean
            M2 += delta * delta2

    if n < 2:
        return mean, 0.0
    variance = M2 / (n - 1)
    return mean, math.sqrt(variance)

def find_outliers(filename, mean, stdev, outfile, chrom_map):
    """Second pass: find depth > mean + 12*stdev, report ranges using mapped chrom names"""
    threshold = mean + 12 * stdev
    with open(filename) as f, open(outfile, "w") as out:
        out.write(f"File: {filename}\n")
        out.write(f"Mean depth: {mean:.2f}\n")
        out.write(f"Stdev: {stdev:.2f}\n")
        out.write(f"Threshold (mean+12*stdev): {threshold:.2f}\n\n")
        out.write("Outlier positions (depth > mean + 12*stdev):\n")

        current_chrom = None
        start_pos = None
        prev_pos = None
        ranges_found = 0

        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            orig_chrom, pos_str, depth_str = parts
            chrom = chrom_map.get(orig_chrom, orig_chrom)

            try:
                pos = int(pos_str)
                depth = int(depth_str)
            except ValueError:
                continue

            if depth > threshold:
                if current_chrom != chrom:
                    if start_pos is not None:
                        if start_pos == prev_pos:
                            out.write(f"{current_chrom}:{start_pos}\n")
                        else:
                            out.write(f"{current_chrom}:{start_pos}-{prev_pos}\n")
                        ranges_found += 1
                    current_chrom = chrom
                    start_pos = pos
                    prev_pos = pos
                else:
                    if start_pos is None:
                        start_pos = pos
                        prev_pos = pos
                    elif pos == prev_pos + 1:
                        prev_pos = pos
                    else:
                        if start_pos == prev_pos:
                            out.write(f"{chrom}:{start_pos}\n")
                        else:
                            out.write(f"{chrom}:{start_pos}-{prev_pos}\n")
                        ranges_found += 1
                        start_pos = pos
                        prev_pos = pos
            else:
                if start_pos is not None:
                    if start_pos == prev_pos:
                        out.write(f"{chrom}:{start_pos}\n")
                    else:
                        out.write(f"{chrom}:{start_pos}-{prev_pos}\n")
                    ranges_found += 1
                    start_pos = None
                    prev_pos = None

        if start_pos is not None:
            if start_pos == prev_pos:
                out.write(f"{current_chrom}:{start_pos}\n")
            else:
                out.write(f"{current_chrom}:{start_pos}-{prev_pos}\n")
            ranges_found += 1

        if ranges_found == 0:
            out.write("None found.\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python find_outliers_with_chrommap.py <depth_file1> [depth_file2 ...]")
        print("Set environment variable CHROMREF=/path/to/chromref.tsv or place chromref.tsv in cwd.")
        sys.exit(1)

    chromref = os.environ.get('CHROMREF', '/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/chromref.tsv')
    chrom_map = {}
    if os.path.exists(chromref):
        chrom_map = load_chrom_map(chromref)
        print(f"Loaded chromosome mapping from {chromref} (entries: {len(chrom_map)})")
    else:
        print(f"No chromref found at '{chromref}'. Using raw accession names.")

    for filename in sys.argv[1:]:
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue

        print(f"Processing {filename}...")
        mean, stdev = compute_mean_stdev(filename)
        outfile = os.path.splitext(filename)[0] + "_outliers12x.txt"
        find_outliers(filename, mean, stdev, outfile, chrom_map)
        print(f"Results written to {outfile}")

if __name__ == "__main__":
    main()