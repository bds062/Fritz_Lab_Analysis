#!/usr/bin/env python3
import sys
import os

def parse_outlier_file(filename):
    """Return a dict: {chrom: set(positions)} for one outlier file"""
    positions = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("File:") or line.startswith("Mean") or line.startswith("Stdev") or line.startswith("Threshold") or line.startswith("Outlier"):
                continue
            # Format: chrom:start-end OR chrom:pos
            if ":" not in line:
                continue
            chrom, coords = line.split(":")
            if "-" in coords:
                start, end = map(int, coords.split("-"))
                pos_range = range(start, end + 1)
            else:
                pos_range = [int(coords)]

            if chrom not in positions:
                positions[chrom] = set()
            positions[chrom].update(pos_range)
    return positions


def collapse_positions_to_ranges(chrom, pos_set):
    """Convert a set of positions back into sorted ranges"""
    if not pos_set:
        return []

    sorted_positions = sorted(pos_set)
    ranges = []
    start = prev = sorted_positions[0]

    for pos in sorted_positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            if start == prev:
                ranges.append(f"{chrom}:{start}")
            else:
                ranges.append(f"{chrom}:{start}-{prev}")
            start = prev = pos
    # flush last range
    if start == prev:
        ranges.append(f"{chrom}:{start}")
    else:
        ranges.append(f"{chrom}:{start}-{prev}")
    return ranges


def main():
    if len(sys.argv) < 3:
        print("Usage: python intersect_outliers.py <outlier_file1> <outlier_file2> [outlier_file3 ...]")
        sys.exit(1)

    # Parse first file
    common_positions = parse_outlier_file(sys.argv[1])

    # Intersect with each subsequent file
    for filename in sys.argv[2:]:
        print(f"Intersecting with {filename}...")
        file_positions = parse_outlier_file(filename)
        new_common = {}
        for chrom in common_positions:
            if chrom in file_positions:
                new_common[chrom] = common_positions[chrom].intersection(file_positions[chrom])
        common_positions = new_common  # keep only shared chroms/positions

    # Write output
    outfile = "shared_outliers.txt"
    with open(outfile, "w") as out:
        out.write("Shared high-depth positions across all files:\n\n")
        if not common_positions:
            out.write("None found.\n")
        else:
            for chrom in sorted(common_positions.keys()):
                ranges = collapse_positions_to_ranges(chrom, common_positions[chrom])
                for r in ranges:
                    out.write(r + "\n")

    print(f"Intersection complete. Results written to {outfile}")


if __name__ == "__main__":
    main()