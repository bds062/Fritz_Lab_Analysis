#!/usr/bin/env python3
"""
Filter .cov.gz files based on copy number (CN) threshold.
Keeps only rows where CN > mean + 6*stdev.
"""

import gzip
import sys
import statistics
from pathlib import Path


def process_file(input_path):
    """Process a single .cov.gz file and create filtered output."""
    input_path = Path(input_path)
    output_path = input_path.parent / f"filtered_3.0_{input_path.name}"
    
    print(f"Processing {input_path.name}...")
    
    # First pass: read CN values to calculate mean and stdev
    cn_values = []
    with gzip.open(input_path, 'rt') as f:
        header = f.readline()  # Skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            cn_values.append(float(fields[-1]))
    
    # Calculate threshold
    mean_cn = statistics.mean(cn_values)
    stdev_cn = statistics.stdev(cn_values)
    # threshold = mean_cn + 6 * stdev_cn
    threshold = 3    
    print(f"  Mean: {mean_cn:.4f}, Stdev: {stdev_cn:.4f}, Threshold: {threshold:.4f}")
    
    # Second pass: filter and write output
    rows_kept = 0
    with gzip.open(input_path, 'rt') as f_in, \
         gzip.open(output_path, 'wt') as f_out:
        # Write header
        header = f_in.readline()
        f_out.write(header)
        
        # Filter and write rows
        for line in f_in:
            fields = line.rstrip('\n').split('\t')
            cn_value = float(fields[-1])
            if cn_value > threshold:
                f_out.write(line)
                rows_kept += 1
    
    print(f"  Kept {rows_kept} / {len(cn_values)} rows ({100*rows_kept/len(cn_values):.2f}%)")
    print(f"  Output: {output_path.name}\n")


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py *.cov.gz")
        sys.exit(1)
    
    input_files = sys.argv[1:]
    
    print(f"Processing {len(input_files)} file(s)...\n")
    
    for input_file in input_files:
        try:
            process_file(input_file)
        except Exception as e:
            print(f"Error processing {input_file}: {e}\n")
            continue
    
    print("Done!")


if __name__ == "__main__":
    main()