#!/usr/bin/env python3
"""
Filter depth file to find positions with mean_depth > genome_wide_mean + 6*genome_wide_stdev
Uses streaming statistics for memory efficiency with large files.
"""

import sys
import math

def welford_update(count, mean, m2, new_value):
    """
    Update running statistics using Welford's online algorithm.
    Returns updated count, mean, and M2 (sum of squared differences).
    """
    count += 1
    delta = new_value - mean
    mean += delta / count
    delta2 = new_value - mean
    m2 += delta * delta2
    return count, mean, m2

def main():
    if len(sys.argv) != 2:
        print("Usage: python depth_filter.py <input_depth_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = f"filtered_1.5xmean_{input_file}"
    
    print(f"Reading input file: {input_file}")
    print("Pass 1: Calculating genome-wide statistics...")
    
    # First pass: calculate mean and stdev using Welford's algorithm
    count = 0
    mean = 0.0
    m2 = 0.0
    
    with open(input_file, 'r') as f:
        header = f.readline().strip()
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) != 4:
                continue
            
            try:
                mean_depth = float(parts[2])
                count, mean, m2 = welford_update(count, mean, m2, mean_depth)
            except ValueError:
                continue
    
    # Calculate standard deviation
    if count < 2:
        print("Error: Not enough valid data points")
        sys.exit(1)
    
    variance = m2 / (count - 1)
    stdev = math.sqrt(variance)
    
    print(f"\nGenome-wide statistics:")
    print(f"  Total positions: {count:,}")
    print(f"  Mean depth: {mean:.2f}")
    print(f"  Stdev depth: {stdev:.2f}")
    
    # Calculate threshold
    # threshold = mean + 6 * stdev
    threshold = mean*1.5
    print(f"  Threshold (1.5*mean): {threshold:.2f}")
    
    print("\nPass 2: Filtering positions above threshold...")
    
    # Second pass: filter and write output
    filtered_count = 0
    
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = f_in.readline().strip()
        f_out.write(header + '\n')
        
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) != 4:
                continue
            
            try:
                mean_depth = float(parts[2])
                if mean_depth > threshold:
                    f_out.write(line + '\n')
                    filtered_count += 1
            except ValueError:
                continue
    
    print(f"\nFiltered positions (above threshold): {filtered_count:,}")
    print(f"Percentage of positions above threshold: {(filtered_count/count)*100:.4f}%")
    print(f"\nOutput written to: {output_file}")

if __name__ == "__main__":
    main()