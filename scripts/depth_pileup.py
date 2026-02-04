#!/usr/bin/env python3
"""
Calculate mean and standard deviation of depth across multiple samtools depth files.
Handles missing positions (zero coverage) correctly.
Usage: python depthPileup.py *_depth.tsv
"""

import sys
import numpy as np
from collections import defaultdict

def read_next_position(file_handle):
    """Read next line and return (chrom, pos, depth) or None if EOF."""
    line = file_handle.readline()
    if not line:
        return None
    parts = line.strip().split('\t')
    if len(parts) != 3:
        return None
    return (parts[0], int(parts[1]), int(parts[2]))

def main():
    if len(sys.argv) < 2:
        print("Usage: python depthPileup.py *_depth.tsv", file=sys.stderr)
        sys.exit(1)
    
    files = sys.argv[1:]
    n_files = len(files)
    
    # Open all files
    try:
        file_handles = [open(fname, 'r') for fname in files]
    except Exception as e:
        print(f"Error opening files: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Read first position from each file
    current_positions = [read_next_position(f) for f in file_handles]
    
    # Output header
    print("chromosome\tposition\tmean_depth\tstdev_depth")
    
    try:
        while any(pos is not None for pos in current_positions):
            # Find the minimum position across all files
            valid_positions = [pos for pos in current_positions if pos is not None]
            if not valid_positions:
                break
            
            min_chrom = min(pos[0] for pos in valid_positions)
            min_pos = min(pos[1] for pos in valid_positions if pos[0] == min_chrom)
            
            # Collect depths for this position (0 if file doesn't have this position)
            depths = []
            for i, pos in enumerate(current_positions):
                if pos is not None and pos[0] == min_chrom and pos[1] == min_pos:
                    depths.append(pos[2])
                    # Read next position from this file
                    current_positions[i] = read_next_position(file_handles[i])
                else:
                    # This file doesn't have this position (zero coverage)
                    depths.append(0)
            
            # Calculate statistics
            depths_arr = np.array(depths)
            mean_depth = np.mean(depths_arr)
            stdev_depth = np.std(depths_arr, ddof=1) if n_files > 1 else 0.0
            
            print(f"{min_chrom}\t{min_pos}\t{mean_depth:.2f}\t{stdev_depth:.2f}")
    
    finally:
        # Close all file handles
        for f in file_handles:
            f.close()

if __name__ == "__main__":
    main()