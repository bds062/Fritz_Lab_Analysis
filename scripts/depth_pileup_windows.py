#!/usr/bin/env python3
"""
Identify windows of high depth coverage from filtered depth file.
Windows must be at least 1000 bp long with at least 80% high coverage positions.
Uses fast heuristic approach for large datasets.
"""

import sys
import re
from collections import defaultdict

# GFF file path
GFF_FILE = "/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/Liftoff_Analysis/test1.gff"

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
            
            # Extract Note field
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
        # Check if window and gene overlap
        if not (window_end < gene_start or window_start > gene_end):
            overlapping.append(gene_name)
    
    return overlapping


def find_windows_heuristic(high_cov_positions, min_length=1000, min_percent=0.80):
    """
    Fast heuristic algorithm: only start windows at high coverage positions.
    
    Args:
        high_cov_positions: set of positions with high coverage
        min_length: minimum window length in bp
        min_percent: minimum percentage of high coverage positions in window
    
    Returns:
        List of tuples: (start, end, high_cov_count, percent)
    """
    if not high_cov_positions:
        return []
    
    positions_list = sorted(high_cov_positions)
    
    if len(positions_list) < min_length * min_percent:
        return []  # Not enough positions to form any valid window
    
    windows = []
    i = 0
    n = len(positions_list)
    
    while i < n:
        start_pos = positions_list[i]
        
        # Quick density check: are there enough nearby positions?
        # Count positions within next min_length bp
        j = i
        while j < n and positions_list[j] <= start_pos + min_length - 1:
            j += 1
        
        positions_in_min_window = j - i
        
        # If minimum window doesn't have enough coverage, skip ahead
        if positions_in_min_window < min_length * min_percent:
            i += 1
            continue
        
        # Extend window as far as possible while maintaining threshold
        # Use positions_list indices for fast counting
        best_end_idx = j - 1
        best_end_pos = positions_list[best_end_idx]
        
        # Try extending further
        k = j
        while k < n:
            test_end_pos = positions_list[k]
            window_length = test_end_pos - start_pos + 1
            high_cov_count = k - i + 1
            
            if high_cov_count / window_length >= min_percent:
                best_end_idx = k
                best_end_pos = test_end_pos
                k += 1
            else:
                # Try a few more positions in case we can recover
                # This handles cases where a gap appears but we can still maintain threshold
                lookahead = min(50, n - k)  # Look ahead up to 50 positions
                found_better = False
                
                for m in range(k + 1, k + lookahead):
                    if m >= n:
                        break
                    test_end_pos_2 = positions_list[m]
                    window_length_2 = test_end_pos_2 - start_pos + 1
                    high_cov_count_2 = m - i + 1
                    
                    if high_cov_count_2 / window_length_2 >= min_percent:
                        best_end_idx = m
                        best_end_pos = test_end_pos_2
                        k = m + 1
                        found_better = True
                        break
                
                if not found_better:
                    break
        
        # Check if we have a valid window
        final_length = best_end_pos - start_pos + 1
        final_count = best_end_idx - i + 1
        
        if final_length >= min_length and final_count / final_length >= min_percent:
            windows.append((start_pos, best_end_pos, final_count, final_count / final_length))
            
            # Skip to the position after this window to avoid overlaps
            i = best_end_idx + 1
        else:
            i += 1
    
    return windows

def main():
    if len(sys.argv) != 2:
        print("Usage: python window_finder.py <filtered_depth_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = f"windows_{input_file}"
    
    print(f"Reading filtered depth file: {input_file}")
    
    # Read high coverage positions by chromosome
    chrom_positions = defaultdict(list)
    total_positions = 0
    
    with open(input_file, 'r') as f:
        header = f.readline().strip()
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) != 4:
                continue
            
            chromosome = parts[0]
            try:
                position = int(parts[1])
                chrom_positions[chromosome].append(position)
                total_positions += 1
            except ValueError:
                continue
    
    # Sort positions for each chromosome
    for chrom in chrom_positions:
        chrom_positions[chrom].sort()
    
    print(f"Total high coverage positions: {total_positions:,}")
    print(f"Found on {len(chrom_positions)} chromosome(s)")
    
    # Parse GFF file for gene annotations
    print(f"\nParsing GFF file: {GFF_FILE}")
    genes_by_chrom = parse_gff_genes(GFF_FILE)
    total_genes = sum(len(genes) for genes in genes_by_chrom.values())
    print(f"Loaded {total_genes:,} genes from GFF file")
    
    # Find windows for each chromosome
    print("\nFinding windows (min 1000 bp, min 80% high coverage)...")
    
    all_windows = []
    
    for chrom in sorted(chrom_positions.keys()):
        num_positions = len(chrom_positions[chrom])
        print(f"  Processing {chrom} ({num_positions:,} positions)... ", end='', flush=True)
        
        # Convert to set for the function
        windows = find_windows_heuristic(set(chrom_positions[chrom]))
        
        for start, end, high_cov_count, percent_high in windows:
            window_length = end - start + 1
            
            # Find overlapping genes
            overlapping_genes = find_overlapping_genes(chrom, start, end, genes_by_chrom)
            
            # Only report windows that have an associated gene
            if overlapping_genes:
                # Create a separate row for each overlapping gene
                for gene_name in overlapping_genes:
                    all_windows.append({
                        'chromosome': chrom,
                        'window_start': start,
                        'window_end': end,
                        'window_length': window_length,
                        'high_coverage_count': high_cov_count,
                        'percent_high_coverage': percent_high,
                        'gene_name': gene_name
                    })
        
        print(f"{len(windows)} windows found")
    
    # Write output
    print(f"\nTotal windows found: {len(all_windows)}")
    print(f"Writing output to: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("chromosome\twindow_start\twindow_end\twindow_length\thigh_coverage_count\tpercent_high_coverage\tgene_name\n")
        
        for window in all_windows:
            f.write(f"{window['chromosome']}\t{window['window_start']}\t{window['window_end']}\t"
                   f"{window['window_length']}\t{window['high_coverage_count']}\t"
                   f"{window['percent_high_coverage']:.4f}\t{window['gene_name']}\n")
    
    print("Done!")

if __name__ == "__main__":
    main()