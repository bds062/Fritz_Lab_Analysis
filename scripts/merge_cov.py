#!/usr/bin/env python3
"""
Find intersecting windows across filtered .cov.gz files.
Outputs windows where at least 50% of files have overlapping regions,
annotated with overlapping gene names from GFF file.
"""

import gzip
import sys
import re
from collections import defaultdict
from pathlib import Path

# Hardcoded GFF file path
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
        if not (window_end <= gene_start or window_start >= gene_end):
            overlapping.append(gene_name)
    
    return overlapping


def read_windows(file_path):
    """Read windows (chr, start, end) from a .cov.gz file."""
    windows = []
    with gzip.open(file_path, 'rt') as f:
        f.readline()  # Skip header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            chr_name = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            windows.append((chr_name, start, end))
    return windows


def find_overlaps(windows_by_file, min_files):
    """
    Find overlapping windows across files.
    Returns intersection coordinates where at least min_files have overlap.
    """
    # Collect all unique breakpoints per chromosome
    breakpoints = defaultdict(set)
    
    for windows in windows_by_file:
        for chr_name, start, end in windows:
            breakpoints[chr_name].add(start)
            breakpoints[chr_name].add(end)
    
    # For each chromosome, create intervals between breakpoints
    intersecting_windows = []
    
    for chr_name in sorted(breakpoints.keys()):
        sorted_points = sorted(breakpoints[chr_name])
        
        # Check each interval between consecutive breakpoints
        for i in range(len(sorted_points) - 1):
            interval_start = sorted_points[i]
            interval_end = sorted_points[i + 1]
            
            # Count how many files have a window overlapping this interval
            overlapping_files = 0
            for windows in windows_by_file:
                has_overlap = any(
                    chr_name == w_chr and 
                    w_start < interval_end and 
                    w_end > interval_start
                    for w_chr, w_start, w_end in windows
                )
                if has_overlap:
                    overlapping_files += 1
            
            # If at least min_files overlap, keep this interval
            if overlapping_files >= min_files:
                intersecting_windows.append((chr_name, interval_start, interval_end))
    
    return intersecting_windows


def merge_adjacent_windows(windows):
    """Merge adjacent windows on the same chromosome."""
    if not windows:
        return []
    
    merged = []
    current_chr, current_start, current_end = windows[0]
    
    for chr_name, start, end in windows[1:]:
        if chr_name == current_chr and start == current_end:
            # Adjacent window, extend the current one
            current_end = end
        else:
            # Non-adjacent, save current and start new
            merged.append((current_chr, current_start, current_end))
            current_chr, current_start, current_end = chr_name, start, end
    
    # Add the last window
    merged.append((current_chr, current_start, current_end))
    
    return merged


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py filtered_*.cov.gz")
        sys.exit(1)
    
    input_files = sys.argv[1:]
    num_files = len(input_files)
    min_files = max(1, int(num_files * 0.5 + 0.5))  # At least 50% (round up)
    
    print(f"Processing {num_files} file(s)...")
    print(f"Minimum files for overlap: {min_files} ({min_files/num_files*100:.0f}%)\n")
    
    # Parse GFF file
    print(f"Parsing GFF file: {GFF_FILE}...")
    genes_by_chrom = parse_gff_genes(GFF_FILE)
    total_genes = sum(len(genes) for genes in genes_by_chrom.values())
    print(f"  Found {total_genes} genes across {len(genes_by_chrom)} chromosomes\n")
    
    # Read all windows from all files
    windows_by_file = []
    for file_path in input_files:
        print(f"Reading {Path(file_path).name}...")
        windows = read_windows(file_path)
        windows_by_file.append(windows)
        print(f"  {len(windows)} windows")
    
    print("\nFinding intersecting windows...")
    intersecting = find_overlaps(windows_by_file, min_files)
    
    print("Merging adjacent windows...")
    merged = merge_adjacent_windows(intersecting)
    
    # Annotate with genes and filter for windows with genes
    print("Annotating with gene names...")
    annotated_windows = []
    for chr_name, start, end in merged:
        overlapping_genes = find_overlapping_genes(chr_name, start, end, genes_by_chrom)
        # Only include windows that overlap with at least one gene
        for gene_name in overlapping_genes:
            annotated_windows.append((chr_name, start, end, gene_name))
    
    # Write output
    output_file = "merged_3.0.tsv"
    print(f"\nWriting output to {output_file}...")
    with open(output_file, 'w') as f:
        f.write("chr\tstart\tend\tgene_name\n")
        for chr_name, start, end, gene_name in annotated_windows:
            f.write(f"{chr_name}\t{start}\t{end}\t{gene_name}\n")
    
    print(f"Done! Found {len(annotated_windows)} windows with associated genes.")
    print(f"({len(merged) - len(set((w[0], w[1], w[2]) for w in annotated_windows))} windows without genes were excluded)")


if __name__ == "__main__":
    main()