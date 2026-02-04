#!/usr/bin/env python3

import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
import sys
import os

# Configuration
CHROMOSOME = "CM040209.1"
BAM_DIR = "."
MERGE_CNVR_FILE = "./mergeCNVR"

# Get current directory name for output files
current_dir = os.path.basename(os.getcwd())
OUTPUT_TSV = f"{current_dir}_coverage_cnv_{CHROMOSOME}.tsv"
OUTPUT_PLOT = f"{current_dir}_coverage_cnv_{CHROMOSOME}.png"

def get_individual_from_filename(filename):
    """Extract SRR identifier from BAM filename."""
    match = re.search(r'RG_trimmed_(\w+)\.bam', filename)
    if match:
        return match.group(1)
    return None

def get_coverage(bam_file, chromosome=None):
    """Get mean depth coverage averaged across all main chromosomes using samtools."""
    try:
        # Run samtools coverage without chromosome filter to get all contigs
        cmd = f"samtools coverage {bam_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"ERROR: samtools failed for {bam_file}: {result.stderr}", file=sys.stderr)
            return None
        
        lines = result.stdout.strip().split('\n')
        
        # Skip header line and parse data
        data_lines = [line for line in lines if not line.startswith('#')]
        
        if len(data_lines) == 0:
            print(f"ERROR: No coverage data for {bam_file}", file=sys.stderr)
            return None
        
        # Filter for main chromosomes (CM0*) and extract meandepth values
        meandepths = []
        for line in data_lines:
            fields = line.split('\t')
            if len(fields) >= 7:
                contig_name = fields[0]
                # Only include main chromosomes (CM0*.1)
                if contig_name.startswith('CM0'):
                    try:
                        meandepth = float(fields[6])
                        meandepths.append(meandepth)
                    except ValueError:
                        continue
        
        if len(meandepths) == 0:
            print(f"ERROR: No CM0* chromosome data found for {bam_file}", file=sys.stderr)
            return None
        
        # Return the average meandepth across all main chromosomes
        average_coverage = sum(meandepths) / len(meandepths)
        return average_coverage
        
    except Exception as e:
        print(f"ERROR: Failed to get coverage for {bam_file}: {e}", file=sys.stderr)
        return None

def get_cnv_data(merge_cnvr_file, chromosome):
    """Parse mergeCNVR file and extract CNV data for the specified chromosome."""
    try:
        df = pd.read_csv(merge_cnvr_file, sep='\t')
        
        # Filter for the specified chromosome
        chrom_data = df[df['chr'] == chromosome]
        
        if len(chrom_data) == 0:
            print(f"ERROR: No CNV data found for chromosome {chromosome} in {merge_cnvr_file}", file=sys.stderr)
            return None
        
        return chrom_data
        
    except Exception as e:
        print(f"ERROR: Failed to read mergeCNVR file: {e}", file=sys.stderr)
        return None

def calculate_average_cnv(cnv_data, individual):
    """Calculate average CNV for a specific individual across all windows."""
    if cnv_data is None:
        return None
    
    if individual not in cnv_data.columns:
        print(f"ERROR: Individual {individual} not found in mergeCNVR file", file=sys.stderr)
        return None
    
    try:
        avg_cnv = cnv_data[individual].mean()
        return avg_cnv
    except Exception as e:
        print(f"ERROR: Failed to calculate average CNV for {individual}: {e}", file=sys.stderr)
        return None

def main():
    # Get the current directory name (last folder in path)
    location = os.path.basename(os.getcwd())
    
    # Find all BAM files
    bam_files = glob.glob(f"{BAM_DIR}/RG_trimmed_*.bam")
    
    if len(bam_files) == 0:
        print(f"ERROR: No BAM files found matching pattern 'RG_trimmed_*.bam' in {BAM_DIR}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(bam_files)} BAM files")
    
    # Load CNV data once
    cnv_data = get_cnv_data(MERGE_CNVR_FILE, CHROMOSOME)
    
    # Process each BAM file
    results = []
    
    for bam_file in bam_files:
        individual = get_individual_from_filename(bam_file)
        
        if individual is None:
            print(f"ERROR: Could not extract individual ID from filename {bam_file}", file=sys.stderr)
            continue
        
        print(f"Processing {individual}...")
        
        # Get coverage
        coverage = get_coverage(bam_file, CHROMOSOME)
        if coverage is None:
            coverage = "NA"
        
        # Get average CNV
        avg_cnv = calculate_average_cnv(cnv_data, individual)
        if avg_cnv is None:
            avg_cnv = "NA"
        
        results.append({
            'location': location,
            'individual': individual,
            'coverage': coverage,
            'CNV': avg_cnv
        })
    
    # Create DataFrame and save to TSV
    df_results = pd.DataFrame(results)
    df_results.to_csv(OUTPUT_TSV, sep='\t', index=False)
    print(f"\nResults saved to {OUTPUT_TSV}")
    
    # Create scatterplot (only for valid numeric values)
    df_plot = df_results[(df_results['coverage'] != "NA") & (df_results['CNV'] != "NA")].copy()
    df_plot['coverage'] = pd.to_numeric(df_plot['coverage'])
    df_plot['CNV'] = pd.to_numeric(df_plot['CNV'])
    
    if len(df_plot) == 0:
        print("ERROR: No valid data points to plot", file=sys.stderr)
        sys.exit(1)
    
    plt.figure(figsize=(10, 6))
    plt.scatter(df_plot['coverage'], df_plot['CNV'])
    plt.xlabel('Coverage (Mean Depth)')
    plt.ylabel('Average CNV')
    plt.title(f'Coverage vs CNV for {CHROMOSOME}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"Plot saved to {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()