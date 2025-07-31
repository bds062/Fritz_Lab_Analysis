import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
import os
from collections import defaultdict
import csv 

# CONFIGURATION
user_chromosome = "chr9" # Set to specific chromosome (e.g., "chr1") or None for all chromosomes
bin_size = 10_000
low_depth_threshold = 5

chromref_path = "/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/chromref.tsv"

def load_chromosome_map(ref_file):
    """Load mapping from friendly chromosome name (e.g., chr2) to GenBank accession (e.g., CM040204.1)"""
    chrom_map = {}
    with open(ref_file, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            human_name = f"chr{row['Chromosome name']}"
            accession = row['GenBank seq accession']
            chrom_map[human_name] = accession
    return chrom_map

def process_file(input_file, target_chrom=None):
    """Process a single depth file and return binned data"""
    bin_sums = defaultdict(int)
    bin_counts = defaultdict(int)
    
    print(f"Processing {input_file}...")
    
    with open(input_file) as f:
        for line in f:
            chrom, pos_str, depth_str = line.strip().split()
            
            if target_chrom and chrom != target_chrom:
                continue
                
            pos = int(pos_str)
            depth = int(depth_str)
            bin_id = (chrom, pos // bin_size)
            bin_sums[bin_id] += depth
            bin_counts[bin_id] += 1
    
    # Compute average depths
    binned_avg = {bin_id: bin_sums[bin_id] / bin_counts[bin_id] for bin_id in bin_sums}
    return binned_avg

def create_plot(binned_avg, output_filename, plot_chrom_name=None):
    """Create and save the depth plot"""
    if not binned_avg:
        print(f"No data found for {plot_chrom_name if plot_chrom_name else 'any chromosome'}")
        return
    
    sorted_bins = sorted(binned_avg.keys(), key=lambda x: (x[0], x[1]))
    
    x_labels = [f"{chrom}:{bin_num}Mb" for chrom, bin_num in sorted_bins]
    x = np.arange(len(sorted_bins))
    y = [binned_avg[bin_id] for bin_id in sorted_bins]
    colors = ['crimson' if val < low_depth_threshold else 'steelblue' for val in y]
    
    # Plot
    plt.figure(figsize=(18, 6))
    plt.bar(x, y, color=colors)
    plt.xticks(x[::max(len(x)//20, 1)], x_labels[::max(len(x_labels)//20, 1)], rotation=45)
    plt.xlabel(f"Genomic Bin ({bin_size}bp)")
    plt.ylabel("Avg Read Depth")
    
    if plot_chrom_name:
        plt.title(f"Binned Read Depth - {plot_chrom_name}")
    else:
        plt.title("Binned Read Depth Across Genome")
    
    plt.grid(True)
    
    # Legend
    import matplotlib.patches as mpatches
    normal_patch = mpatches.Patch(color='steelblue', label=f'> {low_depth_threshold}x Depth')
    low_patch = mpatches.Patch(color='crimson', label=f'< {low_depth_threshold}x Depth')
    plt.legend(handles=[normal_patch, low_patch])
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()
    
    print(f"Plot saved to '{output_filename}'")

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <file_pattern>")
        sys.exit(1)
    
    # Load chromosome map
    chrom_map = load_chromosome_map(chromref_path)
    
    # Translate human-readable chromosome name to actual accession
    if user_chromosome in chrom_map:
        target_chromosome = chrom_map[user_chromosome]
    else:
        print(f"Error: Chromosome name '{user_chromosome}' not found in {chromref_path}")
        sys.exit(1)

    # Get all input files
    input_files = []
    for arg in sys.argv[1:]:
        if '*' in arg:
            matched_files = glob.glob(arg)
            input_files.extend(matched_files)
        elif os.path.exists(arg):
            input_files.append(arg)
        else:
            print(f"File not found: {arg}")
    
    if not input_files:
        print("No valid input files found!")
        sys.exit(1)
    
    print(f"Found {len(input_files)} files to process")
    print(f"Target chromosome: {user_chromosome} → {target_chromosome}")
    print(f"Bin size: {bin_size:,} bp")
    
    for input_file in input_files:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        if base_name.endswith('_depth'):
            base_name = base_name[:-6]
        
        binned_data = process_file(input_file, target_chromosome)
        
        output_filename = f"{base_name}_{user_chromosome}_{bin_size//1000}kbp_depth.png"
        create_plot(binned_data, output_filename, user_chromosome)

if __name__ == "__main__":
    main()