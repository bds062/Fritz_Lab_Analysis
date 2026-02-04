import pandas as pd
import numpy as np
import re
from collections import defaultdict

# Hardcoded input filenames
input_file = "mergeCNVR"
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

# Read the TSV file
df = pd.read_csv(input_file, sep='\t')

# Calculate mean and standard deviation of the 'average' column
mean_avg = df['average'].mean()
std_avg = df['average'].std()

# Calculate threshold
# threshold = mean_avg + 3 * std_avg
threshold = 3

# Print statistics to stdout
print(f"Mean of 'average' column: {mean_avg:.4f}")
print(f"Standard deviation of 'average' column: {std_avg:.4f}")
print(f"Threshold (3): {threshold:.4f}")
print(f"Number of rows above threshold: {(df['average'] > threshold).sum()}")

# Filter rows above threshold
filtered_df = df[df['average'] > threshold]

# Parse GFF file
print("\nParsing GFF file...")
genes_by_chrom = parse_gff_genes(GFF_FILE)

# Find overlapping genes for each filtered row
rows_with_genes = []
for idx, row in filtered_df.iterrows():
    chrom = row['chr']
    window_start = row['start']
    window_end = row['end']
    
    overlapping_genes = find_overlapping_genes(chrom, window_start, window_end, genes_by_chrom)
    
    # Create one row per overlapping gene
    for gene_name in overlapping_genes:
        row_dict = row.to_dict()
        row_dict['gene_name'] = gene_name
        rows_with_genes.append(row_dict)

# Create output dataframe
if rows_with_genes:
    output_df = pd.DataFrame(rows_with_genes)
    print(f"Number of rows with overlapping genes: {len(output_df)}")
else:
    # Create empty dataframe with correct columns if no rows match
    output_df = filtered_df.copy()
    output_df['gene_name'] = []
    output_df = output_df[0:0]  # Empty dataframe with columns
    print("Number of rows with overlapping genes: 0")

# Write filtered results to new file
output_file = f"filtered_{input_file}"
output_df.to_csv(output_file, sep='\t', index=False)

print(f"\nFiltered results written to: {output_file}")