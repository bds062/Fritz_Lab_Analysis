#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from scipy import stats

OUTPUT_PLOT = "coverage_cnv_CM040210.1.png"

def extract_directory_name(tsv_filename):
    """Extract directory name from TSV filename (e.g., 'Ja_coverage_cnv.tsv' -> 'Ja')."""
    basename = os.path.basename(tsv_filename)
    # Remove '_coverage_cnv.tsv' suffix
    dir_name = basename.replace('_coverage_cnv.tsv', '')
    return dir_name

def load_tsv_data(tsv_file):
    """Load TSV file and return valid numeric data points."""
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Filter out NA values and convert to numeric
        df_valid = df[(df['coverage'] != "NA") & (df['CNV'] != "NA")].copy()
        df_valid['coverage'] = pd.to_numeric(df_valid['coverage'])
        df_valid['CNV'] = pd.to_numeric(df_valid['CNV'])
        
        if len(df_valid) == 0:
            print(f"ERROR: No valid data points in {tsv_file}", file=sys.stderr)
            return None
        
        return df_valid
        
    except FileNotFoundError:
        print(f"ERROR: File not found: {tsv_file}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"ERROR: Failed to load {tsv_file}: {e}", file=sys.stderr)
        return None

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <tsv_file1> <tsv_file2> ...", file=sys.stderr)
        sys.exit(1)
    
    tsv_files = sys.argv[1:]
    
    # Load all data
    all_data = []
    colors = plt.cm.tab10(np.linspace(0, 1, len(tsv_files)))
    
    for idx, tsv_file in enumerate(tsv_files):
        print(f"Loading {tsv_file}...")
        df = load_tsv_data(tsv_file)
        
        if df is not None:
            dir_name = extract_directory_name(tsv_file)
            all_data.append({
                'data': df,
                'dir_name': dir_name,
                'color': colors[idx]
            })
    
    if len(all_data) == 0:
        print("ERROR: No valid data to plot", file=sys.stderr)
        sys.exit(1)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Collect all points for line of best fit
    all_coverage = []
    all_cnv = []
    
    # Plot each dataset with different color
    for dataset in all_data:
        df = dataset['data']
        dir_name = dataset['dir_name']
        color = dataset['color']
        
        ax.scatter(df['coverage'], df['CNV'], 
                  label=dir_name, color=color, alpha=0.6, s=50)
        
        all_coverage.extend(df['coverage'].tolist())
        all_cnv.extend(df['CNV'].tolist())
    
    # Calculate line of best fit for all data
    all_coverage = np.array(all_coverage)
    all_cnv = np.array(all_cnv)
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(all_coverage, all_cnv)
    r_squared = r_value ** 2
    
    # Plot line of best fit
    x_line = np.array([all_coverage.min(), all_coverage.max()])
    y_line = slope * x_line + intercept
    
    ax.plot(x_line, y_line, 'k--', linewidth=2, 
            label=f'Best Fit: y = {slope:.4f}x + {intercept:.4f}\n$R^2$ = {r_squared:.4f}')
    
    # Labels and formatting
    ax.set_xlabel('Coverage (Mean Depth) on Chr 8', fontsize=12)
    ax.set_ylabel('Average CNV on Chr 8', fontsize=12)
    ax.set_title('Coverage vs CNV', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"\nCombined plot saved to {OUTPUT_PLOT}")
    print(f"Line of best fit: y = {slope:.4f}x + {intercept:.4f}")
    print(f"R² = {r_squared:.4f}")

if __name__ == "__main__":
    main()