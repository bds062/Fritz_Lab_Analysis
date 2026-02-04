import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ============ USER CONFIGURATION ============
CHROMOSOME = "CM040209.1"  # Change this to plot different chromosomes
MERGECNVR_FILE = "mergeCNVR"
DEPTH_FILE = "merged_depth.tsv"
# ============================================

# Read the CNV data
cnv_data = pd.read_csv(MERGECNVR_FILE, sep='\t')

# Filter for the specified chromosome
cnv_data = cnv_data[cnv_data['chr'] == CHROMOSOME].copy()

# Calculate the middle position of each CNV window
cnv_data['mid_position'] = (cnv_data['start'] + cnv_data['end']) / 2

# Read the depth data efficiently (only the chromosome we need)
# Since the file is sorted by chromosome, we can stop reading once we've passed our chromosome
print(f"Reading depth data for {CHROMOSOME}...")
depth_chunks = []
chunksize = 100000  # Read 100k rows at a time

for chunk in pd.read_csv(DEPTH_FILE, sep='\t', chunksize=chunksize):
    # Filter for our chromosome
    chr_data = chunk[chunk['chromosome'] == CHROMOSOME]
    
    if len(chr_data) > 0:
        depth_chunks.append(chr_data)
    
    # If we've found our chromosome and this chunk doesn't contain it,
    # we've passed it (since file is sorted), so stop reading
    if len(depth_chunks) > 0 and len(chr_data) == 0:
        print(f"Finished reading {CHROMOSOME} data")
        break

# Combine all chunks for our chromosome
if depth_chunks:
    depth_data = pd.concat(depth_chunks, ignore_index=True)
    print(f"Loaded {len(depth_data)} depth measurements for {CHROMOSOME}")
else:
    print(f"No data found for {CHROMOSOME}")
    depth_data = pd.DataFrame(columns=['chromosome', 'position', 'mean_depth', 'stdev_depth'])

# Calculate average depth for each CNV window
averaged_depths = []
window_positions = []

for idx, row in cnv_data.iterrows():
    start = row['start']
    end = row['end']
    mid = row['mid_position']
    
    # Get all depth values within this window
    window_depths = depth_data[(depth_data['position'] >= start) & 
                               (depth_data['position'] <= end)]
    
    if len(window_depths) > 0:
        avg_depth = window_depths['mean_depth'].mean()
        averaged_depths.append(avg_depth)
        window_positions.append(mid)
    else:
        # If no depth data in this window, append NaN
        averaged_depths.append(np.nan)
        window_positions.append(mid)

# Add averaged depths to the dataframe for correlation analysis
cnv_data['averaged_depth'] = averaged_depths

# Calculate correlation between CNV and coverage depth
# Remove NaN values for correlation calculation
valid_data = cnv_data[['average', 'averaged_depth']].dropna()

if len(valid_data) > 1:
    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(valid_data['average'], valid_data['averaged_depth'])
    
    # Spearman correlation (rank-based, more robust to outliers)
    spearman_r, spearman_p = stats.spearmanr(valid_data['average'], valid_data['averaged_depth'])
    
    print(f"\nCorrelation Analysis:")
    print(f"Pearson correlation: r = {pearson_r:.4f}, p-value = {pearson_p:.4e}")
    print(f"Spearman correlation: ρ = {spearman_r:.4f}, p-value = {spearman_p:.4e}")
else:
    pearson_r, pearson_p = np.nan, np.nan
    spearman_r, spearman_p = np.nan, np.nan
    print("\nNot enough valid data points for correlation analysis")

# Create the plot with dual y-axes
fig, ax1 = plt.subplots(figsize=(14, 6))

# Plot CNV data on the left y-axis
color1 = 'tab:blue'
ax1.set_xlabel('Genomic Position (bp)', fontsize=12)
ax1.set_ylabel('Average CNV', color=color1, fontsize=12)
ax1.plot(cnv_data['mid_position'], cnv_data['average'], 
         color=color1, linewidth=1.5, label='CNV', marker='o', markersize=3)
ax1.tick_params(axis='y', labelcolor=color1)
ax1.grid(True, alpha=0.3)

# Create second y-axis for coverage depth
ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel('Mean Coverage Depth', color=color2, fontsize=12)
ax2.plot(window_positions, averaged_depths, 
         color=color2, linewidth=1.5, label='Coverage Depth', marker='s', markersize=3)
ax2.tick_params(axis='y', labelcolor=color2)

# Add title
plt.title(f'CNV and Coverage Depth across {CHROMOSOME}', fontsize=14, fontweight='bold')

# Add vertical line at specific position for CM040209.1
if CHROMOSOME == "CM040209.1":
    ax1.axvline(x=8235200, color='green', linestyle='--', linewidth=2, alpha=0.7, label='8,235,200 bp')

# Add correlation annotation to the plot
if not np.isnan(pearson_r):
    textstr = f'Pearson r = {pearson_r:.3f} (p = {pearson_p:.2e})\nSpearman ρ = {spearman_r:.3f} (p = {spearman_p:.2e})'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=props)

# Add legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.tight_layout()
plt.savefig(f'{CHROMOSOME}_cnv_coverage_plot.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Plot saved as {CHROMOSOME}_cnv_coverage_plot.png")
print(f"Number of CNV windows plotted: {len(cnv_data)}")
print(f"Number of depth measurements: {len(depth_data)}")