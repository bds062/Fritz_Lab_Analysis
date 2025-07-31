import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import defaultdict
import os

# === PARAMETERS ===
bin_size = 1_000_000  # 1 Mb
output_plot = "contaminated_alignment_density.png"
highlight_threshold = 150  # updated threshold for avg count

# === FUNCTION TO COUNT READS BY POSITION BIN ===
def get_read_counts_by_bin(sam_path, bin_size):
    binned_counts = defaultdict(int)
    try:
        with pysam.AlignmentFile(sam_path, "r") as samfile:
            for read in samfile.fetch(until_eof=True):
                if not read.is_unmapped:
                    chrom = read.reference_name
                    pos = read.reference_start
                    bin_id = (chrom, pos // bin_size)
                    binned_counts[bin_id] += 1
    except Exception as e:
        print(f"Error reading {sam_path}: {e}")
    return binned_counts

# === COLLECT ALL FILES ===
sam_files = [f for f in os.listdir() if f.startswith("contaminated_S") and f.endswith(".sam")]
print(f"Found {len(sam_files)} contaminated SAM files.")

# === AGGREGATE BIN COUNTS ===
all_bin_counts = defaultdict(list)

for sam_file in sam_files:
    print(f"Processing {sam_file}...")
    bin_counts = get_read_counts_by_bin(sam_file, bin_size)
    for bin_id, count in bin_counts.items():
        all_bin_counts[bin_id].append(count)

# === SORT BINS ===
sorted_bins = sorted(all_bin_counts.keys(), key=lambda x: (x[0], x[1]))
x_labels = [f"{chrom}:{bin_id * bin_size // 1_000_000}Mb" for chrom, bin_id in sorted_bins]
x_values = np.arange(len(sorted_bins))
y_values = [np.mean(all_bin_counts[bin_id]) for bin_id in sorted_bins]

# === HIGHLIGHTED BINS ===
highlight_x = []
highlight_y = []
highlight_err = []
highlight_stats = []

print("\n=== BINS ABOVE 150 AVERAGE MAPPED READS ===")
for i, bin_id in enumerate(sorted_bins):
    avg = np.mean(all_bin_counts[bin_id])
    if avg > highlight_threshold:
        std = np.std(all_bin_counts[bin_id])
        chrom, bin_num = bin_id
        start_bp = bin_num * bin_size
        end_bp = start_bp + bin_size
        print(f"{chrom}:{start_bp}-{end_bp}\tAvg Count: {avg:.2f}\tStd Dev: {std:.2f}")
        highlight_x.append(x_values[i])
        highlight_y.append(avg)
        highlight_err.append(std)

# === PLOT ===
plt.figure(figsize=(18, 6))
plt.bar(x_values, y_values, width=1.0, color='darkred', label="Avg Mapped Reads per File")

# Add error bars for highlighted bins
if highlight_x:
    plt.errorbar(highlight_x, highlight_y, yerr=highlight_err, fmt='o', ecolor='black', capsize=3, label=">150 (mean ± std)")

plt.xticks(
    x_values[::max(len(x_values) // 10, 1)],
    x_labels[::max(len(x_labels) // 10, 1)],
    rotation=45
)
plt.ylabel("Avg Mapped Read Count per File")
plt.xlabel("Genomic Position (1Mb bins)")
plt.title("Read Density Across Genome - Contaminated Samples")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.close()

print(f"\nSaved plot to {output_plot}")