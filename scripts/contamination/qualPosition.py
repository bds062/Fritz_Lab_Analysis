# mapq_plotter.py
import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import defaultdict
import sys

# === READ SRR ID FROM COMMAND LINE ===
if len(sys.argv) != 2:
    print("Usage: python3 mapq_plotter.py SRR_ID")
    sys.exit(1)

srr_id = sys.argv[1]

# === FILE NAMES BASED ON ID ===
clean_file = f"{srr_id}.sam"
contaminated_file = f"trimmed_{srr_id}.sam"
output_plot = f"{srr_id}_mapq.png"
bin_size = 1_000_000  # 1 Mb bins

# === FUNCTION TO PARSE AND BIN MAPQ BY POSITION ===
def get_binned_mapq(sam_path, bin_size):
    binned_mapq = defaultdict(list)
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for read in samfile.fetch(until_eof=True):
            if not read.is_unmapped:
                chrom = read.reference_name
                pos = read.reference_start
                bin_id = (chrom, pos // bin_size)
                binned_mapq[bin_id].append(read.mapping_quality)
    return binned_mapq

def average_binned_mapq(binned_mapq):
    return {bin_id: np.mean(values) for bin_id, values in binned_mapq.items()}

# === PROCESS FILES ===
cont_binned = get_binned_mapq(contaminated_file, bin_size)
clean_binned = get_binned_mapq(clean_file, bin_size)

cont_avg = average_binned_mapq(cont_binned)
clean_avg = average_binned_mapq(clean_binned)

# === MERGE AND SORT BIN POSITIONS ===
all_bins = sorted(set(cont_avg.keys()).union(clean_avg.keys()), key=lambda x: (x[0], x[1]))
x_labels = [f"{chrom}:{bin_id*bin_size//1_000_000}Mb" for chrom, bin_id in all_bins]
x_indices = np.arange(len(all_bins))

# === EXTRACT PLOT DATA ===
cont_values = [cont_avg.get(bin_id, np.nan) for bin_id in all_bins]
clean_values = [clean_avg.get(bin_id, np.nan) for bin_id in all_bins]

# === PLOT ===
plt.figure(figsize=(15, 6))
plt.plot(x_indices, cont_values, label="Contaminated", color="orange")
plt.plot(x_indices, clean_values, label="Clean", color="blue")
plt.xticks(x_indices[::max(len(x_indices)//10, 1)], x_labels[::max(len(x_labels)//10, 1)], rotation=45)
plt.ylabel("Mean MAPQ per 1Mb bin")
plt.xlabel("Genomic Position (binned)")
plt.title(f"Mapping Quality over Genomic Position - {srr_id}")
plt.legend()
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.close()

print(f"Plot saved as {output_plot}")