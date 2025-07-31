import os
import pysam
import csv

output_file = "mapped_read_summary.tsv"

# Prepare output TSV
with open(output_file, "w", newline="") as tsvfile:
    writer = csv.writer(tsvfile, delimiter="\t")
    writer.writerow(["Sample", "Total_Reads", "Mapped_Reads", "Percent_Mapped"])

    # Loop through SAM files
    for filename in os.listdir("."):
        if filename.startswith("contaminated_S") and filename.endswith(".sam"):
            sample_id = filename.replace("contaminated_", "").replace(".sam", "")
            try:
                with pysam.AlignmentFile(filename, "r") as samfile:
                    total = 0
                    mapped = 0
                    for read in samfile.fetch(until_eof=True):
                        total += 1
                        if not read.is_unmapped:
                            mapped += 1
                    percent = (mapped / total) * 100 if total > 0 else 0
                    writer.writerow([sample_id, total, mapped, f"{percent:.2f}"])
            except Exception as e:
                print(f"Error processing {filename}: {e}")