import os
import statistics
import pandas as pd

# Directory containing the SAM files
directory = "./"  # Change if needed

# Output file
output_file = "mapping_quality_summary.tsv"

# Helper to extract MAPQ values and count reads
def extract_mapq(file_path):
    mapq_values = []
    mapped_count = 0
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                fields = line.strip().split('\t')
                try:
                    mapq = int(fields[4])
                    if mapq > 0:  # Only include mapped reads
                        mapq_values.append(mapq)
                        mapped_count += 1
                except (IndexError, ValueError):
                    continue
    return mapq_values, mapped_count

# Collect SAM files and IDs
sam_files = [f for f in os.listdir(directory) if f.endswith('.sam')]
sample_ids = set(f.replace('trimmed_', '').replace('.sam', '') for f in sam_files)

# Extract MAPQ and compute stats
records = []
for sample_id in sample_ids:
    clean_file = os.path.join(directory, f"{sample_id}.sam")
    contaminated_file = os.path.join(directory, f"trimmed_{sample_id}.sam")

    clean_mean = clean_median = cont_mean = cont_median = ''
    clean_reads = cont_reads = ''

    if os.path.exists(clean_file):
        vals, clean_reads = extract_mapq(clean_file)
        if vals:
            clean_mean = round(statistics.mean(vals), 2)
            clean_median = round(statistics.median(vals), 2)

    if os.path.exists(contaminated_file):
        vals, cont_reads = extract_mapq(contaminated_file)
        if vals:
            cont_mean = round(statistics.mean(vals), 2)
            cont_median = round(statistics.median(vals), 2)

    records.append({
        'ID': sample_id,
        'Contaminated Mean MAPQ': cont_mean,
        'Contaminated Median MAPQ': cont_median,
        'Contaminated Reads': cont_reads,
        'Clean Mean MAPQ': clean_mean,
        'Clean Median MAPQ': clean_median,
        'Clean Reads': clean_reads
    })

# Save output
df = pd.DataFrame(records)
df.to_csv(output_file, sep='\t', index=False)
