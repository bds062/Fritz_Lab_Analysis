#!/bin/bash
set -euo pipefail

DEST="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/sandboxes/rdviewer_full/full_chrom/"
# REGION="CM040209.1:8244001-8266000"
REGION="CM040209.1:0-13341776"
OUTPUT_TSV="$DEST/meandepth_summary.tsv"

# Verify destination exists
if [[ ! -d "$DEST" ]]; then
    echo "ERROR: Destination directory does not exist: $DEST"
    exit 1
fi

# Check for BAM files
shopt -s nullglob
bam_files=(*.bam)
if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "ERROR: No .bam files found in current directory: $(pwd)"
    exit 1
fi

echo "Found ${#bam_files[@]} BAM file(s):"
for f in "${bam_files[@]}"; do echo "  $f"; done
echo ""

# Step 1: Subset each BAM to the target region and index it
for bam in "${bam_files[@]}"; do
    echo "[1/2] Subsetting $bam -> $DEST/$bam"
    samtools view -b "$bam" "$REGION" > "$DEST/$bam"
    echo "      Indexing..."
    samtools index "$DEST/$bam"
done

echo ""
echo "All BAMs subsetted and indexed."
echo ""

# Step 2: Run samtools coverage on each trimmed BAM and collect meandepth
echo -e "ID\tmeandepth" > "$OUTPUT_TSV"

for bam in "$DEST"/*.bam; do
    id=$(basename "$bam" .bam)
    meandepth=$(samtools coverage -r "$REGION" "$bam" | awk 'NR==2 {print $7}')
    echo -e "$id\t$meandepth" >> "$OUTPUT_TSV"
    echo "[2/2] $id -> meandepth: $meandepth"
done

echo ""
echo "Done! Summary TSV written to:"
echo "  $OUTPUT_TSV"
echo ""
cat "$OUTPUT_TSV"