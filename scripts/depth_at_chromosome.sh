#!/bin/bash
set -euo pipefail

DEST="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/sandboxes/rdviewer_full/full_chrom/"
CHR="CM040209.1"
OUTPUT_TSV="$DEST/meandepth_chromosome.tsv"

# Verify destination exists
if [[ ! -d "$DEST" ]]; then
    echo "ERROR: Destination directory does not exist: $DEST"
    exit 1
fi

# Find BAM files in current directory
shopt -s nullglob
bam_files=(*.bam)
if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "ERROR: No .bam files found in current directory: $(pwd)"
    exit 1
fi

echo "Found ${#bam_files[@]} BAM file(s):"
for f in "${bam_files[@]}"; do echo "  $f"; done
echo ""

# Write header
echo -e "ID	meandepth" > "$OUTPUT_TSV"

# For each BAM, run samtools coverage and extract the mean depth for the chromosome
for bam in "${bam_files[@]}"; do
    id=$(basename "$bam" .bam)
    echo "[*] Computing coverage for $id on chromosome $CHR"
    # samtools coverage prints a header then one line per reference; match the reference name
    meandepth=$(samtools coverage "$bam" | awk -v chr="$CHR" '$1==chr {print $7}')

    # If samtools didn't report the chromosome, record NA
    if [[ -z "${meandepth:-}" ]]; then
        meandepth="NA"
    fi

    echo -e "$id	$meandepth" >> "$OUTPUT_TSV"
    echo "    -> meandepth: $meandepth"
done


echo ""
echo "Done! Summary TSV written to:"
echo "  $OUTPUT_TSV"
echo ""
cat "$OUTPUT_TSV"
