#!/bin/bash

# Script to seperate interleaved files in a directory
# Requires: BBMAP or BBTOOLS installation
# BDS 3/13/2024

#SBATCH --job-name=InterleaveSplit
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=540:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./interleave_split.txt
#SBATCH -e ./interleave_split.txt

# Define the path to BBMap reformat tool
BBMAP_PATH="../../programs/BBMap/bbmap/reformat.sh"

# Loop through all matching FASTQ files in the current directory
for file in unmapped_MM_*.fastq; do
    # Extract the SRR_ID from the filename (assuming format: unmapped_MM_[SRR_ID].fastq)
    SRR_ID=$(echo "$file" | sed -E 's/unmapped_MM_(.+)\.fastq/\1/')

    # Define output filenames
    OUT1="clean_${SRR_ID}_1.fastq"
    OUT2="clean_${SRR_ID}_2.fastq"

    # Run BBMap reformat.sh
    echo "Processing $file -> $OUT1 and $OUT2"
    $BBMAP_PATH in="$file" out1="$OUT1" out2="$OUT2"

    # Optional: Check for success and report
    if [[ $? -eq 0 ]]; then
        echo "Successfully processed $file"
    else
        echo "Error processing $file"
    fi
done

echo "All files processed."
