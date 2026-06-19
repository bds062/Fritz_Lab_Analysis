#!/bin/bash
#SBATCH --job-name=fastqc

#SBATCH --time=24:00:00
#SBATCH --output=./logs/fastqc2.log

#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb

#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu

# set -euo pipefail

# Path to FastQC executable (relative to the directory where you run the script)
FASTQC_BIN="../../programs/FastQC/fastqc"

# Loop over all .fastq.gz files in the current directory
for fq in trimmed*.fastq.gz; do
    # Skip if no files match
    [[ -e "$fq" ]] || continue

    # Expected FastQC output zip (FastQC uses input basename + "_fastqc.zip")
    base=$(basename "$fq" .fastq.gz)
    out_zip="${base%.*}_fastqc.zip"      # removes one extension (.gz)
    # If you want to remove both .fastq.gz explicitly:
    # out_zip="${base%.fastq.gz}_fastqc.zip"
    echo "Checking for $out_zip ..."
    if [[ -f "$out_zip" ]]; then
        echo "Skipping $fq (found $out_zip)"
    else
        echo "Running FastQC on $fq ..."
        "$FASTQC_BIN" -t 16 "$fq"
    fi
done
