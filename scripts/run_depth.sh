#!/bin/sh
# Script to turn all bam files into depth
# Requires: samtools-1.21
# BDS 4/23/2025

#SBATCH --job-name=samtobam
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/depth.out
#SBATCH -e ./logs/depth.out

# Set the path to samtools
SAMTOOLS="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/samtools-1.21/samtools"

for bam in *.bam; do
    # Skip if no BAM files are found
    [[ -e "$bam" ]] || { echo "No BAM files found."; exit 0; }

    # Strip the .bam extension to get the ID
    id="${bam%.bam}"
    echo "running $id"
    # Run samtools depth and write to ID_depth.tsv
    "$SAMTOOLS" depth "$bam" > "${id}_depth.tsv"
done
