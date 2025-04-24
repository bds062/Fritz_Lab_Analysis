#!/bin/sh
# Script to turn all .sam in a directory to .bam
# Requires: samtools-1.21
# BDS 4/23/2025

#SBATCH --job-name=samtobam
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=1-00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./sam_to_bam_output.txt
#SBATCH -e ./sam_to_bam_output.txt

# Set the path to samtools
SAMTOOLS="../../programs/samtools-1.21/samtools"

# Loop through each .sam file in the current directory
for sam_file in *.sam; do
  # Remove the .sam extension to get the base name
  base_name="${sam_file%.sam}"

  # Run samtools to convert .sam to .bam
  echo "Converting $sam_file to $base_name.bam..."
  $SAMTOOLS view -b -S "$sam_file" > "${base_name}.bam"
done

echo "All SAM files have been converted to BAM."