#!/bin/bash

#SBATCH --job-name=extractContamReads
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./extractContamReads.txt
#SBATCH -e ./extractContamReads.txt

# Path to BBMap's filterbyname.sh
FILTER="../../programs/BBMap/bbmap/filterbyname.sh"

# Loop through each ID in IDs.txt
while read -r ID; do
  echo "Processing $ID..."

  $FILTER \
    in1="trimmed_${ID}_1_paired.fastq" \
    in2="trimmed_${ID}_2_paired.fastq" \
    out1="contaminated_${ID}_1.fastq" \
    out2="contaminated_${ID}_2.fastq" \
    names="clean_${ID}_1.fastq,clean_${ID}_2.fastq" \
    include=f

done < IDs.txt

echo "All SRR IDs processed."