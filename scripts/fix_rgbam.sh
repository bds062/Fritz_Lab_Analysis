#!/bin/bash
#!/bin/sh
# Script to turn all .sam in a directory to .bam
# Requires: samtools-1.21
# BDS 4/23/2025

#SBATCH --job-name=fix_rg
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=4:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/fix_rg.txt
#SBATCH -e ./logs/fix_rg.txt
source activate
conda init
conda activate /nfshomes/bds062/miniconda3/envs/delly2/

set -euo pipefail



SAMTOOLS=samtools

for f in RG_*.bam; do
    echo "Fixing $f"

    tmp="${f}.tmp.bam"

    # Convert SAM → BAM (BGZF compressed)
    $SAMTOOLS view -@ 16 -b -h "$f" -o "$tmp"

    # Replace original file
    mv "$tmp" "$f"

    # Index the BAM
    $SAMTOOLS index "$f"
done

echo "All files converted to true BAM format."