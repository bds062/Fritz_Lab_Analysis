#!/bin/bash
#SBATCH --job-name=Depth
#SBATCH --array=0-91%20
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=12:00:00
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/depth_%A_%a.out
#SBATCH -e ./logs/depth_%A_%a.out

set -euo pipefail

SAMTOOLS="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/samtools-1.21/samtools"

# -----------------------------
# Build array of BAM files
# -----------------------------
bam_array=($(ls *.bam))

bam=${bam_array[$SLURM_ARRAY_TASK_ID]}
base_name=$(basename "$bam" .bam)

depth_out="${base_name}_depth.tsv"
bam_index="${bam}.bai"

# -----------------------------
# Function to clean partial outputs
# -----------------------------
cleanup_partial() {
    echo "Cleaning partial depth for $base_name..."
    rm -f "$depth_out"
}

# -----------------------------
# Check BAM integrity
# -----------------------------
bam_ok=false
if [[ -s "$bam" ]] && $SAMTOOLS quickcheck "$bam" >/dev/null 2>&1; then
    bam_ok=true
fi

if ! $bam_ok; then
    echo "Skipping $base_name: BAM is missing or corrupted."
    exit 0
fi

# -----------------------------
# Check if depth already exists
# -----------------------------
if [[ -s "$depth_out" ]]; then
    echo "Skipping $base_name: depth file already exists."
    exit 0
fi

# If partial depth exists, remove it
if [[ -f "$depth_out" ]]; then
    cleanup_partial
fi

# -----------------------------
# Ensure BAM is indexed
# -----------------------------
if [[ ! -f "$bam_index" ]]; then
    echo "Indexing BAM for $base_name..."
    $SAMTOOLS index "$bam"
fi

# -----------------------------
# Run samtools depth
# -----------------------------
echo "Computing depth for $base_name..."
$SAMTOOLS depth -@ 16 -a "$bam" > "$depth_out"

echo "Done $base_name!"