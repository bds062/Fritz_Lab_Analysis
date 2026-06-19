#!/bin/bash
#SBATCH --job-name=Reads_to_bam_HzeaFiles
#SBATCH --array=0-10%11
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=5:00:00
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/align_pipeline_%A_%a.out
#SBATCH -e ./logs/align_pipeline_%A_%a.out

set -euo pipefail

BOWTIE="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/bowtie2/bowtie2"
TAGS="-p 16 -t -x /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/zea"
SAMTOOLS="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/samtools-1.21/samtools"

# Build array of R1 files
file_array=($(ls trimmed_*_R1_paired.fastq | sort))

# Current sample
file1=${file_array[$SLURM_ARRAY_TASK_ID]}
base_name=$(basename "$file1" | sed -E 's/trimmed_(.*)_R1_paired.fastq/\1/')
file2="trimmed_${base_name}_R2_paired.fastq"
bam_out="trimmed_${base_name}.bam"
bam_index="${bam_out}.bai"
depth_out="trimmed_${base_name}_depth.tsv"
# sam_out="trimmed_${base_name}.sam"

# Only run if both reads exist
if [[ ! -f "$file1" || ! -f "$file2" ]]; then
    echo "Skipping $base_name: missing paired files."
    exit 0
fi

# Function to clean partial outputs
cleanup_partial() {
    echo "Cleaning partial outputs for $base_name..."
    rm -f "$bam_out" "$bam_index" "$depth_out" "trimmed_${base_name}.temp"* "${bam_out}.tmp" #"$sam_out" 1
}

# Check if BAM and depth are both good
bam_ok=false
depth_ok=false

if [[ -s "$bam_out" ]] && $SAMTOOLS quickcheck "$bam_out" >/dev/null 2>&1; then
    bam_ok=true
else
    bam_ok=false
fi

if [[ -s "$depth_out" ]]; then
    depth_ok=true
else
    depth_ok=false
fi

# If all required outputs are good, skip
if $bam_ok && $depth_ok; then
    echo "Skipping $base_name: valid BAM and depth file already exist."
    exit 0
fi

# If some, but not all, outputs exist, clean them up before rerun
if ! $bam_ok || ! $depth_ok; then
    cleanup_partial
fi

# echo "Running Bowtie2 for $base_name..."
# $BOWTIE $TAGS -1 "$file1" -2 "$file2" --very-sensitive -S "$sam_out"
# if [[ ! -s "$sam_out" ]]; then
#     echo "ERROR: SAM not created for $base_name"
#     exit 1
# fi

# echo "Sorting and converting SAM to BAM..."

# bam_out_temp="trimmed_${base_name}_temp.bam"

# $SAMTOOLS view -@ 16 -h -b "$sam_out" | \
# $SAMTOOLS sort -@ 16 -T "trimmed_${base_name}.temp" -o "$bam_out_temp"

# mv "$bam_out_temp" "$bam_out"

echo "Running Bowtie2 → BAM pipeline for $base_name..."

bam_out_temp="trimmed_${base_name}_temp.bam"

$BOWTIE $TAGS -1 "$file1" -2 "$file2" --very-sensitive 2> "trimmed_${base_name}.bowtie2.log" | \
$SAMTOOLS view -@ 16 -b -h | \
$SAMTOOLS sort -@ 16 -T "trimmed_${base_name}.temp" -o "$bam_out_temp"

mv "$bam_out_temp" "$bam_out"

if ! $SAMTOOLS quickcheck "$bam_out"; then
    echo "ERROR: BAM is invalid for $base_name"
    exit 1
fi

echo "Adding read group to BAM..."
$SAMTOOLS addreplacerg \
    -@ 16 \
    -O bam \
    -r "@RG\tID:${base_name}\tSM:${base_name}\tPL:ILLUMINA" \
    -o "${bam_out}.tmpRG" \
    "$bam_out"

if ! $SAMTOOLS quickcheck "${bam_out}.tmpRG"; then
    echo "ERROR: BAM corrupted after read group step for $base_name"
    exit 1
fi

mv "${bam_out}.tmpRG" "$bam_out"


echo "Indexing BAM..."
$SAMTOOLS index "$bam_out"

echo "Computing depth..."
depth_out_temp="trimmed_${base_name}_temp_depth.tsv"
$SAMTOOLS depth -@ 16 -o "$depth_out_temp" "$bam_out"
mv "$depth_out_temp" "$depth_out"

echo "Done $base_name!"
