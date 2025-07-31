#!/bin/bash
#SBATCH --job-name=Reads_to_bam_Tillard
#SBATCH --array=0-17%10
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=00-10:00:00
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/align_%A_%a.out
#SBATCH -e ./logs/align_%A_%a.err

# === Setup ===
BOWTIE="../../programs/bowtie2/bowtie2"
TAGS="-t -x ../other/zea"
SAMTOOLS="../../programs/samtools-1.21/samtools"

# Ensure logs directory exists
mkdir -p logs

# === Get list of paired read files ===
file_array=($(ls trimmed_*_1_paired.fastq))

# === Get current sample based on SLURM_ARRAY_TASK_ID ===
file1=${file_array[$SLURM_ARRAY_TASK_ID]}
base_name=$(basename "$file1" | sed -E 's/trimmed_(SRR[0-9]+)_1_paired.fastq/\1/')
file2="trimmed_${base_name}_2_paired.fastq"

# === Alignment and BAM Conversion ===
if [[ -f "$file1" && -f "$file2" ]]; then
    echo "Running Bowtie2 for $base_name..."
    $BOWTIE $TAGS -1 "$file1" -2 "$file2" --very-sensitive -S "trimmed_${base_name}.sam"

    echo "Sorting and converting SAM to BAM..."
    $SAMTOOLS view -@ 4 -S -h -u "trimmed_${base_name}.sam" | \
    $SAMTOOLS sort -@ 4 -T "trimmed_${base_name}.temp" -o "trimmed_${base_name}.bam"

    $SAMTOOLS index "trimmed_${base_name}.bam"
    $SAMTOOLS stats "trimmed_${base_name}.bam" > "trimmed_${base_name}.stats"
else
    echo "Skipping $base_name: Missing paired files."
fi