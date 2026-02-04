#!/bin/bash
#SBATCH --job-name=Reads_to_bam_HzeaFiles
#SBATCH --array=0-19%10
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=2:00:00
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/align_%A_%a.out
#SBATCH -e ./logs/align_%A_%a.out

# === Setup ===
BOWTIE="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/bowtie2/bowtie2"
TAGS="-p 16 -t -x /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/zea"
SAMTOOLS="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/samtools-1.21/samtools"

# Ensure logs directory exists
mkdir -p logs

#normal .fastq
# # === Get list of paired read files ===
# file_array=($(ls trimmed_*_R1_paired.fastq))

# # === Get current sample based on SLURM_ARRAY_TASK_ID ===
# file1=${file_array[$SLURM_ARRAY_TASK_ID]}
# base_name=$(basename "$file1" | sed -E 's/trimmed_(.*)_R1_paired.fastq/\1/')
# file2="trimmed_${base_name}_R2_paired.fastq"

# === Get list of paired read files ===
file_array=($(ls trimmed_*_R1_paired.fastq))

# === Get current sample based on SLURM_ARRAY_TASK_ID ===
file1=${file_array[$SLURM_ARRAY_TASK_ID]}
base_name=$(basename "$file1" | sed -E 's/trimmed_(.*)_R1_paired.fastq/\1/')
file2="trimmed_${base_name}_R2_paired.fastq"

# === Alignment and BAM Conversion ===
if [[ -f "$file1" && -f "$file2" ]]; then
    echo "Running Bowtie2 for $base_name..."
    $BOWTIE $TAGS -1 "$file1" -2 "$file2" --very-sensitive -S "trimmed_${base_name}.sam"

    echo "Sorting and converting SAM to BAM..."
    $SAMTOOLS view -@ 16 -S -h -u "trimmed_${base_name}.sam" | \
    $SAMTOOLS sort -@ 16 -T "trimmed_${base_name}.temp" -o "trimmed_${base_name}.bam"
    $SAMTOOLS index "trimmed_${base_name}.bam"
    $SAMTOOLS depth -@ 16 -o "trimmed_${base_name}_depth.tsv" "trimmed_${base_name}.bam"
    # $SAMTOOLS stats "trimmed_${base_name}.bam" > "trimmed_${base_name}.stats"
    rm trimmed_${base_name}.sam
else
    echo "Skipping $base_name: Missing paired files."
fi

echo "Done!"