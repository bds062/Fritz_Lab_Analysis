#!/bin/bash

# Script to run Bowtie2 on trimmed paired-end files using SLURM job arrays
# Each array task aligns one sample
# Requires: Bowtie2, Trimmomatic trimmed files

#TO USE: 
#ls trimmed_*_1_paired.fastq | wc -l
#sbatch --array=0-23 allign_trimmed.sh

# BDS 3/22/25 (array version)

#SBATCH --job-name=Bowtie2ContaminatedTrims
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=6:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o logs/bowtie_%A_%a.out
#SBATCH -e logs/bowtie_%A_%a.err

# ---------- CONFIG ----------
BOWTIE="../../programs/bowtie2/bowtie2"
INDEX="../other/zea"
THREADS=16
# ----------------------------

# Make sure log directory exists
mkdir -p logs

# Build sample list ONCE (sorted for reproducibility)
SAMPLES=($(ls trimmed_*_1_paired.fastq | \
           sed -E 's/trimmed_(SRR[0-9]+)_1_paired.fastq/\1/' | \
           sort))

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

if [[ -z "$SAMPLE" ]]; then
    echo "No sample found for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

R1="trimmed_${SAMPLE}_1_paired.fastq"
R2="trimmed_${SAMPLE}_2_paired.fastq"
OUT="trimmed_${SAMPLE}.sam"

echo "[$(date)] Running Bowtie2 for ${SAMPLE}"
echo "R1=${R1}"
echo "R2=${R2}"

$BOWTIE \
    -x ${INDEX} \
    -1 ${R1} \
    -2 ${R2} \
    --very-sensitive \
    -p ${THREADS} \
    -S ${OUT}

echo "[$(date)] Finished ${SAMPLE}"