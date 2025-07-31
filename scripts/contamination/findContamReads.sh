#!/bin/bash

#SBATCH --job-name=extractContamReads
#SBATCH --array=0-20%10  # Adjust range to number of IDs minus 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/extractContamReads_%A_%a.out
#SBATCH -e ./logs/extractContamReads_%A_%a.err

# Load BBMap if needed (optional)
# module load BBMap

# Path to BBMap filter tool
FILTER="../../programs/BBMap/bbmap/filterbyname.sh"

# Get the SRR ID for this task
ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ids.txt)

echo "Processing $ID..."

$FILTER \
  in1="trimmed_${ID}_1_paired.fastq" \
  in2="trimmed_${ID}_2_paired.fastq" \
  out1="contaminated_${ID}_1.fastq" \
  out2="contaminated_${ID}_2.fastq" \
  names="clean_${ID}_1.fastq,clean_${ID}_2.fastq" \
  include=f

echo "Done with $ID"