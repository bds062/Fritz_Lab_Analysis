#!/bin/bash

#SBATCH --job-name=run_read_depth_alg
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=01:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH --array=0-13%10   # Replace <N-1> with number of FASTA files minus one
#SBATCH -o ./logs/read_depth_alg_%A_%a.out
#SBATCH -e ./logs/read_depth_alg_%A_%a.err

# Make sure logs directory exists
mkdir -p logs

# Get list of FASTA files
FILES=(*depth.tsv)

# Get the FASTA file for this job array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Run the Python script on that file
python /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/Fritz_Lab_Analysis/scripts/read_depth_alg.py "$FILE"