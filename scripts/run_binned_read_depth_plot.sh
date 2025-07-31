#!/bin/bash

#SBATCH --job-name=run_binned_read_depth_plot_array
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=00:02:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH --array=0-13%10   # Replace <N-1> with number of FASTA files minus one
#SBATCH -o ./logs/binned_read_depth_%A_%a.out
#SBATCH -e ./logs/binned_read_depth_%A_%a.err

# Make sure logs directory exists
mkdir -p logs

# Get list of FASTA files
FILES=(*depth.tsv)

# Get the FASTA file for this job array task
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Run the Python script on that file
python /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/Fritz_Lab_Analysis/scripts/binned_read_depth_plot.py "$FILE"