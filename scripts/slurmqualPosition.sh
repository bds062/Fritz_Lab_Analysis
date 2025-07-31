#!/bin/bash

#SBATCH --job-name=qualPositionPlotting
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./slurmqualPosition.txt
#SBATCH -e ./slurmqualPosition.txt

source activate base
# Load Python module if necessary
module load Python3/3.13.2  # Change to your HPC's Python version with matplotlib/pysam
pip3 install matplotlib
pip3 install pysam

# Activate your virtualenv if needed

# Loop through SRR IDs
while read srr_id; do
    echo "Processing $srr_id..."
    python3 qualPosition.py "$srr_id"
done < ids.txt