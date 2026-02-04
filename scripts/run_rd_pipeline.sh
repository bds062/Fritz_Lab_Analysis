#!/bin/sh
# Script to run trimmomatic on the directories specified and output details
# into a specified file
# Requires: Trimmomatic installed 
# BDS 10/28/2024

#SBATCH --job-name=rd_pipeline
#SBATCH -c 2
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=0-15:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/rd_pipeline.txt
#SBATCH -e ./logs/rd_pipeline.txt

# List of specific directories to process
# DIRECTORIES=("Ja" "LA" "MA" "MO" "NA" "NC" "Thrall" "Tillard" "WA" "WB")
# DIRECTORIES=("Ja" "LA" "MA" "MO" "NA" "NC" "Thrall" "WA" "WB") #without Tillard

# Loop through each specified directory
# for DIR in "${DIRECTORIES[@]}"; do
    # Check if the directory exists
    # if [ -d "$DIR" ]; then
        # Change into the directory
        # cd "$DIR" || continue
        echo "Processing directory: $PWD"
        # Run the Trimmomatic script, redirecting output to [DIRECTORYNAME]_trimmomatic_output.txt
        python ../../Fritz_Lab_Analysis/scripts/depth_pileup.py trimmed_*depth.tsv > merged_depth.tsv
        python ../../Fritz_Lab_Analysis/scripts/depth_pileup_filter.py merged_depth.tsv 
        python ../../Fritz_Lab_Analysis/scripts/depth_pileup_windows.py filtered_1.5xmean_merged_depth.tsv 

        cd ..
    # else
        # echo "Directory $DIR does not exist."
    # fi
# done
