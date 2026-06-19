#!/bin/sh
# Script to run trimmomatic on the directories specified and output details
# into a specified file
# Requires: Trimmomatic installed 
# BDS 10/28/2024

#SBATCH --job-name=2019RD
#SBATCH -c 2
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=0-24:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/rd_pipeline_%j.txt
#SBATCH -e ./logs/rd_pipeline_%j.txt

# List of specific directories to process
YEARS=("2019")
# DIRECTORIES=("Ja" "LA" "MA" "MO" "NA" "NC" "Thrall" "WA" "WB") #without Tillard

# Loop through each specified directory
for YEAR in "${YEARS[@]}"; do
    # Check if the directory exists
    # if [ -d "$DIR" ]; then
        # Change into the directory
        # cd "$DIR" || continue
        echo "Processing Year: $YEAR"
        # Run the Trimmomatic script, redirecting output to [DIRECTORYNAME]_trimmomatic_output.txt
        python /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/Fritz_Lab_Analysis/scripts/depth_pileup.py trimmed_Hz_LA_${YEAR}*depth.tsv > merged_depth_${YEAR}.tsv
        python /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/Fritz_Lab_Analysis/scripts/depth_pileup_filter.py merged_depth_${YEAR}.tsv
        python /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/Fritz_Lab_Analysis/scripts/depth_pileup_windows.py filtered_1.5xmean_merged_depth_${YEAR}.tsv

        # cd ..
    # else
        # echo "Directory $DIR does not exist."
    # fi
done
