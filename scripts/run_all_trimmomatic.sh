#!/bin/sh
# Script to run trimmomatic on the directories specified and output details
# into a specified file
# Requires: Trimmomatic installed 
# BDS 10/28/2024

#SBATCH --job-name=trim_all2
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=540:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./trimmomatic_output.txt
#SBATCH -e ./trimmomatic_output.txt

# List of specific directories to process
DIRECTORIES=("Ja" "LA" "MA" "MO" "NA" "NC" "Thrall" "Tillard" "WA" "WB")

# Loop through each specified directory
for DIR in "${DIRECTORIES[@]}"; do
    # Check if the directory exists
    if [ -d "$DIR" ]; then
        # Change into the directory
        cd "$DIR" || continue
        
        # Run the Trimmomatic script, redirecting output to [DIRECTORYNAME]_trimmomatic_output.txt
        ../../scripts/run_trimmomatic.sh > "${DIR}_trimmomatic_output.txt" 2>&1
        
        # Go back to the main directory
        cd ..
    else
        echo "Directory $DIR does not exist."
    fi
done
