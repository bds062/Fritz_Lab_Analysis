#!/bin/sh

# Script to run Kraken1 on all Run #1  files beginning with "trimmed_"
# Requires: kraken1 and Kraken database path
# BDS 11/8/2024

#SBATCH --job-name=Kraken1_News
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=540:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./kraken_output.txt
#SBATCH -e ./kraken_output.txt

# Paths to Kraken and database
KRAKEN_PATH="../../programs/kraken1/kraken_files/kraken"
KRAKEN_TRANSLATE_PATH="../../programs/kraken1/kraken_files/kraken-translate"
KRAKEN_REPORT_PATH="../../programs/kraken1/kraken_files/kraken-report"
KRAKEN_DB="../../programs/kraken1/database/minikraken_20171013_4GB/"

# Directory for contaminated files list
CURRENT_DIR=$(basename "$PWD")
CONTAMINATED_FILE="${CURRENT_DIR}_contaminated_files.txt"

# Initialize contaminated files list
echo "List of contaminated files (>50% classified sequences):" > "$CONTAMINATED_FILE"

# Find all files beginning with "trimmed_" and ending in ".fastq"
for INPUT_FILE in *_1.fastq; do
    # Extract the base filename (remove directory and extension)
    BASE_NAME=$(basename "$INPUT_FILE" .fastq)
    
    # Define output filenames
    OUTPUT_KRAKEN="${BASE_NAME}.kraken"
    OUTPUT_LABELS="${BASE_NAME}.labels"
    OUTPUT_REPORT="${BASE_NAME}.report"
    
    echo "Running Kraken on $INPUT_FILE"
    
    # Run Kraken
    KRAKEN_OUTPUT=$($KRAKEN_PATH --db $KRAKEN_DB "$INPUT_FILE" > "$OUTPUT_KRAKEN" 2>&1)
    echo "$KRAKEN_OUTPUT"  # Display Kraken output for debugging

    # Extract classification percentage from Kraken output
    CLASSIFIED_PERCENT=$(echo "$KRAKEN_OUTPUT" | grep -oP '(?<=classified \().*?(?=%)' | head -n 1)
    
    # Check if classification percentage is above 50%
    if (( $(echo "$CLASSIFIED_PERCENT > 50" | bc -l) )); then
        echo "$INPUT_FILE" >> "$CONTAMINATED_FILE"
    fi

    # Run Kraken-translate on the Kraken output file
    $KRAKEN_TRANSLATE_PATH --db $KRAKEN_DB "$OUTPUT_KRAKEN" > "$OUTPUT_LABELS"
    
    # Run Kraken-report on the Kraken output file
    $KRAKEN_REPORT_PATH --db $KRAKEN_DB "$OUTPUT_KRAKEN" > "$OUTPUT_REPORT"
done

# Check if the contaminated files list is empty (other than the header)
if [ $(wc -l < "$CONTAMINATED_FILE") -eq 1 ]; then
    echo "No files in $CURRENT_DIR found to be contaminated above 50%." >> "$CONTAMINATED_FILE"
fi
