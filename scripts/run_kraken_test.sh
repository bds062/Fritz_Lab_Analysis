#!/bin/sh

#Script to run kraken1 on three sets of random files that are given
#requires: kraken1 and kraken database path
#BDS 11/8/2024

#SBATCH --job-name=Kraken1_Ja
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=540:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./Ja_kraken_output.txt
#SBATCH -e ./Ja_kraken_output.txt

# List of base filenames
FILES=("SRR18709155" "SRR18709332" "SRR18709365")

# Paths to Kraken and database
KRAKEN_PATH="../../programs/kraken1/kraken_files/kraken"
KRAKEN_TRANSLATE_PATH="../../programs/kraken1/kraken_files/kraken-translate"
KRAKEN_REPORT_PATH="../../programs/kraken1/kraken_files/kraken-report"

KRAKEN_DB="../../programs/kraken1/database/minikraken_20171013_4GB/"

# Loop through each base filename and run Kraken and Kraken-translate
for FILE in "${FILES[@]}"; do
    # Process both _1 and _2 fastq files for each base filename
    for SUFFIX in "_1" "_2"; do
        # Define the input and output filenames
        INPUT_FILE="trimmed_${FILE}${SUFFIX}.fastq"
        OUTPUT_KRAKEN="${FILE}${SUFFIX}.kraken"
        OUTPUT_LABELS="${FILE}${SUFFIX}.labels"
        OUTPUT_REPORT="${FILE}${SUFFIX}.report"
        
	echo "Running kraken on ${FILE}${SUFFIX}"
        # Run Kraken on the fastq file
        $KRAKEN_PATH --db $KRAKEN_DB "$INPUT_FILE" > "$OUTPUT_KRAKEN"
        
        # Run Kraken-translate on the Kraken output file
        $KRAKEN_TRANSLATE_PATH --db $KRAKEN_DB "$OUTPUT_KRAKEN" > "$OUTPUT_LABELS"

        # Run Kraken-repot on the Kraken output file
        $KRAKEN_REPORT_PATH --db $KRAKEN_DB "$OUTPUT_KRAKEN" > "$OUTPUT_REPORT"
    done
done
