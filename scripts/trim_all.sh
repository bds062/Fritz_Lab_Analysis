#!/bin/bash
# Script to run Trimmomatic on all paired-end fastq files in specified directories
# Requires: Trimmomatic installed and adapter sequences provided
# BDS 01/14/2025

#SBATCH --job-name=Trimmomatic_All
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=1-00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./trimmomatic_output.txt
#SBATCH -e ./trimmomatic_output.txt

# Path to Trimmomatic JAR file
TRIMMOMATIC_JAR="../../programs/Trimmomatic-0.39/trimmomatic-0.39.jar"

# Path to the adapter file
ADAPTERS="../../programs/Trimmomatic-0.39/adapters/adapt_seq.fa"

# Trimming parameters
TRIM_PARAMS="ILLUMINACLIP:${ADAPTERS}:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:50"

# Directories to process
DIRECTORIES=("MO" "NA" "NC" "Thrall" "Tillard" "WA" "WB")

# Loop through each directory
for DIR in "${DIRECTORIES[@]}"; do
    echo "Processing directory: $DIR"
    
    # Change to the directory
    cd "$DIR" || { echo "Failed to enter directory $DIR. Skipping..."; continue; }
    
    # Loop through all _1.fastq files in the directory
    for INPUT1 in *_1.fastq; do
        # Extract the base name (e.g., SRR123 from SRR123_1.fastq)
        BASENAME=$(basename "$INPUT1" "_1.fastq")
        
        # Define the corresponding _2.fastq file
        INPUT2="${BASENAME}_2.fastq"
        
        # Ensure the matching _2.fastq file exists
        if [[ -f "$INPUT2" ]]; then
            # Define output file paths
            OUTPUT1_PAIRED="trimmed_${BASENAME}_1_paired.fastq"
            OUTPUT1_UNPAIRED="trimmed_${BASENAME}_1_unpaired.fastq"
            OUTPUT2_PAIRED="trimmed_${BASENAME}_2_paired.fastq"
            OUTPUT2_UNPAIRED="trimmed_${BASENAME}_2_unpaired.fastq"
            
            # Run Trimmomatic for paired-end reads
            java -jar "$TRIMMOMATIC_JAR" PE -phred33 \
                "$INPUT1" "$INPUT2" \
                "$OUTPUT1_PAIRED" "$OUTPUT1_UNPAIRED" \
                "$OUTPUT2_PAIRED" "$OUTPUT2_UNPAIRED" \
                $TRIM_PARAMS
            
            # Check if the command was successful
            if [ $? -eq 0 ]; then
                echo "Processed $BASENAME successfully."
            else
                echo "Failed to process $BASENAME."
            fi
        else
            echo "Matching file for $INPUT1 not found. Skipping..."
        fi
    done
    
    # Return to the parent directory
    cd ..
done
