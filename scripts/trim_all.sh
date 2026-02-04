#!/bin/bash
# Script to run Trimmomatic on all paired-end fastq files in specified directories
# Skips Trimmomatic if trimmed files are already present
# BDS 01/14/2025

#SBATCH --job-name=Trimmomatic_All
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=0-08:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/trimmomatic_output.txt
#SBATCH -e ./logs/trimmomatic_output.txt

# Path to Trimmomatic JAR file
TRIMMOMATIC_JAR="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/Trimmomatic-0.39/trimmomatic-0.39.jar"

# Path to the adapter file
ADAPTERS="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/Trimmomatic-0.39/adapters/adapt_seq.fa"

# Trimming parameters
TRIM_PARAMS="ILLUMINACLIP:${ADAPTERS}:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:50"

# Directories to process
DIRECTORIES=("/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/MDTaylor")

# Loop through each directory
for DIR in "${DIRECTORIES[@]}"; do
    echo "Processing directory: $DIR"
    
    cd "$DIR" || { echo "Failed to enter directory $DIR. Skipping..."; continue; }
    
    for INPUT1 in *_1.fastq; do
        BASENAME=$(basename "$INPUT1" "_1.fastq")
        INPUT2="${BASENAME}_2.fastq"
        
        if [[ -f "$INPUT2" ]]; then
            OUTPUT1_PAIRED="trimmed_${BASENAME}_R1_paired.fastq"
            OUTPUT1_UNPAIRED="trimmed_${BASENAME}_R1_unpaired.fastq"
            OUTPUT2_PAIRED="trimmed_${BASENAME}_R2_paired.fastq"
            OUTPUT2_UNPAIRED="trimmed_${BASENAME}_R2_unpaired.fastq"

            # Check if all output files already exist
            if [[ -f "$OUTPUT1_PAIRED" && -f "$OUTPUT2_PAIRED" ]]; then
                echo "Trimmed files for $BASENAME already exist. Skipping..."
                continue
            fi


            java -jar "$TRIMMOMATIC_JAR" PE -phred33 \
                "$INPUT1" "$INPUT2" \
                "$OUTPUT1_PAIRED" "$OUTPUT1_UNPAIRED" \
                "$OUTPUT2_PAIRED" "$OUTPUT2_UNPAIRED" \
                $TRIM_PARAMS
            
            if [ $? -eq 0 ]; then
                echo "Processed $BASENAME successfully."
            else
                echo "Failed to process $BASENAME."
            fi
        else
            echo "Matching file for $INPUT1 not found. Skipping..."
        fi
    done
    
    cd ..
done
