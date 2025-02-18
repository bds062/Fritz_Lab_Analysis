#!/bin/bash
# Script to run trimmomatic on all fastq files in a directory with output
# redirected to a specified file
# Requires: Trimmomatic installed
# BDS 10/25/2024

# Path to Trimmomatic JAR file
TRIMMOMATIC_JAR="../../programs/Trimmomatic-0.39/trimmomatic-0.39.jar"

# Path to the adapter file
ADAPTERS="../../programs/Trimmomatic-0.39/adapters/adapt_seq.fa"

# Trimming parameters
PHRED="PE -phred33"
TRIM_PARAMS="ILLUMINACLIP:${ADAPTERS}:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:50"

# Loop through each .fastq file in the directory
for INPUTFILE in *.fastq; do
    # Define the output filename
    OUTPUTFILE="trimmed_${INPUTFILE}"
    
    # Run Trimmomatic
    java -jar $TRIMMOMATIC_JAR $PHRED "$INPUTFILE" "$OUTPUTFILE" $TRIM_PARAMS
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Processed $INPUTFILE successfully."
    else
        echo "Failed to process $INPUTFILE."
    fi
done
