#!/bin/bash
# Script to find these specific SRR files in a group of fastq files
# BDS 12/1/24

# List of SRR strings to search for in file names
SRR_LIST=("SRR18709194" "SRR18709232" "SRR18709361" "SRR18709343")

# Start from the current directory
START_DIR=$(pwd)

# Iterate through all files in subdirectories
find . -type f | while read -r file; do
    # Get the base name of the file
    base_name=$(basename "$file")
    
    # Check if the file name contains any of the SRR strings
    for srr in "${SRR_LIST[@]}"; do
        if [[ "$base_name" == *"$srr"* ]]; then
            # Print the relative path of the file
            echo "Match found in: ${file#$START_DIR/}"
            break
        fi
    done
done
