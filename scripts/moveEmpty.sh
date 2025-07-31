#!/bin/bash

# Create the directory if it doesn't exist
mkdir -p empty

# Loop through all .fastq files
for file in *.fastq; do
    # Check if the file exists and is empty
    if [ -f "$file" ] && [ ! -s "$file" ]; then
        echo "Moving empty file: $file"
        mv "$file" empty/
    fi
done

echo "Done. All empty .fastq files moved to ./empty"