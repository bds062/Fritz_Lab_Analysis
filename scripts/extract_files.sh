#!/bin/bash

# Create the destination directory if it does not exist
mkdir -p contaminated

# List of SRR values
srr_ids=(
  "SRR18709133" "SRR18709131" "SRR18709241" "SRR18709246" "SRR18709340"
  "SRR18709342" "SRR18709344" "SRR18709348" "SRR18709355" "SRR18709356"
  "SRR18709218" "SRR18709220" "SRR18709221" "SRR18709269" "SRR18709189"
  "SRR18709192" "SRR18709205" "SRR18709206" "SRR18709212" "SRR18709318"
  "SRR18709326"
)

# Iterate over each directory in the current directory
for dir in */; do
  # Change into the directory
  cd "$dir" || continue
  
  echo "Searching in directory: $dir"
  
  # Iterate through each SRR_ID and copy matching files
  for srr_id in "${srr_ids[@]}"; do
    # Check for matching files
    for file in *"${srr_id}"*.fastq; do
      if [ -f "$file" ]; then
        echo "Found file: $file for SRR_ID: $srr_id"
        cp "$file" ../contaminated/
      fi
    done
  done
  
  # Return to the parent directory
  cd ..
done

echo "Copying complete. Contaminated files are in the ./contaminated directory."
