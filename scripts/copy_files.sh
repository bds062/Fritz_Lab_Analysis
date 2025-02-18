#!/bin/bash

# List of SRR IDs to search for
srr_ids=(
  "SRR18709133" "SRR18709131" "SRR18709241" "SRR18709246"
  "SRR18709340" "SRR18709342" "SRR18709344" "SRR18709348"
  "SRR18709355" "SRR18709356" "SRR18709218" "SRR18709220"
  "SRR18709221" "SRR18709269" "SRR18709189" "SRR18709192"
  "SRR18709205" "SRR18709206" "SRR18709212" "SRR18709318"
  "SRR18709326"
)

# Create the contaminated directory if it doesn't exist

# Loop through each SRR ID and find matching files
for srr in "${srr_ids[@]}"; do
  # Find files containing the SRR ID in their name within subdirectories
  find . -type f -name "*${srr}*" -exec cp {} ./contaminated/ \;
done

echo "Copying complete. Contaminated files are in ./contaminated"

