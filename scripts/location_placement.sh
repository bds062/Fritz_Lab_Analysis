#!/bin/bash
#This code takes a data file of SRR values and their corresponding locations and places each file into its
#location directory.
#October 12, 2024 BDS

# Define the input SRRInfo.txt file
info_file="../SRRInfo.txt"

# Check if the SRRInfo.txt file exists
if [ ! -f "$info_file" ]; then
  echo "Error: $info_file not found."
  exit 1
fi

# Loop through each line of the SRRInfo.txt file
while IFS=$'\t' read -r _ _ srr_id location _; do
  # Check if the directory for the location already exists
  if [ ! -d "$location" ]; then
    mkdir -p "$location"
    echo "Created directory: $location"
  fi
  
  # Move the files matching the SRR ID to the corresponding directory
  mv ${srr_id}_*.fastq "$location"/
  
done < "$info_file"

echo "Files organized into directories successfully."
  fi
  
  # Move the files matching the SRR ID to the corresponding directory
  mv ${srr_id}_*.fastq "$location"/
  
done < "$info_file"
# Confirm program
echo "Files organized into directories successfully."
