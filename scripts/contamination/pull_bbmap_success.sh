#!/bin/bash

# Input file
input_file="./bbmap_output.txt"
# Output file
output_file="./bbmap_success.txt"

# Ensure the output file is empty or create a new one
> "$output_file"

# Parse the input file
while IFS= read -r line; do
  # Look for lines with "Processing file:"
  if [[ $line == Processing\ file:* ]]; then
    # Extract the SRRID_# (trimmed filename without prefix)
    srr_id=$(echo "$line" | sed -E 's/.*trimmed_([^ ]+).fastq.*/\1/')
  fi
  
  # Look for lines with "Percent mapped:"
  if [[ $line == Percent\ mapped:* ]]; then
    # Extract the percent mapped value
    percent_mapped=$(echo "$line" | awk '{print $3}')
    
    # Append the SRRID_# and percent mapped value to the output file
    echo -e "$srr_id\t$percent_mapped" >> "$output_file"
  fi
done < "$input_file"

echo "bbmap_success.txt has been created."
