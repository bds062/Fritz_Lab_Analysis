#!/bin/bash

# Ensure the output file exists
output_file="../contamination_data.txt"
if [ ! -f "$output_file" ]; then
  # Create the file and add a header row
  echo -e "Directory\tSRR_ID\tSecond_Line_First_Value\tLast_Value" > "$output_file"
fi

# Loop through all files ending in _1.report
for file in *_1.report; do
  # Extract directory name
  directory=$(basename "$PWD")
  
  # Extract SRR ID from the file name
  srr_id=$(echo "$file" | sed -E 's/.*(SRR[0-9]+).*/\1/')

  # Extract the required values from the file
  if [ -f "$file" ]; then
    # Get the first number from the second line
    second_line_first_value=$(sed -n '2p' "$file" | awk '{print $1}')

    # Find the first instance of "S" in the 4th column and extract the species name
    species=$(awk '$4 == "S" {print substr($0, index($0, $6)); exit}' "$file")

    # Append the extracted values to the output file
    echo -e "$directory\t$srr_id\t$second_line_first_value\t$species" >> "$output_file"
  fi
done

echo "Data analysis complete. Results appended to $output_file."
