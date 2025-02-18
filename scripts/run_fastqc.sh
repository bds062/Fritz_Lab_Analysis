#!/bin/bash

# Script to run FastQC on 3 sampl .fastq files in a directory and zips
# their returned .fastq files
# Requires: FastQC installed
# BDS 10/8/2024

# Check if at least 3 parameters are passed
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <SRR_ID1> <SRR_ID2> <SRR_ID3> <SRR_ID4>"
  exit 1
fi

# Extract the SRR IDs from the arguments
srr1=$1
srr2=$2
srr3=$3
srr4=$4

# List of SRR IDs to process
srr_ids=($srr1 $srr2 $srr3 $srr4)

# Loop through the SRR IDs and run fastqc on the corresponding files
for srr in "${srr_ids[@]}"; do
  # Check if both files (_1.fastq and _2.fastq) exist for each SRR ID
  if [[ -f "$mapped_{srr}_1.fastq" && -f "$mapped_{srr}_2.fastq" ]]; then
    echo "Running fastqc on ${srr}_1.fastq and ${srr}_2.fastq"
    ../../FastQC/fastqc "$mapped_{srr}_1.fastq" "$mapped_{srr}_2.fastq"
  else
    echo "Warning: Files for ${srr} not found or incomplete."
  fi
done

echo "Done Running FASTQC Analysis"

# Check for the generated *_fastqc.html files and zip them
html_files=()

for srr in "${srr_ids[@]}"; do
  # Add fastqc result files for _1 and _2 to the list
  if [[ -f "$mapped_{srr}_1_fastqc.html" && -f "$mapped_{srr}_2_fastqc.html" ]]; then
    html_files+=("$mapped_{srr}_1_fastqc.html" "$mapped_{srr}_2_fastqc.html")
  else
    echo "Warning: fastqc HTML files for ${srr} not found."
  fi
done

# Zip the HTML files if there are 6 of them
if [ "${#html_files[@]}" -eq 8 ]; then
  zip fastqc_results.zip "${html_files[@]}"
  echo "Zipped all 8 HTML files into fastqc_results.zip"
else
  echo "Error: Not all HTML files were generated, zip not created."
fi

echo "Done Code"
