#!/bin/bash

# Script to run FastQC on any number of mapped SRR .fastq files
# Requires: FastQC installed
# BDS 1/3/2025

# Check if at least one parameter is passed
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <SRR_ID1> [<SRR_ID2> ... <SRR_IDN>]"
  exit 1
fi

# Loop through the provided SRR IDs and process each one
for srr in "$@"; do
  # Check if both mapped files (_1.fastq and _2.fastq) exist for each SRR ID
  if [[ -f "mapped_${srr}_1.fastq" && -f "mapped_${srr}_2.fastq" ]]; then
    echo "Running fastqc on mapped_${srr}_1.fastq and mapped_${srr}_2.fastq"
    ../../programs/FastQC/fastqc "mapped_${srr}_1.fastq" "mapped_${srr}_2.fastq"
  else
    echo "Warning: Files for mapped_${srr} not found or incomplete."
  fi
done

# Gather the generated *_fastqc.html files for zipping
html_files=()

for srr in "$@"; do
  if [[ -f "mapped_${srr}_1_fastqc.html" && -f "mapped_${srr}_2_fastqc.html" ]]; then
    html_files+=("mapped_${srr}_1_fastqc.html" "mapped_${srr}_2_fastqc.html")
  else
    echo "Warning: fastqc HTML files for mapped_${srr} not found."
  fi
done

# Zip the HTML files if any are found
if [ "${#html_files[@]}" -gt 0 ]; then
  zip fastqc_results.zip "${html_files[@]}"
  echo "Zipped all FastQC HTML files into fastqc_results.zip"
else
  echo "Error: No HTML files generated, zip not created."
fi

echo "Done Code"
