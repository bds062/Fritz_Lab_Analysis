#!/bin/sh

# Script to run BBMap on all previously identified contaminated files
# Requires: BBMap Installation
# BDS 1/3/25

#SBATCH --job-name=Contamination_mapping
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=540:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./bbmap_output.txt
#SBATCH -e ./bbmap_output.txt

# Path to the BBMap program
bbmap_program="../../programs/BBMap/bbmap/bbmap.sh"

# Iterate through each .fastq file in the directory
for file in trimmed_*.fastq; do
  # Extract the base filename without the "trimmed_" prefix
  base_name=${file#trimmed_}
  
  # Construct the output file name for the first mapping command
  output_file="mapped_${base_name}"
  
  # Run the BBMap mapping command
  echo "Processing file: $file"
  $bbmap_program in="$file" out="$output_file"
  
  # Check if the mapping command succeeded
  if [ $? -eq 0 ]; then
    echo "Successfully processed: $file -> $output_file"
    
    # Run the second BBMap command to generate coverage statistics
    $bbmap_program in="$file" \
      covstats="covstats_${base_name}.txt" \
      covhist="covhist_${base_name}.txt" \
      basecov="basecov_${base_name}.txt" \
      bincov="bincov_${base_name}.txt"
    
    # Check if the second command succeeded
    if [ $? -eq 0 ]; then
      echo "Coverage stats generated for: $file"
    else
      echo "Error generating coverage stats for: $file" >&2
    fi
  else
    echo "Error processing: $file" >&2
  fi
done

echo "All files processed."
