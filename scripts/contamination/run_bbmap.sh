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
#SBATCH --begin=now
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./bbmap_output_MM.txt
#SBATCH -e ./bbmap_output_MM.txt

# Path to the BBMap program
BBMAP="../../programs/BBMap/bbmap/bbmap.sh"

FASTQ_DIR="."
# Iterate through each .fastq file in the directory
# Species tag for the output file (edit this before running)
PREV="HZN"
SPECIES="MM"

../../programs/BBMap/bbmap/bbmap.sh ref=ref_files/MorganellaMorganii.fna

# Loop through paired FASTQ files
for file in ${FASTQ_DIR}/trimmed_*_1_paired.fastq; do
    # Extract the SRR ID
    base_name=$(basename "$file" | sed -E 's/trimmed_(SRR[0-9]+)_1_paired.fastq/\1/')
    
    # Define the matching read 2 file
    read2_file="${FASTQ_DIR}/trimmed_${base_name}_2_paired.fastq"
    
    # Check if read 2 file exists
    if [[ -f "$read2_file" ]]; then
        echo "Running BBMap for $base_name..."
        # $BBMAP in="$file" outm="mapped_${SPECIES}_${base_name}_1.fastq" outu="unmapped_${SPECIES}_${base_name}_1.fastq"
        # $BBMAP in="$read2_file" outm="mapped_${SPECIES}_${base_name}_2.fastq" outu="unmapped_${SPECIES}_${base_name}_2.fastq"
        # $BBMAP in1="$file" in2="$read2_file" out1m="mapped_${SPECIES}_${base_name}_1.fastq" out1u="unmapped_${SPECIES}_${base_name}_1.fastq" out2m="mapped_${SPECIES}_${base_name}_2.fastq" out2u="unmapped_${SPECIES}_${base_name}_2.fastq"
        $BBMAP in="unmapped_${PREV}_${base_name}.fastq" outm="mapped_${SPECIES}_${base_name}.fastq" outu="unmapped_${SPECIES}_${base_name}.fastq"
    else
        echo "Skipping $read2_file: Read 2 file not found."
    fi
done

echo "BBMap processing complete."
