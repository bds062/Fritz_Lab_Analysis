#!/bin/sh

#Script to run Bowtie2 and allign all "clean" files in a directory
#Requires: Bowtie2 Installation
#BDS 3/21/25

#SBATCH --job-name=Bowtie2Cleans
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=540:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./cleanAllignment.txt
#SBATCH -e ./cleanAllignment.txt

BOWTIE="../../programs/bowtie2/bowtie2"
TAGS="-t -x zea"

for file1 in ./clean_*_1.fastq; do
    # Extract the SRR ID
    base_name=$(basename "$file1" | sed -E 's/clean_(SRR[0-9]+)_1.fastq/\1/')
    file2="clean_${base_name}_2.fastq"

    if [[ -f "$file2" ]]; then
        echo "Running Bowtie2 for $base_name..."
        $BOWTIE $TAGS -1 "$file1" -2 "$file2" --very-sensitive -S "${base_name}.sam"
        else
        echo "Skipping ${file1}: Read 2 file not found."
    fi
done