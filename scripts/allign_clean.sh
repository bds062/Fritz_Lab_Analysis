#!/bin/sh

#Script to run Bowtie2 and allign all "clean" files in a directory
#Requires: Bowtie2 Installation
#BDS 3/21/25

#SBATCH --job-name=Bowtie2Test
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=540:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./Ja_testBowtie.txt
#SBATCH -e ./Ja_testBowtie.txt

../../programs/bowtie2/bowtie2 -t -x zea -1 clean_SRR18709356_1.fastq -2 clean_SRR18709356_2.fastq --very-sensitive -S test56.sam
../../programs/bowtie2/bowtie2 -t -x zea -1 clean_SRR18709355_1.fastq -2 clean_SRR18709355_2.fastq --very-sensitive -S test55.sam
../../programs/bowtie2/bowtie2 -t -x zea -1 clean_SRR18709348_1.fastq -2 clean_SRR18709348_2.fastq --very-sensitive -S test48.sam

BOWTIE="../../programs/bowtie2/bowtie2"
TAGS="-t -x zea"

for file1 in ${FASTQ_DIR}/clean_*_1.fastq; do
    # Extract the SRR ID
    base_name=$(basename "$file" | sed -E 's/clean_(SRR[0-9]+)_1.fastq/\1/')
    file2="clean_${base_name}_2.fastq"

    if [[ -f "$read2_file" ]]; then
        echo "Running Bowtie2 for $base_name..."
        $BOWTIE $TAGS -1 $file1 -2 $file2 --very-sensitive -S $base_name.sam
        else
        echo "Skipping $file2: Read 2 file not found."
    fi