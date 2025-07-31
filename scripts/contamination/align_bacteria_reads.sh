#!/bin/bash

#SBATCH --job-name=Bowtie2ContaminatedReads
#SBATCH -c 16
#SBATCH --mem=30G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/contam_bowtie2_all.out
#SBATCH -e ./logs/contam_bowtie2_all.err

# Paths
BOWTIE="../../programs/bowtie2/bowtie2"
TAGS="-t -x zea"

# Ensure logs directory exists
mkdir -p logs

# Loop through each ID in ids.txt
while read -r ID; do
    R1="contaminated_${ID}_1.fastq"
    R2="contaminated_${ID}_2.fastq"
    OUT="contaminated_${ID}.sam"

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "Running Bowtie2 for $ID..."
        $BOWTIE $TAGS -1 "$R1" -2 "$R2" --very-sensitive -S "$OUT"
    else
        echo "Skipping $ID: Paired FASTQ files not found."
    fi
done < ids.txt

echo "All Bowtie2 alignments completed."