#!/bin/sh

#Script to run Bowtie2 and allign all "clean" files in a directory
#Requires: Bowtie2 Installation
#BDS 5/27/25

#SBATCH --job-name=BactoZea
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=1-00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./bacteriaToZeaAlignment.txt
#SBATCH -e ./bacteriaToZeaAlignment.txt

BBMAP="../../../programs/BBMap/bbmap/bbmap.sh"
BOWTIE="../../../programs/bowtie2/bowtie2"
TAGS="-t -x ../zea"

# Path to reference genome
REF="ZeaRef.fna"

# List of bacterial genome files
BACTERIAL_GENOMES=(
    "BacillusCereus.fna"
    "EnterobacterCloacae.fna"
    "HeliothisZeaNudivirus.fna"
    "MorganellaMorganii.fna"
)

# Loop through bacterial genomes and align to Zea reference
for GENOME in "${BACTERIAL_GENOMES[@]}"; do
    if [[ -f "$GENOME" ]]; then
        BASE=$(basename "$GENOME" .fna)
        echo "Aligning $BASE to Zea reference..."
        $BBMAP in="$GENOME" ref="$REF" out="${BASE}_to_zea.sam"
    else
        echo "File $GENOME not found. Skipping..."
    fi
done

echo "BBMap alignments complete."