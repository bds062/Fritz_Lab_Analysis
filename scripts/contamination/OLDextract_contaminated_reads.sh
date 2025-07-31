#!/bin/sh

#SBATCH --job-name=ExtractContaminatedReads
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=05:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./extract_contaminated_reads.txt
#SBATCH -e ./extract_contaminated_reads.txt

# Path to filterbyname.sh
FILTERBYNAME="../../programs/BBMap/bbmap/filterbyname.sh"

# Directory containing FASTQ files (change if needed)
FASTQ_DIR="./"

# Loop through paired FASTQ files
for file in ${FASTQ_DIR}/trimmed_*_1_paired.fastq; do
    # Extract the SRR ID
    base_name=$(basename "$file" | sed -E 's/trimmed_(SRR[0-9]+)_1_paired.fastq/\1/')
    
    # Define the matching read 2 file
    read2_file="${FASTQ_DIR}/trimmed_${base_name}_2_paired.fastq"
    
    # Define clean file names
    clean1="clean_${base_name}_1.fastq"
    clean2="clean_${base_name}_2.fastq"
    
    # Check if all required files exist
    if [[ -f "$read2_file" && -f "$clean1" && -f "$clean2" ]]; then
        echo "Running filterbyname.sh for $base_name..."
        $FILTERBYNAME in1="$file" in2="$read2_file" out="contaminated_reads_${base_name}.fastq" \
    names="${clean1},${clean2}" include=t -Xmx16g
    else
        echo "Skipping $base_name: Missing one or more required files."
    fi
done

echo "filterbyname.sh processing complete."