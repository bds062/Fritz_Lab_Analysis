#!/bin/sh
# Script to turn all .sam in a directory to .bam
# Requires: samtools-1.21
# BDS 4/23/2025

#SBATCH --job-name=samtobam
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=4:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/samtobam_output.txt
#SBATCH -e ./logs/samtobam_output.txt

# Set the path to samtools
SAMTOOLS="../../../programs/samtools-1.21/samtools"

for sam_file in *.sam; do
  INFILE=$sam_file
  TEMPFILE=$(echo $INFILE | sed 's/.sam/.temp/')
  OUTFILE=$(echo $INFILE | sed 's/.sam/.bam/')
  STATSFILE=$(echo $INFILE | sed 's/.sam/_stats.txt/')
  echo the sam file $INFILE was sorted and converted to the sorted bam file $OUTFILE

  $SAMTOOLS view -@ 16 -S -h -u $INFILE | \
  $SAMTOOLS sort -@ 16 -T $TEMPFILE -> $OUTFILE

  $SAMTOOLS index $OUTFILE

  $SAMTOOLS stats $OUTFILE > $STATSFILE

  $SAMTOOLS depth $OUTFILE > $(echo $INFILE | sed 's/.sam/_depth.tsv/')

  rm $sam_file
done
