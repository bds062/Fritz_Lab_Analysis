#!/bin/sh
# Script to genotype a group of .bam files
# Requires: bcftools
# BDS 4/23/2025

#SBATCH --job-name=bcftools
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=1-00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./bcftools_output.txt
#SBATCH -e ./bcftools_output.txt

# Set the path to bcftools
BCFTOOLS="../../programs/bcftools/bcftools"

$BCFTOOLS mpileup -f ./ref_files/ZeaRef.fna -b ./clean_bam.txt -I -O u -o ./clean.bcf
# $BCFTOOLS mpileup -f ./ref_files/ZeaRef.fna -b ./contaminated_bam.txt -I -O u -o ./contaminated.bcf

# $BCFTOOLS call -vmO v ./clean.bcf -o ./clean.vcf 
# $BCFTOOLS call -vmO v ./contaminated.bcf -o ./contaminated.vcf 