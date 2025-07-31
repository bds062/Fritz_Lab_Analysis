#!/bin/sh

#Script to make vcf for contaminated files
#Requires: Samtools-1.21 installed
#BDS 5/13/25

#SBATCH --job-name=Contaminated_vcf
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=1-00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu


# ../../programs/bcftools/bcftools mpileup -f ./ref_files/ZeaRef.fna -b ./clean_bam.txt -I --threads 4 -O u -o ./clean.bcf
# ../../programs/bcftools/bcftools mpileup -f ./ref_files/ZeaRef.fna -b ./contaminated_bam.txt -I --threads 4 -O u -o ./contaminated.bcf
# ../../programs/bcftools/bcftools call -vmO v ./clean.bcf -o ./clean.vcf 
../../programs/bcftools/bcftools call -vmO v ./contaminated.bcf -o ./contaminated.vcf 
# ../../programs/samtools-1.21/bgzip clean.vcf
# ../../programs/samtools-1.21/bgzip contaminated.vcf
# ../../programs/bcftools/tabix -p vcf clean.vcf.gz
# ../../programs/bcftools/tabix -p vcf contaminated.vcf.gz
# mkdir isec_output
# bcftools isec -p isec_output clean.vcf.gz contaminated.vcf.gz