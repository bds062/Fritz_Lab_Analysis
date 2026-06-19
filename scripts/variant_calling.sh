#!/bin/sh

#Script to make vcf for contaminated files
#Requires: BCFTools installed
#BDS 9/19/25

#SBATCH --job-name=HzeaSequencerVariantCalling
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH --qos=huge-long
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=5-00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/variantCalling_BenParams.out
#SBATCH -e ./logs/variantCalling_BenParams.out

BCFTOOLS="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/bcftools/bcftools"
REF="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/ZeaRef.fna"

ls \
/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/Hzea_WGS_TimeSeries_RAW2/*.bam \
/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/Taylor2021/2002samples/*.bam \
/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/Taylor2021/2012samples/*.bam \
/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/Taylor2021/2017samples/*.bam \
> bam_list.txt
# "$BCFTOOLS" mpileup -Ou -f "$REF" --bam-list bam_list.txt --threads 32 | "$BCFTOOLS" call -vmO v --threads 32 -o ./Hzea_WGS_TimeSeries.redo.vcf

"$BCFTOOLS" mpileup \
-f "$REF" \
--threads 32 \
-b bam_list.txt \
-a FORMAT/DP,FORMAT/AD \
-Ou | \
"$BCFTOOLS" call \
-mv \
--threads 32 \
-Ou | \
"$BCFTOOLS" view \
-v snps -m2 -M2 \
--threads 32 \
-Oz -o biallelic_snps.vcf.gz
