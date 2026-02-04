#!/bin/bash
#SBATCH --job-name=Taylor_delly_cnv

#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=./logs/delly_cnv.log

#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb

#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu

# source activate
# conda init
# conda activate /nfshomes/bds062/miniconda3/envs/delly2/

# bwa mem -t 8 ZeaRef.fna read1.fq.gz read2.fq.gz | samtools sort -o srt.bam
# samtools index srt.bam

# dicey mappability2 srt.bam
# gunzip map.fa.gz && bgzip map.fa && samtools faidx map.fa.gz 

# samtools sort -o trimmed_SRR18709187.srt.bam trimmed_SRR18709187.bam
# samtools index trimmed_SRR18709187.srt.bam

# delly cov -g ZeaRef.fna trimmed_SRR18709187.srt.bam > out.cov.gz
# delly mg -g ZeaRef.fna trimmed_SRR18709187.srt.bam > map.txt

# genmap index -F ZeaRef.fna -I ZeaRef_index
# bedtools map -a windows.bed \
            #  -b genmap_output/genmap.bedgraph \
            #  -c 4 -o mean > ZeaRef.map



for file in *.bam; do
  # Skip if no BAM files match
  [ -e "$file" ] || continue
  # samtools view -@ 16 -b $file > $file
  # Get just the basename (no path)
  base="${file##*/}"        # sample.srt.bam

  # Strip .srt.bam and append new suffixes
  outC="${base%.srt.bam}.cov.gz"
  outO="${base%.srt.bam}.bcf"

  echo "CNV Calling $file -> $outC"
#   samtools sort -o "$out" "$bam"
  # delly cnv \
  #  -a \
  #  -i 4000000000 \
  #  -g /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/delly2/ZeaRef.fna \
  #  -m /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/delly2/map.fa.gz \
  #  -c ./"$outC" \
  #  -o ./"$outO" \
  #  "$file"
  /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/delly/src/delly cnv \
   -a -n \
   -g /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/delly2/ZeaRef.fna \
   -m /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/delly2/map.fa.gz \
   -c ./"$outC" \
   -o ./"$outO" \
   "$file"
done