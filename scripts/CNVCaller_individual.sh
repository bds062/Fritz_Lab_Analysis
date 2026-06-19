#!/usr/bin/env bash

#SBATCH --job-name=MD_CNVCaller_individual_process
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./logs/CNVCaller_individual_process_%a.out
#SBATCH -e ./logs/CNVCaller_individual_process_%a.out
#SBATCH --array=0-126

samtools="samtools"

# DO THIS:
# cp /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/Ja/window.link ./"
# cp /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/Ja/referenceDB.800 ./"

# Build array of all BAM files once per task
bam_files=(trimmed_*.bam)
bam_file="${bam_files[$SLURM_ARRAY_TASK_ID]}"

# Exit cleanly if index is out of range (handles over-estimated array size)
if [[ -z "$bam_file" || ! -f "$bam_file" ]]; then
    echo "No BAM file for task index $SLURM_ARRAY_TASK_ID — exiting."
    exit 0
fi

srrid="${bam_file#trimmed_}"
srrid="${srrid%.bam}"
echo "Task $SLURM_ARRAY_TASK_ID → $srrid"

bash /fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/programs/CNVcaller/Individual.Process.sh \
    -b "/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/Hzea_WGS_TimeSeries_RAW2/${bam_file}" \
    -h "${srrid}" \
    -d "./window.link" \
    -s Z \
    > "./logs/${srrid}_CNVCaller.out" 2>&1