#!/bin/sh

#Script to plot mapq scores for contaminated and clean files
#Requires: Clean and Contaminated sam files
#BDS 6/23/25

#SBATCH --job-name=MapqSummary
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --qos=high
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT
#SBATCH --mail-user=bds062@terpmail.umd.edu
#SBATCH -o ./mapqPlot.txt
#SBATCH -e ./mapqPlot.txt

source activate base
module load Python3/3.13.2 
pip install pandas
python3 ../../Fritz_Lab_Analysis/scripts/mapq_summary.py