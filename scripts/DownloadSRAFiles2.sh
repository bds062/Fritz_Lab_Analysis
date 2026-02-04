#!/bin/bash

module load sra

OUTDIR="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/MDTaylor"
TMPDIR="/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/MDTaylor/tmp"

mkdir -p $OUTDIR $TMPDIR

for i in $(seq 29490973 29490992); do
    SRR="SRR$i"
    echo "Downloading $SRR"

    fasterq-dump $SRR \
        --outdir $OUTDIR \
        --temp $TMPDIR \
        --split-files \
        --threads 8

done