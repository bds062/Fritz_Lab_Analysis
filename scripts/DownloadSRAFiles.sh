#! /bin/bash

#I used this file to download all North et al. files from the SRA on 9/19/2024

module load sra

while read SRRlist; do
  echo "$SRRlist"
  fasterq-dump "$SRRlist"
done <SRR_Numbers.txt


