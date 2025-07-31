#!/usr/bin/env python3
import os
import re
import glob

# Step 1: Get all SRR IDs from trimmed SAM files
sam_files = glob.glob("trimmed*.sam")
srr_ids = sorted(set(re.findall(r'SRR\d+', ' '.join(sam_files))))

# Step 2: Parse alignment files
def parse_alignment_file(filepath):
    alignment_rates = {}
    current_srr = None
    pattern_srr = re.compile(r'Running Bowtie2 for (SRR\d+)')
    pattern_rate = re.compile(r'([\d.]+)% overall alignment rate')

    with open(filepath, 'r') as f:
        for line in f:
            srr_match = pattern_srr.search(line)
            if srr_match:
                current_srr = srr_match.group(1)
            elif current_srr:
                rate_match = pattern_rate.search(line)
                if rate_match:
                    alignment_rates[current_srr] = rate_match.group(1)
                    current_srr = None  # Reset for next SRR block
    return alignment_rates

# Step 3: Load clean and contaminated alignments
clean_alignments = parse_alignment_file("cleanAllignment.txt")
contaminated_alignments = parse_alignment_file("trimmedAllignment.txt")

# Step 4: Write to TSV
with open("alignment_summary.tsv", "w") as out:
    out.write("SRR\tClean_Alignment\tContaminated_Alignment\n")
    for srr in srr_ids:
        clean = clean_alignments.get(srr, "NA")
        contaminated = contaminated_alignments.get(srr, "NA")
        out.write(f"{srr}\t{clean}\t{contaminated}\n")

print("alignment_summary.tsv created.")