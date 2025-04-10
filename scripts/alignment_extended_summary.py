#!/usr/bin/env python3
import re
import csv

def extract_data(filename, prefix):
    data = {}
    with open(filename, 'r') as f:
        lines = f.readlines()

    current_srr = None
    for i, line in enumerate(lines):
        srr_match = re.match(r"Running Bowtie2 for (SRR\d+)\.\.\.", line)
        if srr_match:
            current_srr = srr_match.group(1)
            data[current_srr] = {}
            continue

        if current_srr:
            if "reads; of these:" in line:
                data[current_srr][f"{prefix}TotalReads"] = line.strip().split()[0]
            elif "aligned concordantly 0 times" in line:
                data[current_srr][f"{prefix}ConcordantlyZero"] = line.strip().split()[0]
            elif "aligned concordantly exactly 1 time" in line:
                data[current_srr][f"{prefix}ConcordantlyOne"] = line.strip().split()[0]
            elif "aligned concordantly >1 times" in line:
                data[current_srr][f"{prefix}ConcordantlyGOne"] = line.strip().split()[0]
            elif "aligned discordantly 1 time" in line:
                data[current_srr][f"{prefix}DisconcordantlyOne"] = line.strip().split()[0]
            elif "mates make up the pairs" in line:
                data[current_srr][f"{prefix}ExactMateReads"] = line.strip().split()[0]
            elif "aligned exactly 1 time" in line and f"{prefix}ExactOne" not in data[current_srr]:
                data[current_srr][f"{prefix}ExactOne"] = line.strip().split()[0]
            elif "aligned >1 times" in line and f"{prefix}ExactG1" not in data[current_srr]:
                data[current_srr][f"{prefix}ExactG1"] = line.strip().split()[0]
            elif "overall alignment rate" in line:
                data[current_srr][f"{prefix}AlignmentRate"] = line.strip().split('%')[0]
    return data

# Extract from both files
co_data = extract_data("trimmedAllignment.txt", "Co")
cl_data = extract_data("cleanAllignment.txt", "Cl")

# Combine all SRRs
all_srrs = set(co_data.keys()).union(set(cl_data.keys()))

# Headers
fields = ["SRR"]
for prefix in ["Co", "Cl"]:
    fields.extend([
        f"{prefix}TotalReads", f"{prefix}ConcordantlyZero", f"{prefix}ConcordantlyOne",
        f"{prefix}ConcordantlyGOne", f"{prefix}DisconcordantlyOne", f"{prefix}ExactMateReads",
        f"{prefix}ExactOne", f"{prefix}ExactG1", f"{prefix}AlignmentRate"
    ])

# Write to TSV
with open("alignment_extended_summary.tsv", "w", newline='') as out_f:
    writer = csv.DictWriter(out_f, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    for srr in sorted(all_srrs):
        row = {"SRR": srr}
        if srr in co_data:
            row.update(co_data[srr])
        else:
            row.update({f: "NA" for f in fields if f.startswith("Co")})
        if srr in cl_data:
            row.update(cl_data[srr])
        else:
            row.update({f: "NA" for f in fields if f.startswith("Cl")})
        writer.writerow(row)
