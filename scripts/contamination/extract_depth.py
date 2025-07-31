# User-defined file paths
input_file = "cleanvcf_stats.txt"
output_file = "cleandepth.tsv"

previous_line = ""
header = []
data_rows = []

with open(input_file, "r") as infile:
    for line in infile:
        line = line.strip()

        # First data line
        if line.startswith("DP") and not header:
            # Parse the line before as header
            if previous_line.startswith("# DP"):
                raw_header_fields = previous_line.split("\t")[1:]  # skip '# DP '
                header = [field.split("]", 1)[-1].strip() if "]" in field else field for field in raw_header_fields]
            else:
                raise ValueError("Expected header line before first DP entry was not found.")
            
            # Process the current data line
            data_rows.append(line.split("\t")[1:])  # skip 'DP'

        elif line.startswith("DP") and header:
            data_rows.append(line.split("\t")[1:])  # skip 'DP'

        previous_line = line  # Keep track of last line for next loop

# Write to output file
with open(output_file, "w") as outfile:
    outfile.write("\t".join(header) + "\n")
    for row in data_rows:
        outfile.write("\t".join(row) + "\n")

print(f"✅ Data written to {output_file}")