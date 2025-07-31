import re

# User-defined sample size
num_samples = 21
dp_threshold = num_samples * 3

# Input/output filenames
vcf_file = "contaminated.vcf"
output_file = "contaminated_filtered3x.tsv"

with open(vcf_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith("##"):
            continue  # Skip metadata
        elif line.startswith("#CHROM"):
            header = line.strip().split('\t')
            sample_count = len(header) - 9  # Number of sample columns
            outfile.write("\t".join(header) + "\n")
        else:
            fields = line.strip().split('\t')
            info_field = fields[7]
            genotype_fields = fields[9:]

            # Count number of valid genotypes (not missing)
            valid_gt_count = sum(1 for gt in genotype_fields if not gt.startswith('./.'))

            # Extract DP value from INFO field
            match = re.search(r'DP=(\d+)', info_field)
            if not match:
                continue  # Skip if no DP found

            dp_value = int(match.group(1))

            # Apply both filters
            if valid_gt_count / num_samples >= 0.75 and dp_value >= dp_threshold:
                outfile.write("\t".join(fields) + "\n")