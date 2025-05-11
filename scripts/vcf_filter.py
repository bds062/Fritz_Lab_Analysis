# Set your input VCF file name here
vcf_file = "clean.vcf"
output_file = "clean_filtered.tsv"

with open(vcf_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith("##"):
            continue  # Skip metadata
        elif line.startswith("#CHROM"):
            header = line.strip().split('\t')
            sample_count = len(header) - 9  # Number of sample columns
            outfile.write("\t".join(header) + "\n")  # Write header to output
        else:
            fields = line.strip().split('\t')
            genotypes = fields[9:]  # Sample columns
            valid_gt_count = sum(1 for gt in genotypes if not gt.startswith('./.'))

            if valid_gt_count / sample_count >= 0.75:
                outfile.write("\t".join(fields) + "\n")