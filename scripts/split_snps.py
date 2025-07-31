import pandas as pd

print("Script started")

# ---- USER SETTINGS ----
clean_file = "clean_filtered3x.tsv"
contam_file = "contaminated_filtered3x.tsv"
output_shared = "shared_snps3x.tsv"
output_clean = "clean_unique_snps3x.tsv"
output_contam = "contaminated_unique_snps3x.tsv"
# ------------------------

# Read the clean and contaminated SNPs (skip comment lines)
clean_df = pd.read_csv(clean_file, sep='\t', comment=None)
contam_df = pd.read_csv(contam_file, sep='\t', comment=None)

# Create keys for comparison (CHROM + POS as tuple)
clean_df["key"] = list(zip(clean_df["#CHROM"], clean_df["POS"]))
contam_df["key"] = list(zip(contam_df["#CHROM"], contam_df["POS"]))

# Set lookups
clean_keys = set(clean_df["key"])
contam_keys = set(contam_df["key"])

# Get keys
shared_keys = clean_keys & contam_keys
clean_only_keys = clean_keys - contam_keys
contam_only_keys = contam_keys - clean_keys

# Filter dataframes by keys
shared_df = clean_df[clean_df["key"].isin(shared_keys)]
clean_unique_df = clean_df[clean_df["key"].isin(clean_only_keys)]
contam_unique_df = contam_df[contam_df["key"].isin(contam_only_keys)]

# Drop helper column
for df in [shared_df, clean_unique_df, contam_unique_df]:
    df.drop(columns=["key"], inplace=True)

# Write output
shared_df.to_csv(output_shared, sep='\t', index=False)
clean_unique_df.to_csv(output_clean, sep='\t', index=False)
contam_unique_df.to_csv(output_contam, sep='\t', index=False)

print("SNP separation complete:")
print("- Shared: {}".format(len(shared_df)))
print("- Clean-only: {}".format(len(clean_unique_df)))
print("- Contaminated-only: {}".format(len(contam_unique_df)))
# print(f"- Clean-only: {len(clean_unique_df)}")
# print(f"- Contaminated-only: {len(contam_unique_df)}")