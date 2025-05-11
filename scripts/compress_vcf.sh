../../programs/samtools-1.21/bgzip clean.vcf
../../programs/samtools-1.21/bgzip contaminated.vcf
../../programs/bcftools/tabix -p vcf clean.vcf.gz
../../programs/bcftools/tabix -p vcf contaminated.vcf.gz
mkdir isec_output
bcftools isec -p isec_output clean.vcf.gz contaminated.vcf.gz