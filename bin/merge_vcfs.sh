#!/bin/bash

# Check for correct number of arguments
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <file1.vcf.gz> <file2.vcf.gz> <output_merged.vcf.gz> <output_excluded_1.vcf.gz> <output_excluded_2.vcf.gz>"
    exit 1
fi

# Input VCF files
VCF1="$1"
VCF2="$2"

# Output files
OUTPUT_MERGED="$3"
OUTPUT_EXCLUDED1="$4"
OUTPUT_EXCLUDED2="$5"

# Identify overlapping and unique variants
bcftools query -f '%ID\n' "$VCF1" | sort > ID1
bcftools query -f '%ID\n' "$VCF2" | sort > ID2
comm -12 ID1 ID2 > overlapping_ids.txt
comm -23 ID1 ID2 > vcf1_unique_ids.txt
comm -13 ID1 ID2 > vcf2_unique_ids.txt

# filter both files
bcftools view -i 'ID=@overlapping_ids.txt' -Oz "$VCF1" > vcf1_common.vcf.gz
bcftools view -i 'ID=@overlapping_ids.txt' -Oz "$VCF2" > vcf2_common.vcf.gz
bcftools view -i 'ID=@vcf1_unique_ids.txt' -Oz "$VCF1" > ${OUTPUT_EXCLUDED1}
bcftools view -i 'ID=@vcf2_unique_ids.txt' -Oz "$VCF2" > ${OUTPUT_EXCLUDED2}

# Merge the overlapping variants
tabix -p vcf vcf1_common.vcf.gz
tabix -p vcf vcf2_common.vcf.gz
bcftools merge -m id -Oz -o "$OUTPUT_MERGED" vcf1_common.vcf.gz vcf2_common.vcf.gz

# Make tabix index
tabix -p vcf "$OUTPUT_MERGED"


