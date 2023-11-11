# generate pgen
Use the created VCF and generate pgen file using plink2

```
vcfdir="test/example_data/genotypes/vcf"
out_plink2_ds="test/example_data/genotypes/pgen"
mkdir -p ${out_plink2_ds}

# Run
for vcf in ${vcfdir}/*vcf.gz; do
  bvcf=$(basename $vcf)
  echo "running: $vcf"
  plink2 --vcf ${vcf} dosage=DS --double-id --memory 2500 --threads 5 --out ${out_plink2_ds}/${bvcf%.vcf.gz}
done

```


