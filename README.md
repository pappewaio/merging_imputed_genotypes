# Merging imputation

## Todo
It would be great to prepare a few .ppt slides describing the merge and conversion to plink2.  Counts before, counts after, any issues, the basic workflow, basic commands.

Find the MAF and INFO score data.
We also want to restrict to just EUR subset for the test analysis 

output also the separate common ID sets, which can be useful

## Installation

```
# install nextflow using mamba (requires conda/mamba)
mamba create -n merging_imputed_genotypes --channel bioconda \
  nextflow==23.10.0 \
  plink=1.90b6.21 \
  plink2=2.00a5 \
  bcftools=1.18 \
  tabix=1.11
  

mamba activate merging_imputed_genotypes
```


## Run merger vcf
```
srun --mem=8g --ntasks 1 --cpus-per-task 2 --time=9:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
conda activate merging_imputed_genotypes
nextflow run main.nf --type vcf
```

Use a wrapper to process chromosomes in parallel on a sbatch system
```
for chr in $(seq 1 21);do
  sleep 0.1;
  sbatch --mem=40g --ntasks 1 --cpus-per-task 30 --time=9:00:00 --account ibp_pipeline_cleansumstats --job-name="cl_${chr}" --output="clean_${chr}.out" --error="clean_${chr}.err" --wrap="
  echo ${chr} 
  date
  conda activate merging_imputed_genotypes
  nextflow run main.nf --type vcf --mergevcf.set1 test/example_data/assets/vcf/set1
  echo ${chr}
  date
  "
done
```


## Run merger pgen
```
srun --mem=8g --ntasks 1 --cpus-per-task 2 --time=9:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
conda activate merging_imputed_genotypes
nextflow run main.nf --type pgen
```


