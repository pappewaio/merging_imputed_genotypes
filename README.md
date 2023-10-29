# Merging imputation

## Todo
It would be great to prepare a few .ppt slides describing the merge and conversion to plink2.  Counts before, counts after, any issues, the basic workflow, basic commands.

Find the MAF and INFO score data.
We also want to restrict to just EUR subset for the test analysis 

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


## Run merger
```
srun --mem=8g --ntasks 1 --cpus-per-task 2 --time=9:00:00 --account ibp_pipeline_cleansumstats --pty /bin/bash
conda activate merging_imputed_genotypes
nextflow run main.nf
```


