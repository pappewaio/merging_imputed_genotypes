#!/usr/bin/env bash

set -euo pipefail

test_script="merge_vcfs"
initial_dir=$(pwd)"/${test_script}"
curr_case=""

mkdir "${initial_dir}"
cd "${initial_dir}"

#=================================================================================
# Helpers
#=================================================================================

function _setup {
  mkdir "${1}"
  cd "${1}"
  curr_case="${1}"
}

function _check_results {
  obs=$1
  exp=$2
  if ! diff ${obs} ${exp} &> ./difference; then
    echo "- [FAIL] ${curr_case}"
    cat ./difference 
    exit 1
  fi

}

function _run_script {

  "${test_script}.sh" ./input1.vcf.gz ./input2.vcf.gz output_merged.vcf.gz output_excluded_1.vcf.gz output_excluded_2.vcf.gz

  _check_results <(grep -v "##" <(zcat output_merged.vcf)) ./expected-overlap.tsv
  _check_results <(grep -v "##" <(zcat output_excluded_1.vcf.gz)) ./expected-excluded-1.tsv
  _check_results <(grep -v "##" <(zcat output_excluded_2.vcf.gz)) ./expected-excluded-2.tsv

  echo "- [OK] ${curr_case}"

  cd "${initial_dir}"
}

echo ">> Test ${test_script}"

#=================================================================================
# Cases
#=================================================================================

#---------------------------------------------------------------------------------
# Check that the final output is what we think it is

_setup "merge two vcf files based on overlapping variants"

cat <<EOF | bcftools sort | bgzip -c > ./input1.vcf.gz
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=11032019_15h52m43s
##source=IGSRpipeline
##contig=<ID=10>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the East Asian population">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the European population">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the African population">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the American population">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the South Asian population">
##INFO=<ID=VT,Number=1,Type=String,Description="Variant Type">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=EX_TARGET,Number=0,Type=Flag,Description="indicates whether a variant is within the exon pull down target boundaries">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
10	100157733	snp1	C	T	.	PASS	AC=415;AN=5096;DP=20612;AF=0.08;EAS_AF=0.04;EUR_AF=0.07;AFR_AF=0.01;AMR_AF=0.05;SAS_AF=0.26;VT=SNP;NS=2548	GT	0/1	1/1
10	101966771	snp2	T	C	.	PASS	AC=459;AN=5096;DP=20740;AF=0.09;EAS_AF=0.07;EUR_AF=0.17;AFR_AF=0.01;AMR_AF=0.17;SAS_AF=0.07;VT=SNP;NS=2548	GT	0/1	0/1
10	102814179	snp3	T	C	.	PASS	AC=2793;AN=5096;DP=20193;AF=0.55;EAS_AF=0.4;EUR_AF=0.59;AFR_AF=0.54;AMR_AF=0.58;SAS_AF=0.64;VT=SNP;NS=2548	GT	1/1	0/1
10	104355789	snp4	T	C	.	PASS	AC=763;AN=5096;DP=23195;AF=0.15;EAS_AF=0.08;EUR_AF=0.28;AFR_AF=0.06;AMR_AF=0.18;SAS_AF=0.2;VT=SNP;NS=2548	GT	0/1	1/1
10	10574522	snp25	T	C	.	PASS	AC=1180;AN=5096;DP=20409;AF=0.23;EAS_AF=0.33;EUR_AF=0.18;AFR_AF=0.24;AMR_AF=0.11;SAS_AF=0.25;VT=SNP;NS=2548	GT	1/1	0/1
10	105905360	snp6	A	G	.	PASS	AC=3254;AN=5096;DP=16605;AF=0.64;EAS_AF=0.7;EUR_AF=0.46;AFR_AF=0.88;AMR_AF=0.59;SAS_AF=0.48;VT=SNP;NS=2548	GT	0/1	0/1
10	106322897	snp7	C	T	.	PASS	AC=695;AN=5096;DP=19936;AF=0.14;EAS_AF=0;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.08;SAS_AF=0.1;VT=SNP;NS=2548	GT	1/1	1/1
10	106371793	snp8	C	A	.	PASS	AC=906;AN=5096;DP=19871;AF=0.18;EAS_AF=0.13;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.21;SAS_AF=0.09;VT=SNP;NS=2548	GT	0/1	1/1
10	106524727	snp23	G	A	.	PASS	AC=2570;AN=5096;DP=18067;AF=0.5;EAS_AF=0.44;EUR_AF=0.37;AFR_AF=0.71;AMR_AF=0.35;SAS_AF=0.54;VT=SNP;NS=2548	GT	1/1	0/1
EOF

tabix -p vcf ./input1.vcf.gz


cat <<EOF | bcftools sort | bgzip -c > ./input2.vcf.gz
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=11032019_15h52m43s
##source=IGSRpipeline
##contig=<ID=10>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the East Asian population">
##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the European population">
##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the African population">
##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the American population">
##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the South Asian population">
##INFO=<ID=VT,Number=1,Type=String,Description="Variant Type">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=EX_TARGET,Number=0,Type=Flag,Description="indicates whether a variant is within the exon pull down target boundaries">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample3	Sample4
10	100157733	snp1	C	T	.	PASS	AC=415;AN=5096;DP=20612;AF=0.08;EAS_AF=0.04;EUR_AF=0.07;AFR_AF=0.01;AMR_AF=0.05;SAS_AF=0.26;VT=SNP;NS=2548	GT	0/1	1/1
10	101966771	snp2	T	C	.	PASS	AC=459;AN=5096;DP=20740;AF=0.09;EAS_AF=0.07;EUR_AF=0.17;AFR_AF=0.01;AMR_AF=0.17;SAS_AF=0.07;VT=SNP;NS=2548	GT	0/0	0/1
10	102814179	snp3	T	C	.	PASS	AC=2793;AN=5096;DP=20193;AF=0.55;EAS_AF=0.4;EUR_AF=0.59;AFR_AF=0.54;AMR_AF=0.58;SAS_AF=0.64;VT=SNP;NS=2548	GT	1/1	1/1
10	104355789	snp4	T	C	.	PASS	AC=763;AN=5096;DP=23195;AF=0.15;EAS_AF=0.08;EUR_AF=0.28;AFR_AF=0.06;AMR_AF=0.18;SAS_AF=0.2;VT=SNP;NS=2548	GT	0/1	0/0
10	10574522	snp5	T	C	.	PASS	AC=1180;AN=5096;DP=20409;AF=0.23;EAS_AF=0.33;EUR_AF=0.18;AFR_AF=0.24;AMR_AF=0.11;SAS_AF=0.25;VT=SNP;NS=2548	GT	1/0	1/1
10	105905360	snp6	A	G	.	PASS	AC=3254;AN=5096;DP=16605;AF=0.64;EAS_AF=0.7;EUR_AF=0.46;AFR_AF=0.88;AMR_AF=0.59;SAS_AF=0.48;VT=SNP;NS=2548	GT	1/1	0/1
10	106322897	snp7	C	T	.	PASS	AC=695;AN=5096;DP=19936;AF=0.14;EAS_AF=0;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.08;SAS_AF=0.1;VT=SNP;NS=2548	GT	0/1	1/0
10	106371793	snp11	C	A	.	PASS	AC=906;AN=5096;DP=19871;AF=0.18;EAS_AF=0.13;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.21;SAS_AF=0.09;VT=SNP;NS=2548	GT	1/1	0/0
10	106524727	snp10	G	A	.	PASS	AC=2570;AN=5096;DP=18067;AF=0.5;EAS_AF=0.44;EUR_AF=0.37;AFR_AF=0.71;AMR_AF=0.35;SAS_AF=0.54;VT=SNP;NS=2548	GT	0/0	1/1
EOF

tabix -p vcf ./input2.vcf.gz

cat <<EOF > ./expected-overlap.tsv
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3	Sample4
10	100157733	snp1	C	T	.	PASS	VT=SNP;NS=2548;DP=41224;AF=0.08;EAS_AF=0.04;EUR_AF=0.07;AFR_AF=0.01;AMR_AF=0.05;SAS_AF=0.26;AN=8;AC=6	GT	0/1	1/1	0/1	1/1
10	101966771	snp2	T	C	.	PASS	VT=SNP;NS=2548;DP=41480;AF=0.09;EAS_AF=0.07;EUR_AF=0.17;AFR_AF=0.01;AMR_AF=0.17;SAS_AF=0.07;AN=8;AC=3	GT	0/1	0/1	0/0	0/1
10	102814179	snp3	T	C	.	PASS	VT=SNP;NS=2548;DP=40386;AF=0.55;EAS_AF=0.4;EUR_AF=0.59;AFR_AF=0.54;AMR_AF=0.58;SAS_AF=0.64;AN=8;AC=7	GT	1/1	0/1	1/1	1/1
10	104355789	snp4	T	C	.	PASS	VT=SNP;NS=2548;DP=46390;AF=0.15;EAS_AF=0.08;EUR_AF=0.28;AFR_AF=0.06;AMR_AF=0.18;SAS_AF=0.2;AN=8;AC=4	GT	0/1	1/1	0/1	0/0
10	105905360	snp6	A	G	.	PASS	VT=SNP;NS=2548;DP=33210;AF=0.64;EAS_AF=0.7;EUR_AF=0.46;AFR_AF=0.88;AMR_AF=0.59;SAS_AF=0.48;AN=8;AC=5	GT	0/1	0/1	1/1	0/1
10	106322897	snp7	C	T	.	PASS	VT=SNP;NS=2548;DP=39872;AF=0.14;EAS_AF=0;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.08;SAS_AF=0.1;AN=8;AC=6	GT	1/1	1/1	0/1	1/0
EOF

cat <<EOF > ./expected-excluded-1.tsv
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
10	10574522	snp25	T	C	.	PASS	AC=1180;AN=5096;DP=20409;AF=0.23;EAS_AF=0.33;EUR_AF=0.18;AFR_AF=0.24;AMR_AF=0.11;SAS_AF=0.25;VT=SNP;NS=2548	GT	1/1	0/1
10	106371793	snp8	C	A	.	PASS	AC=906;AN=5096;DP=19871;AF=0.18;EAS_AF=0.13;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.21;SAS_AF=0.09;VT=SNP;NS=2548	GT	0/1	1/1
10	106524727	snp23	G	A	.	PASS	AC=2570;AN=5096;DP=18067;AF=0.5;EAS_AF=0.44;EUR_AF=0.37;AFR_AF=0.71;AMR_AF=0.35;SAS_AF=0.54;VT=SNP;NS=2548	GT	1/1	0/1
EOF

cat <<EOF > ./expected-excluded-2.tsv
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample3	Sample4
10	10574522	snp5	T	C	.	PASS	AC=1180;AN=5096;DP=20409;AF=0.23;EAS_AF=0.33;EUR_AF=0.18;AFR_AF=0.24;AMR_AF=0.11;SAS_AF=0.25;VT=SNP;NS=2548	GT	1/0	1/1
10	106371793	snp11	C	A	.	PASS	AC=906;AN=5096;DP=19871;AF=0.18;EAS_AF=0.13;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.21;SAS_AF=0.09;VT=SNP;NS=2548	GT	1/1	0/0
10	106524727	snp10	G	A	.	PASS	AC=2570;AN=5096;DP=18067;AF=0.5;EAS_AF=0.44;EUR_AF=0.37;AFR_AF=0.71;AMR_AF=0.35;SAS_AF=0.54;VT=SNP;NS=2548	GT	0/0	1/1
EOF

_run_script

#---------------------------------------------------------------------------------
# Next case

#_setup "valid_rows_missing_afreq"
#
#cat <<EOF > ./acor.tsv
#0	A1	A2	CHRPOS	RSID	EffectAllele	OtherAllele	EMOD
#1	A	G	12:126406434	rs1000000	G	A	-1
#EOF
#
#cat <<EOF > ./stat.tsv
#0	B	SE	Z	P
#1	-0.0143	0.0156	-0.916667	0.3604
#EOF
#
#_run_script
