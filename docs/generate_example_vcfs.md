# Generate example VCFs
The starting material was real data, but now there is no link between ids, chromosomes or positions anymore.

## CHR 10
```
cat <<EOF | bcftools sort | bgzip -c > ./input1_chr10.vcf.gz
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

tabix -p vcf ./input1_chr10.vcf.gz


cat <<EOF | bcftools sort | bgzip -c > ./input2_chr10.vcf.gz
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

tabix -p vcf ./input2_chr10.vcf.gz


```

## CHR 4
```
cat <<EOF | bcftools sort | bgzip -c > ./input1_chr4.vcf.gz
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=11032019_15h52m43s
##source=IGSRpipeline
##contig=<ID=4>
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
4	100157733	snp1	C	T	.	PASS	AC=415;AN=5096;DP=20612;AF=0.08;EAS_AF=0.04;EUR_AF=0.07;AFR_AF=0.01;AMR_AF=0.05;SAS_AF=0.26;VT=SNP;NS=2548	GT	0/1	1/1
4	101966771	snp2	T	C	.	PASS	AC=459;AN=5096;DP=20740;AF=0.09;EAS_AF=0.07;EUR_AF=0.17;AFR_AF=0.01;AMR_AF=0.17;SAS_AF=0.07;VT=SNP;NS=2548	GT	0/1	0/1
4	102814179	snp3	T	C	.	PASS	AC=2793;AN=5096;DP=20193;AF=0.55;EAS_AF=0.4;EUR_AF=0.59;AFR_AF=0.54;AMR_AF=0.58;SAS_AF=0.64;VT=SNP;NS=2548	GT	1/1	0/1
4	104355789	snp4	T	C	.	PASS	AC=763;AN=5096;DP=23195;AF=0.15;EAS_AF=0.08;EUR_AF=0.28;AFR_AF=0.06;AMR_AF=0.18;SAS_AF=0.2;VT=SNP;NS=2548	GT	0/1	1/1
4	10574522	snp25	T	C	.	PASS	AC=1180;AN=5096;DP=20409;AF=0.23;EAS_AF=0.33;EUR_AF=0.18;AFR_AF=0.24;AMR_AF=0.11;SAS_AF=0.25;VT=SNP;NS=2548	GT	1/1	0/1
4	105905360	snp6	A	G	.	PASS	AC=3254;AN=5096;DP=16605;AF=0.64;EAS_AF=0.7;EUR_AF=0.46;AFR_AF=0.88;AMR_AF=0.59;SAS_AF=0.48;VT=SNP;NS=2548	GT	0/1	0/1
4	106322897	snp7	C	T	.	PASS	AC=695;AN=5096;DP=19936;AF=0.14;EAS_AF=0;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.08;SAS_AF=0.1;VT=SNP;NS=2548	GT	1/1	1/1
4	106371793	snp8	C	A	.	PASS	AC=906;AN=5096;DP=19871;AF=0.18;EAS_AF=0.13;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.21;SAS_AF=0.09;VT=SNP;NS=2548	GT	0/1	1/1
4	106524727	snp23	G	A	.	PASS	AC=2570;AN=5096;DP=18067;AF=0.5;EAS_AF=0.44;EUR_AF=0.37;AFR_AF=0.71;AMR_AF=0.35;SAS_AF=0.54;VT=SNP;NS=2548	GT	1/1	0/1
EOF

tabix -p vcf ./input1_chr4.vcf.gz


cat <<EOF | bcftools sort | bgzip -c > ./input2_chr4.vcf.gz
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=11032019_15h52m43s
##source=IGSRpipeline
##contig=<ID=4>
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
4	100157733	snp1	C	T	.	PASS	AC=415;AN=5096;DP=20612;AF=0.08;EAS_AF=0.04;EUR_AF=0.07;AFR_AF=0.01;AMR_AF=0.05;SAS_AF=0.26;VT=SNP;NS=2548	GT	0/1	1/1
4	101966771	snp2	T	C	.	PASS	AC=459;AN=5096;DP=20740;AF=0.09;EAS_AF=0.07;EUR_AF=0.17;AFR_AF=0.01;AMR_AF=0.17;SAS_AF=0.07;VT=SNP;NS=2548	GT	0/0	0/1
4	102814179	snp3	T	C	.	PASS	AC=2793;AN=5096;DP=20193;AF=0.55;EAS_AF=0.4;EUR_AF=0.59;AFR_AF=0.54;AMR_AF=0.58;SAS_AF=0.64;VT=SNP;NS=2548	GT	1/1	1/1
4	104355789	snp4	T	C	.	PASS	AC=763;AN=5096;DP=23195;AF=0.15;EAS_AF=0.08;EUR_AF=0.28;AFR_AF=0.06;AMR_AF=0.18;SAS_AF=0.2;VT=SNP;NS=2548	GT	0/1	0/0
4	10574522	snp5	T	C	.	PASS	AC=1180;AN=5096;DP=20409;AF=0.23;EAS_AF=0.33;EUR_AF=0.18;AFR_AF=0.24;AMR_AF=0.11;SAS_AF=0.25;VT=SNP;NS=2548	GT	1/0	1/1
4	105905360	snp6	A	G	.	PASS	AC=3254;AN=5096;DP=16605;AF=0.64;EAS_AF=0.7;EUR_AF=0.46;AFR_AF=0.88;AMR_AF=0.59;SAS_AF=0.48;VT=SNP;NS=2548	GT	1/1	0/1
4	106322897	snp7	C	T	.	PASS	AC=695;AN=5096;DP=19936;AF=0.14;EAS_AF=0;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.08;SAS_AF=0.1;VT=SNP;NS=2548	GT	0/1	1/0
4	106371793	snp11	C	A	.	PASS	AC=906;AN=5096;DP=19871;AF=0.18;EAS_AF=0.13;EUR_AF=0.16;AFR_AF=0.28;AMR_AF=0.21;SAS_AF=0.09;VT=SNP;NS=2548	GT	1/1	0/0
4	106524727	snp10	G	A	.	PASS	AC=2570;AN=5096;DP=18067;AF=0.5;EAS_AF=0.44;EUR_AF=0.37;AFR_AF=0.71;AMR_AF=0.35;SAS_AF=0.54;VT=SNP;NS=2548	GT	0/0	1/1
EOF

tabix -p vcf ./input2_chr4.vcf.gz

```


