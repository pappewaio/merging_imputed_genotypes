process MERGE_VCFS {
    publishDir "${params.outdir}/merged_files", mode: 'copy'

    input:
    tuple val(chr), path(vcf1), path(vcf2)

    output:
    tuple val(chr), path("${chr}_output_merged.vcf.gz"), path("${chr}_output_excluded_1.vcf.gz"), path("${chr}_output_excluded_2.vcf.gz")

    script:
    """
    merge_vcfs.sh $vcf1 $vcf2 ${chr}_output_merged.vcf.gz ${chr}_output_excluded_1.vcf.gz ${chr}_output_excluded_2.vcf.gz
    """
}

