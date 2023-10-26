process MERGE_VCFS {
    publishDir "${params.outdir}/merged", mode: 'copy', pattern: '*merged*'
    publishDir "${params.outdir}/excluded", mode: 'copy', pattern: '*excluded*'
    publishDir "${params.outdir}/extra", mode: 'copy', pattern: '*ids*'

    input:
    tuple val(chr), path(vcf1), path(vcf2)

    output:
    tuple val(chr), 
      path("${chr}_output_merged.vcf.gz"), 
      path("${chr}_output_excluded_1.vcf.gz"), 
      path("${chr}_output_excluded_2.vcf.gz"),
      path("${chr}_output_merged.vcf.gz.tbi"), 
      path("${chr}_output_excluded_1.vcf.gz.tbi"), 
      path("${chr}_output_excluded_2.vcf.gz.tbi")
    output: path('*ids*')

    script:
    """
    merge_vcfs.sh \
      $vcf1 \
      $vcf2 \
      ${chr}_output_merged.vcf.gz \
      ${chr}_output_excluded_1.vcf.gz \
      ${chr}_output_excluded_2.vcf.gz \
      ${chr}_overlapping_ids.txt \
      ${chr}_vcf1_unique_ids.txt \
      ${chr}_vcf2_unique_ids.txt
    """
}

