process QUERY_VCF {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 2

    input:
    tuple val(chr), val(file), path(vcf)

    output:
    tuple val(chr),  val(file), path("${chr}_${file}_ID")

    script:
    """
    bcftools query -f '%ID\n' "$vcf" | sort > ${chr}_${file}_ID
    """
}

process COMMON_IDS {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(file1), path(ID1), val(file2), path(ID2)

    output:
    tuple val(chr), path("${chr}_overlapping_ids.txt")

    script:
    """
    comm -12 $ID1 $ID2 > ${chr}_overlapping_ids.txt
    """
}

process COMMON_FILTER {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(file), path(vcf), path("overlapping_ids.txt")

    output:
    tuple val(chr), val(file), path("${chr}_${file}_common.vcf.gz")

    script:
    """
    bcftools view -i 'ID=@overlapping_ids.txt' -Oz "$vcf" > ${chr}_${file}_common.vcf.gz
    """
}

process MERGE_VCF {
    publishDir "${params.outdir}", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(file1), path("vcf1_common.vcf.gz"), val(file2), path("vcf2_common.vcf.gz")

    output:
    tuple val(chr), path("${chr}_merged.vcf.gz"), path("${chr}_merged.vcf.gz.tbi")

    script:
    """
    tabix -p vcf vcf1_common.vcf.gz
    tabix -p vcf vcf2_common.vcf.gz
    bcftools merge -m id -Oz -o "${chr}_merged.vcf.gz" vcf1_common.vcf.gz vcf2_common.vcf.gz
    tabix -p vcf "${chr}_merged.vcf.gz"
    """
}


