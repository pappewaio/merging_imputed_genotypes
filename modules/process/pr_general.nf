process QUERY_VCF {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 2

    input:
    tuple val(chr), val(group), val(genofile), path("genodir")

    output:
    tuple val(chr),  val(group), path("${chr}_${group}_ID"), emit: ids
    tuple val(chr),  val(group), path("${chr}_${group}_samples.txt"), emit: samples

    script:
    """
    bcftools query -f '%ID\n' "genodir/${genofile}" | sort > ${chr}_${group}_ID
    bcftools query -l "genodir/${genofile}" > "${chr}_${group}_samples.txt"

    """
}

process QUERY_PGEN {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 2

    input:
    tuple val(chr), val(group), val(genofile), path("genodir")

    output:
    tuple val(chr),  val(group), path("${chr}_${group}_ID")

    script:
    """
    plink2 --pfile "genodir/${genofile}" --write-snplist --out extracted_rsids 
    sort extracted_rsids.snplist > ${chr}_${group}_ID
    """
}

process COMMON_IDS {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(group1), path(ID1), val(group2), path(ID2)

    output:
    tuple val(chr), path("${chr}_overlapping_ids.txt")

    script:
    """
    comm -12 $ID1 $ID2 > ${chr}_overlapping_ids.txt
    """
}

process COMMON_SAMPLES {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(group1), path(ID1), val(group2), path(ID2)

    output:
    tuple val(chr), path("${chr}_overlapping_samples.txt")

    script:
    """
    comm -12 $ID1 $ID2 > ${chr}_overlapping_samples.txt
    """
}


process COMMON_SAMPLES_FILTER_VCF {
    publishDir "${params.outdir}/${group}", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(group), val(genofile), path("genodir"), path("overlapping_samples.txt")

    output:
    tuple val(chr), val(group), path("${chr}_${group}_independent_samples.vcf.gz")

    script:
    """
    bcftools view -S '^overlapping_samples.txt' -Oz "genodir/${genofile}" > ${chr}_${group}_independent_samples.vcf.gz
    """
}

process COMMON_ID_FILTER_VCF {
    publishDir "${params.outdir}/${group}", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(group), path(genofile), path("overlapping_ids.txt")

    output:
    tuple val(chr), val(group), path("${chr}_${group}_common.vcf.gz")

    script:
    """
    bcftools view -i 'ID=@overlapping_ids.txt' -Oz "${genofile}" > ${chr}_${group}_common.vcf.gz
    """
}

process COMMON_FILTER_PGEN {
    publishDir "${params.outdir}/${group}", mode: 'copy'

    cpus 1

    input:
    tuple val(chr), val(group), val(genofile), path("genodir"), path("overlapping_ids.txt")

    output:
    tuple val(chr), val(group), path("${chr}_${group}_common.pgen"), path("${chr}_${group}_common.pvar"), path("${chr}_${group}_common.psam")

    script:
    """
    plink2 --pfile genodir/${genofile} \
           --extract overlapping_ids.txt \
           --make-pgen \
           --sort-vars \
           --out ${chr}_${group}_common
    """
}

process MERGE_VCF {
    publishDir "${params.outdir}", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(group1), path("vcf1_common.vcf.gz"), val(group2), path("vcf2_common.vcf.gz")

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

process MERGE_PGEN {
    publishDir "${params.outdir}", mode: 'copy'

    cpus 1

    input:
    tuple val(chr), val(group1), path("pgen1_common.pgen"), path("pgen1_common.pvar"), path("pgen1_common.psam"), val(group2), path("pgen2_common.pgen"), path("pgen2_common.pvar"), path("pgen2_common.psam")

    output:
    tuple val(chr), path("${chr}_merged.pgen"), path("${chr}_merged.pvar"), path("${chr}_merged.psam")

    script:
    """
    # Create a pmerge_list.txt file containing paths to the files to be merged
    echo "pgen2_common.pgen pgen2_common.pvar pgen2_common.psam" >> pmerge_list.txt

    # Run PLINK2 to merge the PGEN files
    plink2 --pfile "pgen1_common" \
           --pmerge-list pmerge_list.txt \
           --make-pgen \
           --out ${chr}_merged
    """


}





