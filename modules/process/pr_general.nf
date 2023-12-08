process QUERY_VCF_SAMPLES {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 2

    input:
    tuple val(chr), val(group), val(genofile), path("genodir")

    output:
    tuple val(chr),  val(group), path("${chr}_${group}_samples.txt"), emit: samples

    script:
    """
    bcftools query -l "genodir/${genofile}" | sort > "${chr}_${group}_samples.txt"


    """
}

process QUERY_VCF_VARIANTS {
    publishDir "${params.outdir}/intermediates", mode: 'copy'
 
    cpus 2

    input:
    tuple val(chr), val(group), val(split), val(genofile)

    output:
    tuple val(chr),  val(group), val(split), path("${chr}_${group}_ID"), emit: ids

    script:
    """
    bcftools query -f '%ID\n' "${genofile}" | sort > ${chr}_${group}_ID


    """
}

process SPLIT_VCF {
    publishDir "${params.outdir}/intermediates/splits", mode: 'copy'
 
    input:
    tuple val(chr), val(group), path(genofile)

    output:
    tuple val(chr),  val(group), path("*")

    script:
    """
    zcat ${genofile} | split_vcf.sh "${params.maxvariants}"
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
    tuple val(chr), val(split), val(group1), path(ID1), val(group2), path(ID2)

    output:
    tuple val(chr), path("${chr}_${split}_overlapping_ids.txt"), val(split)

    script:
    """
    comm -12 $ID1 $ID2 > ${chr}_${split}_overlapping_ids.txt
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
    publishDir "${params.outdir}/intermediates/${group}", mode: 'copy'
 
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
    publishDir "${params.outdir}/intermediates/${group}", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(split), val(group), path(genofile), path("overlapping_ids.txt")

    output:
    tuple val(chr), val(split), val(group), path("${chr}_${split}_${group}_common.vcf.gz")

    script:
    """
    bcftools view -i 'ID=@overlapping_ids.txt' -Oz "${genofile}" > ${chr}_${split}_${group}_common.vcf.gz
    """
}

process REMOVE_INFO_VCF {
    publishDir "${params.outdir}/intermediates/${group}", mode: 'copy'
 
    cpus 1

    input:
    tuple val(chr), val(split), val(group), path(genofile)

    output:
    tuple val(chr), val(group), val(split), path("${chr}_${group}_${split}_no_info.vcf.gz"), emit: no_info
    tuple val(chr), val(group), val(split), path("${chr}_${group}_${split}_info_data.txt"), emit: info

    script:
    """
    # Take now and later remove all INFO columns, as they are impossible to keep in a merged dataset. 
    bcftools view -h ${genofile} | grep "^##INFO" | cut -d ',' -f1 | cut -d '=' -f3 | awk '{printf("%%" \$1 ";")}' | sed 's/;\$//' > info_fields.txt

    # Use bcftools query to extract the data
    format_string="%ID;\$(cat info_fields.txt)"

    bcftools query -H -f "\${format_string}\n" ${genofile} > ${chr}_${group}_${split}_info_data.txt 

    # Remove INFO
    bcftools annotate --remove \$(cat info_fields.txt | sed 's/%/INFO\\//g' | tr ';' ',') -Oz -o ${chr}_${group}_${split}_no_info.vcf.gz ${genofile}

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
    publishDir "${params.outdir}/intermediates", mode: 'rellink'
 
    cpus 1

    input:
    tuple val(chr), val(split), val(group1), path("vcf1_common.vcf.gz"), val(group2), path("vcf2_common.vcf.gz")

    output:
    tuple val(chr), val(split), path("${chr}_${split}_merged.vcf.gz")

    script:
    """
    tabix -p vcf vcf1_common.vcf.gz
    tabix -p vcf vcf2_common.vcf.gz
    bcftools merge -m id -Oz -o "${chr}_${split}_merged.vcf.gz" vcf1_common.vcf.gz vcf2_common.vcf.gz
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

process CONCAT_VCF {
    publishDir "${params.outdir}/merged", mode: 'copy'
 
    input:
    tuple val(chr), val(files_order), path(files)

    output:
    tuple val(chr), path("${chr}_merged.vcf.gz"), path("${chr}_merged.vcf.gz.tbi")

    script:
    def files_order2 = files_order.join("\n")
    def files2 = files.join("\n")
    """
    echo "${files_order2}" > flo
    echo "${files2}" > fl
    paste flo fl > to_concatenate
    sort -n -k1,1 to_concatenate | cut -f2 > to_concatenate_2

    bcftools concat --naive --file-list to_concatenate_2 -Oz -o "${chr}_merged.vcf.gz"
    tabix -p vcf "${chr}_merged.vcf.gz"
    """
}

process CONCAT_INDEP_VCF {
    publishDir "${params.outdir}/${group}", mode: 'copy'
 
    input:
    tuple val(chr), val(group), val(files_order), path(files)

    output:
    tuple val(chr), val(group), path("${chr}_common_variants.vcf.gz"), path("${chr}_common_variants.vcf.gz.tbi")

    script:
    def files_order2 = files_order.join("\n")
    def files2 = files.join("\n")
    """
    echo "${files_order2}" > flo
    echo "${files2}" > fl
    paste flo fl > to_concatenate
    sort -n -k1,1 to_concatenate | cut -f2 > to_concatenate_2

    bcftools concat --naive --file-list to_concatenate_2 -Oz -o "${chr}_common_variants.vcf.gz"
    tabix -p vcf "${chr}_common_variants.vcf.gz"
    """
}

process CONCAT_INFO {
    publishDir "${params.outdir}/${group}", mode: 'copy'
 
    input:
    tuple val(chr), val(group), val(files_order), path(files)

    output:
    tuple val(chr), val(group), path("${chr}_info.gz")

    script:
    def files_order2 = files_order.join("\n")
    def files2 = files.join("\n")
    """
    echo "${files_order2}" > flo
    echo "${files2}" > fl
    paste flo fl > to_concatenate
    sort -n -k1,1 to_concatenate | cut -f2 > to_concatenate_2

    concat_info.sh "to_concatenate_2" "${chr}_info"
    gzip "${chr}_info"

    """
}

