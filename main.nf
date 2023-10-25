#!/usr/bin/env nextflow

// Read the lists
list1 = file(params.list1).readLines().collect { line ->
    def (chr, file) = line.split(/\s+/)
    return [chr: chr, file: file]
}

list2 = file(params.list2).readLines().collect { line ->
    def (chr, file) = line.split(/\s+/)
    return [chr: chr, file: file]
}

// Pair up files by chromosome
pairs = list1.findAll { map1 ->
    list2.find { map2 ->
        map1.chr == map2.chr
    }
}.collect { map1 ->
    def map2 = list2.find { it.chr == map1.chr }
    return [vcf1: map1.file, vcf2: map2.file]
}

// Create the processing channel
Channel.from(pairs).set { vcf_pairs }

process run_merge_script {
    input:
    set val(vcf1), val(vcf2) from vcf_pairs

    """
    merge_vcfs.sh ./input1.vcf.gz ./input2.vcf.gz output_merged.vcf.gz output_excluded_1.vcf.gz output_excluded_2.vcf.gz
    """
}

workflow.onError {
    println("Oops... Something went wrong!")
}

