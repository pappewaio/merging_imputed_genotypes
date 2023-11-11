#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "results"

include { merge_vcfs } from './modules/subworkflow/wf_merge_vcf.nf'
include { merge_pgens } from './modules/subworkflow/wf_merge_pgen.nf'

// Initialize variables
def set1_path, set2_path, set1, set2

if (params.type == "vcf") {
    set1_path = params.mergevcf.dir1
    set2_path = params.mergevcf.dir2
    set1 = params.mergevcf.set1
    set2 = params.mergevcf.set2
} else if (params.type == "pgen") {
    set1_path = params.mergepgen.dir1
    set2_path = params.mergepgen.dir2
    set1 = params.mergepgen.set1
    set2 = params.mergepgen.set2
} else {
    // Handle unexpected type or set a default
    error "Unknown params.type value: ${params.type}"
}

// Read the files into channels
channel1 = Channel
    .fromPath(set1)
    .splitCsv(sep:'\t', header:false)
    .map { row -> tuple(row[0], "file1", "${row[1]}", file("${set1_path}")) }

channel2 = Channel
    .fromPath(set2)
    .splitCsv(sep:'\t', header:false)
    .map { row -> tuple(row[0], "file2", "${row[1]}", file("${set2_path}")) }

// Mix input channels to allow for parallell processing
channel1.mix(channel2)
.set { ch_mixed_input }

workflow {
  if("${params.type}"=="vcf"){merge_vcfs(ch_mixed_input)}
  if("${params.type}"=="pgen"){merge_pgens(ch_mixed_input)}
}

