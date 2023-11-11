#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "results"

include { merge_vcfs } from './modules/subworkflow/wf_merge_vcf.nf'
include { merge_pgens } from './modules/subworkflow/wf_merge_pgen.nf'

// Read the paths from the configuration
set1_path = params.dir1
set2_path = params.dir2

// Read the files into channels
channel1 = Channel.fromPath(params.set1).splitCsv(sep:'\t', header:false).map { row -> tuple(row[0],"file1", file("${set1_path}/${row[1]}")) }
channel2 = Channel.fromPath(params.set2).splitCsv(sep:'\t', header:false).map { row -> tuple(row[0],"file2", file("${set2_path}/${row[1]}")) }

// Mix input channels to allow for parallell processing
channel1.mix(channel2)
.set { ch_mixed_input }

workflow {
  if("${params.type}"=="vcf"){merge_vcfs(ch_mixed_input)}
  if("${params.type}"=="pgen"){merge_pgens(ch_mixed_input)}
}

