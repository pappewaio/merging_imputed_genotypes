#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "results"

include { MERGE_VCFS } from './modules/merge_vcfs.nf'

// Read the paths from the configuration
set1_path = params.dir1
set2_path = params.dir2

// Read the files into channels
channel1 = Channel.fromPath(params.set1).splitCsv(sep:'\t', header:false).map { row -> tuple(row[0], file("${set1_path}/${row[1]}")) }
channel2 = Channel.fromPath(params.set2).splitCsv(sep:'\t', header:false).map { row -> tuple(row[0], file("${set2_path}/${row[1]}")) }

// Join channels based on chromosome
pairs_channel = channel1.join(channel2, by: 0)

workflow {
    MERGE_VCFS(pairs_channel)
}

