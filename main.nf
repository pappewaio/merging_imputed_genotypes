#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "results"

include { 
  QUERY
  COMMON_IDS
  COMMON_FILTER
  MERGE
} from './modules/merge_vcfs.nf'

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
    // Extract only IDs
    QUERY(ch_mixed_input)

    // Check common IDs, by joining into same using filter
    ch_query_out_file1 = QUERY.out.filter { it[1] == "file1" }
    ch_query_out_file2 = QUERY.out.filter { it[1] == "file2" }
    ch_query_out_file1.join(ch_query_out_file2, by: 0)
    .set { ch_joined_output }
    COMMON_IDS(ch_joined_output)

    // Filter the VCFs based on common IDs
    ch_mixed_input
    .combine(COMMON_IDS.out, by:0)
    .set { ch_to_filter }
    COMMON_FILTER(ch_to_filter)

    // Do the merge
    ch_cfilter_out_file1 = COMMON_FILTER.out.filter { it[1] == "file1" }
    ch_cfilter_out_file2 = COMMON_FILTER.out.filter { it[1] == "file2" }
    ch_cfilter_out_file1.join(ch_cfilter_out_file2, by: 0)
    .set { ch_joined_output_2 }
    MERGE(ch_joined_output_2)
}


