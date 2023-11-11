#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 QUERY_VCF
 COMMON_IDS
 COMMON_FILTER
 MERGE_VCF
} from '../process/pr_general.nf'

workflow merge_vcfs {

    take:
    ch_mixed_input

    main:

    // Extract only IDs
    QUERY_VCF(ch_mixed_input)

    // Check common IDs, by joining into same using filter
    ch_query_out_file1 = QUERY_VCF.out.filter { it[1] == "file1" }
    ch_query_out_file2 = QUERY_VCF.out.filter { it[1] == "file2" }
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
    MERGE_VCF(ch_joined_output_2)
    
}

