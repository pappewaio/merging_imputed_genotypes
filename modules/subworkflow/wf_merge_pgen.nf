#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 QUERY_PGEN
 COMMON_IDS
COMMON_FILTER_PGEN
 MERGE_PGEN
} from '../process/pr_general.nf'

workflow merge_pgens {

    take:
    ch_mixed_input

    main:

    // Extract only IDs
    QUERY_PGEN(ch_mixed_input)

    // Check common IDs, by joining into same using filter
    ch_query_out_file1 = QUERY_PGEN.out.filter { it[1] == "file1" }
    ch_query_out_file2 = QUERY_PGEN.out.filter { it[1] == "file2" }
    ch_query_out_file1.join(ch_query_out_file2, by: 0)
    .set { ch_joined_output }
    COMMON_IDS(ch_joined_output)

    // Filter the VCFs based on common IDs
    ch_mixed_input
    .combine(COMMON_IDS.out, by:0)
    .set { ch_to_filter }
    COMMON_FILTER_PGEN(ch_to_filter)

    // Do the merge
    ch_cfilter_out_file1 = COMMON_FILTER_PGEN.out.filter { it[1] == "file1" }
    ch_cfilter_out_file2 = COMMON_FILTER_PGEN.out.filter { it[1] == "file2" }
    ch_cfilter_out_file1.join(ch_cfilter_out_file2, by: 0)
    .set { ch_joined_output_2 }
    MERGE_PGEN(ch_joined_output_2)
    
}

