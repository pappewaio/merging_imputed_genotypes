#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {  
 QUERY_VCF_VARIANTS
 QUERY_VCF_SAMPLES
 SPLIT_VCF
 COMMON_IDS
 COMMON_SAMPLES
 COMMON_SAMPLES_FILTER_VCF
 COMMON_ID_FILTER_VCF
 REMOVE_INFO_VCF
 MERGE_VCF
 CONCAT_VCF
 CONCAT_INFO
 CONCAT_INDEP_VCF
} from '../process/pr_general.nf'

workflow merge_vcfs {

    take:
    ch_mixed_input

    main:

    // Query only for sample IDs
    QUERY_VCF_SAMPLES(ch_mixed_input)

    // Check common Samples, by joining into same using filter
    ch_query_out_file1 = QUERY_VCF_SAMPLES.out.samples.filter { it[1] == "file1" }
    ch_query_out_file2 = QUERY_VCF_SAMPLES.out.samples.filter { it[1] == "file2" }
    ch_query_out_file1.join(ch_query_out_file2, by: 0)
    .set { ch_joined_output_samples }
    COMMON_SAMPLES(ch_joined_output_samples)

    // Remove the samples which are in both sets
    ch_mixed_input
    .combine(COMMON_SAMPLES.out, by:0)
    .set { ch_to_filter_samples }
    COMMON_SAMPLES_FILTER_VCF(ch_to_filter_samples)

    // chunk the vcf for less waiting time
    SPLIT_VCF(COMMON_SAMPLES_FILTER_VCF.out)
    SPLIT_VCF.out
    .flatMap { chr, group, files ->
      files.collect { file ->
        def splitMatcher = file =~ /split_(\d+)\.vcf$/
        if (splitMatcher) {
          def splitNumber = splitMatcher[0][1]
          return tuple(chr, group, splitNumber, file)
        }
      }
    }
    .set { ch_split_vcf }


    // Query only for variant IDs
    QUERY_VCF_VARIANTS(ch_split_vcf)

    // Check common IDs, by joining into same using filter
    ch_query_out_file1 = QUERY_VCF_VARIANTS.out.ids.filter { it[1] == "file1" }
    ch_query_out_file2 = QUERY_VCF_VARIANTS.out.ids.filter { it[1] == "file2" }
    ch_query_out_file1.join(ch_query_out_file2, by: [0,2])
    .set { ch_joined_output_ids }
    COMMON_IDS(ch_joined_output_ids)

    // Filter the VCFs based on common IDs
  //  ch_split_vcf.view()
  //  COMMON_IDS.out.view()

    ch_split_vcf
    .combine(COMMON_IDS.out, by:[0,2])
    .set { ch_to_filter }
    COMMON_ID_FILTER_VCF(ch_to_filter)

  //  // Remove INFO as it will be misleading in the merged file
    REMOVE_INFO_VCF(COMMON_ID_FILTER_VCF.out)

  //  // Do the merge
    ch_cfilter_out_file1 = REMOVE_INFO_VCF.out.no_info.filter { it[1] == "file1" }
    ch_cfilter_out_file2 = REMOVE_INFO_VCF.out.no_info.filter { it[1] == "file2" }
    ch_cfilter_out_file1.join(ch_cfilter_out_file2, by: [0,2])
    .set { ch_joined_output_2 }
    MERGE_VCF(ch_joined_output_2)
   // MERGE_VCF.out.view()

   // Concatenate per chromosome
   MERGE_VCF.out
   .groupTuple(by: 0)
   .set {ch_to_concatenate}
   CONCAT_VCF(ch_to_concatenate)

   // Concatenate per chromosome (indep vcfs)
   REMOVE_INFO_VCF.out.no_info
   .groupTuple(by: [0,1])
   .set {ch_indep_to_concatenate}
   CONCAT_INDEP_VCF(ch_indep_to_concatenate)

   // Concatenate also info
   REMOVE_INFO_VCF.out.info
   .groupTuple(by: [0,1])
   .set {ch_info_to_concatenate}
   CONCAT_INFO(ch_info_to_concatenate)
   //ch_info_to_concatenate.view()
}


