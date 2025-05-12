/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
import groovy.json.JsonSlurper

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_bifrost_pipeline'

// Include process modules
include { PREPARE_INPUTS } from '../../modules/local/prepare_inputs/main.nf'
include { SPLIT_DATA } from '../../modules/local/split_data/main.nf'
include { CONC_RESPONSE_ANALYSIS } from '../../modules/local/conc_response_analysis/main.nf'
include { COMPRESS_OUTPUT } from '../../modules/local/compress_output/main.nf'

workflow BIFROST {
    take:
    ch_input
    ch_meta_mapper
    ch_counts
    ch_substances_cell_types
    ch_model
    n_cores

    main:
    // Step 1: Prepare input files from meta data, counts, and config
    ch_prepared_inputs = PREPARE_INPUTS(
        ch_input,
        ch_meta_mapper,
        ch_counts,
        ch_substances_cell_types
    )

    // Step 2: Process prepared inputs and probes
    ch_named_prepared_inputs = PREPARE_INPUTS.out.prepared_inputs
        .flatten()
        .map{ tuple([id: it.simpleName], it) }

    // Step 3: Split data for each input file
    SPLIT_DATA(ch_named_prepared_inputs)

    // Step 4: Prepare probes channel for concentration response analysis

    // Derive targz file(s) to pass to concentration response analysis
    // - With batching turned on, there will be one targz file per batch
    // - With batching turned off, there will be one targz file per input file

    ch_probes = SPLIT_DATA.out.probe_files  // meta, [manifests], [targzs]
        .map{meta, manifests, batch_files ->
            // Convert single objects to lists if needed
            def manifestList = manifests instanceof List ? manifests : [manifests]
            def batchList = batch_files instanceof List ? batch_files.sort{it.simpleName} : batch_files
            def batch_nums = manifestList.collect{it.simpleName.split('_batch')[1].toInteger()}
            [meta + [n_batches: manifestList.size()], batch_nums, manifestList.sort{it.simpleName}, batchList]
        }
        .transpose()
        .map{meta, batch_num, manifest, batch_file ->
            [meta + [batch_number: batch_num], manifest, batch_file]
        }
        .splitCsv(elem: 1).map{meta, probes, targzs -> [meta, probes[0], targzs]} // splitCSV always returns a list, and we know we have one probe per line
        .groupTuple().map{meta, probes, targzs -> [meta, probes, targzs[0]]} // We know that there is only one targz per batch
        .view()

    // Step 5: Run concentration response analysis
    CONC_RESPONSE_ANALYSIS(
        ch_model,
        ch_probes
    )

    // Step 6: Compress and output results

    ch_results_for_compression = CONC_RESPONSE_ANALYSIS.out.compressed_fits_files
            .map { meta, fits ->
                tuple(groupKey([id: meta.name], meta.n_batches), fits) // Use groupKey to release results as soon as all batches for a file have been processed
            }
            .groupTuple()

    COMPRESS_OUTPUT(ch_results_for_compression)

}
