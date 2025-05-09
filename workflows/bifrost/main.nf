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
    // - Transpose to group probes by file
    // - Group into chunks based on number of cores
    // - Count the number of chunks per file (for later use with groupKey)
    // - Combine with split data output
    ch_probes = ch_named_prepared_inputs
        .splitJson(path: 'probes')                                          // [name, [probe]]
        .map{ meta, probe -> tuple([id: meta.id], probe) }                  // [name, [batch of probes]]
        .groupTuple(size: n_cores.toInteger(), remainder: true, sort: true) // [name, [[batch of probes], [batch of probes], ...]]
        .groupTuple()
        .map{meta, batches ->
            [meta, batches.size(), batches, (1..batches.size())]            // [name, n_batches, [batches], [batch_numbers]]
        }
        .transpose()                                                        // [name, n_batches, batch, batch_number]
        .combine(SPLIT_DATA.out.all_probe_files, by: 0)                     // [name, n_batches, batch, batch_number, probe_file]
        .map{ meta, n_batches, batch, batch_number, probe_file ->
            [[id: meta.id + "_" + batch_number, name: meta.id, n_batches: n_batches, batch_number: batch_number], batch, probe_file]
        }

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
