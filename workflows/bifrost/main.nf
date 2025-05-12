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

    // Process probe files and targz files for concentration response analysis
    // The input channel contains: meta, [manifests], [targzs] where:
    // - manifests: List of manifest files (one per batch)
    // - targzs: List of targz files (one per batch when batching is on, one per input file when off)
    // The transformations:
    // 1. Convert single objects to lists if needed and sort batch files
    // 2. Extract batch numbers from manifest names
    // 3. Transpose to get per-batch processing
    // 4. Split CSV to get individual probes
    // 5. Group by meta to collect all probes for a batch
    // Note: splitCsv and groupTuple operations ensure proper list handling

    probes_by_tarname = SPLIT_DATA.out.manifest
        .splitCsv(sep: '\t', header: true)
        .map{meta, row ->
            [meta + ['tar_name': row.tar_file], meta, row.probes.split(',')]
        }

    tars_by_tarname = SPLIT_DATA.out.probe_files
        .transpose()
        .map{meta, targz ->
            [meta + ['tar_name': targz.name], targz.simpleName.split('_batch')[1].toInteger(), targz]
        }

    ch_probes = probes_by_tarname
        .combine(tars_by_tarname, by: 0).map{it.tail()}
        .groupTuple()
        .map{meta, batches, batch_nums, batch_files ->
             [meta + [n_batches: batch_nums.size()], batch_nums, batches, batch_files]
        } // meta, batch_nums, batches, batch_file(s)
        .transpose() // meta, batch_num, batch, batch_file
        .map{meta, batch_num, batch, batch_file ->
            [meta + [batch_number: batch_num], batch, batch_file]
        } // meta, batch_num, batch, batch_file

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
