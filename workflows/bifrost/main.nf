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
    // The input channels contain:
    // - SPLIT_DATA.out.manifest: meta, manifest (tab-separated file with tar_file and probes columns)
    // - SPLIT_DATA.out.probe_files: meta, [probe_files] (list of targz files)
    // The transformations:
    // 1. For manifests: Split CSV to get tar_file and probes, then create [meta_tar_name, meta, probes] tuples
    // 2. For probe files: Transpose and extract batch numbers to create [meta_tar_name, batch_num, targz] tuples
    // 3. Combine both channels by tar_name to match probes with their corresponding targz files
    // 4. Group by meta to collect all batches and their files
    // 5. Add batch count to meta and transpose for per-batch processing
    // 6. Final output: [meta + batch_number, batch, batch_file] for each batch
    // Note: The tar_name field is used as a key to match probes with their corresponding targz files

    probes_by_tarname = SPLIT_DATA.out.manifest // meta, manifest
        .splitCsv(sep: '\t', header: true) // meta, [tar_file, probes]
        .map{meta, row ->
            [meta + ['tar_name': row.tar_file], meta, row.probes.split(',')]
        } // [meta_tar_name], meta, probes]

    tars_by_tarname = SPLIT_DATA.out.probe_files //[ meta, [probe_files]]
        .transpose() // [meta, probe_file]
        .map{meta, targz ->
            [meta + ['tar_name': targz.name], targz.simpleName.split('_batch')[1].toInteger(), targz]
        } // [meta_tar_name, batch_num, targz]

    ch_probes = probes_by_tarname
        .combine(tars_by_tarname, by: 0).map{it.tail()} // [meta, probes, batch_num, targz]
        .groupTuple() // [meta, [batch_nums], [batches], [batch_files]]
        .map{meta, batches, batch_nums, batch_files ->
             [meta + [n_batches: batch_nums.size()], batch_nums, batches, batch_files]
        } // meta, [batch_nums], [batches], [batch_files]
        .transpose() // meta, batch_num, batch, batch_file
        .map{meta, batch_num, batch, batch_file ->
            [meta + [batch_number: batch_num], batch, batch_file]
        } // meta, batch, batch_file

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
