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

// Function to extract probe names from a JSON file
def get_probes(String data_file) {
    json_slurper = new JsonSlurper()
    def dataset = json_slurper.parse(file(data_file))
    def all_probes = dataset.probes
    return all_probes
}

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

    // Step 2: Process prepared inputs to create two channels:
    // - inputs: tuples of (name, file) for SPLIT_DATA
    // - probes: tuples of (name, probes) for later use
    ch_named_prepared_inputs = PREPARE_INPUTS.out.prepared_inputs
        .flatten()
        .multiMap{
            inputs: tuple(it.simpleName, it)
            probes: tuple(it.simpleName, get_probes(it.toString()))
        }

    // Step 3: Split data for each input file
    SPLIT_DATA(ch_named_prepared_inputs.inputs)
    SPLIT_DATA.out.all_probe_files

    // Step 4: Prepare probes channel for concentration response analysis
    // - Transpose to group probes by file
    // - Group into chunks based on number of cores
    // - Combine with split data output
    ch_probes = ch_named_prepared_inputs
        .probes
        .transpose()
        .groupTuple(size: n_cores.toInteger(), remainder: true, sort: true)
        .combine(SPLIT_DATA.out.all_probe_files, by: 0)

    // Step 5: Run concentration response analysis
    CONC_RESPONSE_ANALYSIS(
        ch_model,
        ch_probes
    )

    // Step 6: Compress and output results
    COMPRESS_OUTPUT(
        CONC_RESPONSE_ANALYSIS.out.all_fits_files.groupTuple()
    )
}
