/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_bipode_pipeline'

// Include process modules
include { PREPARE_INPUTS } from '../../modules/local/prepare_inputs/main.nf'
include { SPLIT_DATA } from '../../modules/local/split_data/main.nf'
include { COMPILE_STAN_MODEL } from '../../modules/local/compile_stan_model/main.nf'
include { CONC_RESPONSE_ANALYSIS } from '../../modules/local/conc_response_analysis/main.nf'
include { COMPRESS_OUTPUT } from '../../modules/local/compress_output/main.nf'
include { CREATE_MULTIQC_REPORT } from '../../modules/local/create_multiqc_report/main.nf'

workflow BIPODE {
    take:
    ch_input
    ch_meta_mapper
    ch_counts
    ch_bipode_config
    ch_model
    n_cores
    precompile_model

    main:

    // Branch inputs into JSON and raw data
    ch_input_json = ch_input.branch { input_file ->
        json: input_file.name.endsWith('.json')
        raw: true
    }

    // (optional) Step 1: Prepare input files from meta data, counts, and config
    PREPARE_INPUTS(
        ch_input_json.raw,
        ch_meta_mapper,
        ch_counts,
        ch_bipode_config
    )

    // Step 2: Process prepared inputs and probes
    ch_named_prepared_inputs = PREPARE_INPUTS.out.prepared_inputs
        .flatten()
        .mix(ch_input_json.json)
        .map { input_file ->

            def json = new groovy.json.JsonSlurper().parseText(input_file.text)

            [
                [
                    id             : input_file.simpleName,
                    test_substance : json.test_substance,
                    cell_type      : json.cell_type
                ],
                input_file
            ]
        }

    // Step 3: Split data for each input file
    SPLIT_DATA(ch_named_prepared_inputs)

    // Step 4: Compile Stan model once (if precompile_model is true)
    if (precompile_model) {
        COMPILE_STAN_MODEL(ch_model)
        ch_model_for_analysis = COMPILE_STAN_MODEL.out.compiled_model
    }
    else {
        ch_model_for_analysis = ch_model
    }

    // Step 5: Run concentration response analysis using pre-compiled model

    // SPLIT_DATA can package probes into tar.gz files in different ways. It
    // produces a manifest to describe which probes are in which tar.gz. To
    // prepare an input channel for the concentration response analysis, we
    // just need to parse out the probes from the manifest, and count the
    // number of batches per file for the later benefit of groupKey.

    ch_probes = SPLIT_DATA.out.probe_files
        .splitCsv(sep: '\t', header: true, elem: 1)
        .groupTuple()
        .map { meta, rows, targzs ->

            def tar_to_gz = targzs
                .flatten()
                .collectEntries { tar_file ->
                    [(tar_file.name): tar_file]
                }

            [
                meta + [n_batches: rows.size()],
                rows*.batch,
                rows*.probes*.split(','),
                rows*.tar_file.collect { tar_file_name ->
                    tar_to_gz[tar_file_name]
                }
            ]
        }
        .transpose()
        .map { meta, batch_num, batch, batch_file ->

            [
                meta + [
                    id           : "${meta.id}_${batch_num}",
                    name         : meta.id,
                    batch_number : batch_num
                ],
                batch,
                batch_file
            ]
        }

    CONC_RESPONSE_ANALYSIS(
        ch_model_for_analysis,
        ch_probes
    )

    // Step 6: Compress and output results

    ch_results_for_compression = CONC_RESPONSE_ANALYSIS.out.compressed_fits_files
        .map { meta, fits ->

            def grouped_meta = meta.findAll { key, _value ->
                key != 'batch_number'
            }

            tuple(
                groupKey(
                    grouped_meta + [id: meta.name],
                    meta.n_batches
                ),
                fits
            )
        }
        .groupTuple()

    COMPRESS_OUTPUT(ch_results_for_compression)

    // Step 7: Create reports

    ch_compressed_output = COMPRESS_OUTPUT.out.compressed_fits_files
        .map { meta, fits ->
            [
                meta,
                fits,
                meta.cell_type,
                meta.test_substance
            ]
        }

    CREATE_MULTIQC_REPORT(ch_compressed_output)
}