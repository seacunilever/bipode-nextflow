//
// Subworkflow with functionality specific to the seacunilever/bipode-nextflow pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters        } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args // array: List of positional Nextflow CLI arguments
    outdir            // string: Output directory
    input             // string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required, and dump pipeline parameters
    // to a JSON file.
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters against the configured nextflow_schema.json.
    //
    // paramsSummaryLog is deliberately not called here because it fails
    // while processing the global process.container configuration.
    //
    if (validate_params) {
        validateParameters()
    }

    //
    // Check the configuration provided to the pipeline.
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Create the input channel from the supplied file path.
    //
    ch_input = Channel.fromPath(params.input, checkIfExists: true)

    emit:
    input    = ch_input
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir          // path: Output directory where results are published
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:

    workflow.onComplete {
        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error 'Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions,
    // for example:
    // params["run_xyz"] ? "Tool (Foo et al. 2023)" : ""

    def citation_text = [
        'Tools used in the workflow included:',
        '.'
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions,
    // for example:
    // params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : ""

    def reference_text = [
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {

    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    //
    // Pipeline DOI
    //
    if (meta.manifest_map.doi) {
        def temp_doi_ref = ''
        def manifest_doi = meta.manifest_map.doi.tokenize(',')

        manifest_doi.each { doi_ref ->
            def cleaned_doi = doi_ref
                .replace('https://doi.org/', '')
                .replace(' ', '')

            temp_doi_ref += "(doi: https://doi.org/${cleaned_doi}), "
        }

        meta['doi_text'] = temp_doi_ref.substring(
            0,
            temp_doi_ref.length() - 2
        )
    } else {
        meta['doi_text'] = ''
    }

    meta['nodoi_text'] = meta.manifest_map.doi
        ? ''
        : '<li>If available, update the text to include the Zenodo DOI of the pipeline version used.</li>'

    //
    // Tool references
    //
    meta['tool_citations'] = ''
    meta['tool_bibliography'] = ''

    // TODO nf-core: Uncomment when toolCitationText and toolBibliographyText
    // contain pipeline-specific citation information.
    //
    // meta['tool_citations'] = toolCitationText()
    //     .replaceAll(', \\.', '.')
    //     .replaceAll('\\. \\.', '.')
    //
    // meta['tool_bibliography'] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text
    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}