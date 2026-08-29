process CREATE_MULTIQC_REPORT {

    tag "${meta.id}"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(meta), path(input_file), val(cell_type), val(test_substance)

    output:
    path "${task.ext.prefix ?: meta.id}.html", emit: report
    path "${task.ext.prefix ?: meta.id}_data", emit: data
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def report_prefix = task.ext.prefix ?: meta.id

    """
    bipode-httr create-report \
        --summary-file ${input_file} \
        --test-substance "${test_substance}" \
        --cell-type "${cell_type}" \
        --output-name "${report_prefix}" \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
        bipode-httr: \$(bipode-httr --version | sed 's/bipode-httr, version //')
    END_VERSIONS
    """

    stub:
    def report_prefix = task.ext.prefix ?: meta.id

    """
    touch "${report_prefix}.html"
    mkdir -p "${report_prefix}_data"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
        bipode-httr: \$(bipode-httr --version | sed 's/bipode-httr, version //')
    END_VERSIONS
    """
}
