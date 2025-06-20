process CREATE_MULTIQC_REPORT {
    tag "${meta.id}"
    publishDir "${params.outdir}/reports", mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/89/89f6935902ad9f958380b66f7bb6b1e1bbc590907dff9438da89f61d3441bd82/data' :
        'community.wave.seqera.io/library/multiqc_numpy_pandas_python:821220a16cc579d0' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input_file), val(cell_type), val(test_substance)

    output:
    path "${prefix}.html", emit: report
    path "${prefix}_data", emit: data
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bifrost-httr create-report \\
        --summary-file ${input_file} \\
        --test-substance "${test_substance}" \\
        --cell-type "${cell_type}" \\
        --output-name "${prefix}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_multiqc_report.html"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
