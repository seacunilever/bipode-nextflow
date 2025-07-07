process CREATE_MULTIQC_REPORT {
    tag "${meta.id}"
    publishDir "${params.outdir}/reports", mode: 'copy'

    conda "bifrost-httr=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37b8df53325072038a7dce735a434031e56946aeca172a683216dbdbc197fc5d/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.3.0--6161c4cc71c68c4c' }"

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
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_multiqc_report.html"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """
}
