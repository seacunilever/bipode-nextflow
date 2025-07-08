process CREATE_MULTIQC_REPORT {
    tag "${meta.id}"
    publishDir "${params.outdir}/reports", mode: 'copy'

    conda "bifrost-httr=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/76/76e8817651482fe89237efe5d385050d40144519c9f0c9fc5b0f9ee506292428/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.3.1--b4c49de956618921' }"

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
