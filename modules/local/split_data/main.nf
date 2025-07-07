process SPLIT_DATA {
    tag "${meta.id}"

    conda "bifrost-httr=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37b8df53325072038a7dce735a434031e56946aeca172a683216dbdbc197fc5d/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.3.0--6161c4cc71c68c4c' }"

    cpus { 2 * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    tuple val(meta), path(input_data)

    output:
    tuple val(meta), path("${prefix}.manifest.csv"), path("${prefix}_batch*"), emit: probe_files
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    bifrost-httr split-data \\
        --input-file $input_data \\
        --output-dir . \\
        --prefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.manifest.csv
    touch ${prefix}_batch1
    touch ${prefix}_batch2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """
}
