process SPLIT_DATA {
    tag "${meta.id}"

    cpus 2
    memory 1.GB

    input:
    tuple val(meta), path(input_data)

    output:
    tuple val(meta), path("${prefix}.manifest.csv"), path("${prefix}_batch*"), emit: probe_files
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    bipode-httr split-data \\
        --input-file $input_data \\
        --output-dir . \\
        --prefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bipode-httr: \$(bipode-httr --version | sed 's/bipode-httr, version //')
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
        bipode-httr: \$(bipode-httr --version | sed 's/bipode-httr, version //')
    END_VERSIONS
    """
}
