process SPLIT_DATA {
    tag "${meta.id}"

    conda "bifrost-httr=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05fa08012fb11ccb282c05e1b48b53c6220b0853692cafcbaf8829749d6aabc/data' :
        'wave.seqera.io/wt/25fa77f460cd/wave/build:bifrost-httr-0.1.0--2c648d2de87966a9' }"

    cpus { 2 * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    tuple val(meta), path(input_data)

    output:
    tuple val(meta), path("${prefix}.manifest.csv"), path("${prefix}_batch*"), emit: probe_files

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    bifrost-httr split-data \\
        --input-file $input_data \\
        --analysis-dir . \\
        --prefix $prefix \\
        $args
    """
}
