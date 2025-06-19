process SPLIT_DATA {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

    input:
    tuple val(meta), path(input_data)

    output:
    tuple val(meta), path("${prefix}.manifest.csv"), path("${prefix}_batch*"), emit: probe_files

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    split_data.py \\
        --input-file $input_data \\
        --analysis-dir . \\
        --prefix $prefix \\
        $args
    """
}
