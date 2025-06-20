process COMPRESS_OUTPUT {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33499ba0fef01239be35b4b8ffae2f35bc921bd88d622a5f5f5c6ed2edb3eaa0/data' :
        'community.wave.seqera.io/library/python_cmdstanpy_numpy_pandas_pruned:b21b7854a692918a' }"

    input:
    tuple val(meta), path(all_fits_files)

    output:
    tuple val(meta), path("${prefix}.json{,.zip}")

    script:
    def args = task.ext.args ?: ''
    def fits_files = all_fits_files.join(" ")
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_ext = args.contains('--no-compression') ? '.json' : '.json.zip'
    """
    mkdir Fits
    for file in $fits_files; do tar -zxf "\$file" -C Fits/; done

    bifrost-httr compress-output \
        --fits-dir Fits \
        --output ${prefix}${output_ext} \
        $args
    """
}
