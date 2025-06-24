process COMPRESS_OUTPUT {
    tag "${meta.id}"

    conda "bifrost-httr=0.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05fa08012fb11ccb282c05e1b48b53c6220b0853692cafcbaf8829749d6aabc/data' :
        'wave.seqera.io/wt/25fa77f460cd/wave/build:bifrost-httr-0.1.0--2c648d2de87966a9' }"

    cpus { 2 * task.attempt }
    memory { 1.GB * task.attempt }

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
