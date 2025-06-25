process COMPRESS_OUTPUT {
    tag "${meta.id}"

    conda "bifrost-httr=0.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed90af4777d8d7086ed99d0a825f99e20e39278a75c70d3f4f7b6336edf7e210/data' :
        'community.wave.seqera.io/library/bifrost-httr:0.2.0--e8ca5c015e9a6142' }"

    cpus { 2 * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    tuple val(meta), path(all_fits_files)

    output:
    tuple val(meta), path("${prefix}.json{,.zip}")
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    touch ${prefix}.json.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bifrost-httr: \$(bifrost-httr --version | sed 's/bifrost-httr, version //')
    END_VERSIONS
    """
}
