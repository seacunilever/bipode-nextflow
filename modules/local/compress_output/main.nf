process COMPRESS_OUTPUT {
    tag "${meta.id}"

    cpus { 2 * task.attempt }
    memory { 1.GB * task.attempt }

    input:
    tuple val(meta), path(all_fits_files)

    output:
    tuple val(meta), path("${prefix}.json{,.zip}"), emit: compressed_fits_files
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
        --test_substance "${meta.test_substance}" \
        --cell_type "${meta.cell_type}" \
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
